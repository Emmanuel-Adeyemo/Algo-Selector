library(shiny)
library(dplyr)
library(GA)
library(ggplot2)
library(GGally)
library(patchwork)
library(bslib)
library(purrr) 
library(janitor)

# --- Global Data and GA Parameters ---
# pop_dta_raw = read.csv("pop_dta.csv", stringsAsFactors = FALSE)
pop_dta_raw = readRDS('wheat_dta.RDS') %>% clean_names() %>% rename(ht = height)
coa = readRDS('modified_COA.RDS')

secondary_traits = c('fhb', 'don', 'ht')
all_traits = c('yield', secondary_traits)
overuse_weight = 0.8 
pop_size = 70
yield_weight = 2.5

# Helper for handling NULL in bounds
`%||%` = function(x, y) if (is.null(x) || is.na(x)) y else x # Added is.na(x) check


pop_dta_initial_processed = pop_dta_raw
trait_cols_numeric = pop_dta_initial_processed %>% dplyr::select(where(is.numeric)) %>% names()

for(trait in trait_cols_numeric){
  new_col = paste('Norm_', trait, sep="")
  pop_dta_initial_processed[[new_col]] = (pop_dta_initial_processed[[trait]] - mean(pop_dta_initial_processed[[trait]], na.rm=TRUE)) / sd(pop_dta_initial_processed[[trait]], na.rm=TRUE)
}



# --- UI ---
ui = page_sidebar(
  title = "Selection Optimizer",
  
  sidebar = sidebar(
    width = 350,
    
    selectInput("primary_trait", "Primary Trait", choices = trait_cols_numeric, selected = "yield"),
    
    uiOutput('primary_control'),
 
    fluidRow(
      column(4, numericInput("n_sel", "Number to Select", min = 10, max = 10000, value = 100)),
      # column(4, numericInput("pop_size", "Population Size", value = 60)),
      column(4, numericInput("parent_limit", "Max Parent Usage", value = 5)),
      column(4, numericInput("max_iter", "Max Iterations", value = 120))
    ),
    
    
    sliderInput("div_target", "Diversity Target", min = 0, max = 1, value = 0.70),
    
    uiOutput("trait_controls"),
    #hr(),
    actionButton("submit", "Run Optimizer", class = "btn-primary")
  ),
  
  tabsetPanel(
    type = 'pill',
    
    tabPanel("Trait Preview",
             br(),
             plotOutput('traitPrimary'),
             plotOutput("traitPreview", height = "600px"),
             br()
    ),
    
    tabPanel("GA Diagnostics",
             
             br(),
             fluidRow(
               column(6, verbatimTextOutput("log")), 
               column(6, plotOutput("fitnessPlot"), height='800px')  
             ),
             br(),
             fluidRow(
               column(6, plotOutput("scatterPlot")), 
               column(6, plotOutput("divYieldPlot")) 
             ),
             br()
    ),
    
    tabPanel("Selection Summary",
             plotOutput("selection_differential_plot", height='600px'),
             br(),
             fluidRow(
               column(12, downloadButton("downloadSummary", "Download Selected Pop")), 
               br(),
               column(12, DT::DTOutput("summaryTable"))
             ),
             br()
    )
  ),
  
  theme = bs_theme(
    version = 5,
    bootswatch = "yeti",
    primary = "#008c3d",
    bg = "white",
    fg = "black",
    base_font = font_google("Open Sans")
  )
)

# --- Server Logic ---
server = function(input, output, session) {
  
  # UI for secondary trait controls
  output$trait_controls = renderUI({
    
    trait_inputs = lapply(secondary_traits, function(trait) {
      tagList(
        h5(paste("Bounds for:", toupper(trait))),
        fluidRow(
          
          column(4, numericInput(paste0(trait, "_lower_prune"), HTML("Lower<br>Prune"), value = NA)),
          column(4, numericInput(paste0(trait, "_lower_thresh"), "Lower Threshold", value = NA)),
          column(4, numericInput(paste0(trait, '_lower_thres_penalty'), 'Lower Thres Penalty', value=NA))
        ),
        fluidRow(
          
          column(4, numericInput(paste0(trait, "_upper_prune"), HTML("Upper<br>Prune"), value = NA)),
          column(4, numericInput(paste0(trait, "_upper_thresh"), "Upper Threshold", value = NA)),
          column(4, numericInput(paste0(trait, '_upper_thres_penalty'), 'Upper Thres Penalty', value=NA))
        ),
        tags$hr()
      )
    })
    do.call(tagList, trait_inputs)
  })
  
  # UI for Primary trait. Added later.
  output$primary_control = renderUI({
    
    prim_trait = input$primary_trait
    
    primary_input = tagList(
      h5(paste0('Prunes for ', toupper(prim_trait))),
      
      fluidRow(
        column(6, numericInput(paste0(prim_trait, '_lower_prune'), HTML('Lower Prune'), value = NA)),
        column(6, numericInput(paste0(prim_trait, '_upper_prune'), HTML('Upper Prune'), value = NA))
      ),
      tags$hr()
    )
  
  })
  
    
  
  
  
  # --- Helper for calculating individual cross diversity ---
  calculate_individual_cross_diversity = function(data_subset, coa_matrix) {
    if (nrow(data_subset) == 0) return(numeric(0))
    
    scores = map_dbl(1:nrow(data_subset), function(i) {
      p1 = data_subset$parent1[i]
      p2 = data_subset$parent2[i]
      parents = na.omit(unique(c(p1, p2)))
      
      
      if (length(parents) == 0 || any(!parents %in% rownames(coa_matrix))) {
        return(NA_real_)
      }
      
      if (length(parents) == 1) {
        # parent with itself = 1
        return(0.8 - coa_matrix[parents, parents]) # should be 1
      } else {
        
        return(0.8 - coa_matrix[parents[1], parents[2]]) # should be 1
      }
    })
    return(scores)
  }
  
  
  # Original Population Individual Cross Diversity
  original_pop_individual_diversity_scores = reactive({
    req(pop_dta_raw) 
    
    scores = calculate_individual_cross_diversity(pop_dta_raw, coa)
    data.frame(Diversity = scores) %>% filter(!is.na(Diversity))
  })
  
  
  # for yield only
  output$traitPrimary = renderPlot({
    req(pop_dta_raw)
    
    prim_trait = input$primary_trait
    
    vals = pop_dta_raw[[prim_trait]]
    if(is.null(vals) || all(is.na(vals))){
      return(ggplot() + labs(title = paste('No data for trait: ', prim_trait)))
    }
    
    lower_prune = input[[paste0(prim_trait, '_lower_prune')]] %||% -Inf
    upper_prune = input[[paste0(prim_trait, '_upper_prune')]] %||% Inf
    
    status = case_when(
      vals < lower_prune | vals > upper_prune ~ 'Pruned',
      vals >= lower_prune & vals <= upper_prune ~ 'Safe'
    )
    status[is.na(status)] = 'Unknown'
    
    dta = data.frame(trait_value = vals, 
                     status = factor(status, levels = c('Pruned', 'Safe', 'Unknown')))
    
    ggplot(dta, aes(x = trait_value, fill = status)) +
      geom_histogram(bins = 40, alpha = 0.8, color = "white") +
      scale_fill_manual(values = c(
        Pruned = "red", Safe = "forestgreen", Unknown = "lightgrey"
      )) +
      labs(title = paste("Primary Trait:", toupper(prim_trait)), x = toupper(prim_trait), y = "Count") +
      theme_minimal() + theme(legend.position = "top")
    
  })
  
  
  output$traitPreview = renderPlot({
    req(pop_dta_raw)
    
    plots = lapply(secondary_traits, function(trait) {
      vals = pop_dta_raw[[trait]]
      if (is.null(vals) || all(is.na(vals))) {
        return(ggplot() + labs(title = paste("No data for trait:", trait)))
      }
      
      lower_thresh = input[[paste0(trait, "_lower_thresh")]] %||% -Inf
      lower_prune  = input[[paste0(trait, "_lower_prune")]]  %||% -Inf
      upper_thresh = input[[paste0(trait, "_upper_thresh")]] %||% Inf
      upper_prune = input[[paste0(trait, "_upper_prune")]]  %||% Inf
      
      status = case_when(
        vals < lower_prune | vals > upper_prune ~ "Pruned",
        (vals >= lower_prune & vals < lower_thresh) | (vals <= upper_prune & vals > upper_thresh) ~ "Penalized",
        vals >= lower_thresh & vals <= upper_thresh ~ "Safe",
        TRUE ~ "Non_Penalized"
      )
      status[is.na(status)] = "Unknown"
      dta = data.frame(trait_value = vals, status = factor(status, levels = c("Pruned", "Penalized", "Safe", "Non_Penalized", "Unknown")))
      
      ggplot(dta, aes(x = trait_value, fill = status)) +
        geom_histogram(bins = 40, alpha = 0.8, color = "white") +
        scale_fill_manual(values = c(
          Pruned = "red", Penalized = "orange", Non_Penalized = "grey70", Safe = "forestgreen", Unknown = "lightgrey"
        )) +
        labs(title = paste("Secondary Trait:", toupper(trait)), x = toupper(trait), y = "Count") +
        theme_minimal() + theme(legend.position = "top")
    })
    
    # Add diversity histogram for original population
    original_div_df = original_pop_individual_diversity_scores()
    if (!is.null(original_div_df) && nrow(original_div_df) > 0) {
      mean_orig_div = mean(original_div_df$Diversity, na.rm = TRUE)
      div_plot = ggplot(original_div_df, aes(x = Diversity)) +
        geom_histogram(bins = 40, fill = "lightblue", color = "black", alpha = 0.7) +
        geom_vline(xintercept = mean_orig_div, color = "blue", linetype = "solid", size = 1) +
        geom_vline(xintercept = input$div_target, color = "forestgreen", linetype = "dashed", size = 1) +
        labs(
          title = "Distribution of Original Population Diversity",
          subtitle = paste0("Mean: ", round(mean_orig_div, 2), ", Target: ", input$div_target),
          x = "Diversity Score (1 - COA)",
          y = "Count"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
      
      plots = c(list(div_plot), plots)
    } else {
      warning("Cannot plot original population diversity histogram: No valid individual diversity scores found.")
    }
    
    patchwork::wrap_plots(plots, ncol = 2)
  })
  
  
  
  

  
  
  # ------------------------------------------------------------------
  # GA Logic 
  # ------------------------------------------------------------------
  
  
  reactive_trait_params = reactive({
    current_trait_bounds = list()
    current_penalty_scale = list()
    
    for (trait in all_traits) {
      current_trait_bounds[[trait]] = list(
        threshold = c(
          lower = input[[paste0(trait, "_lower_thresh")]] %||% -Inf,
          upper = input[[paste0(trait, "_upper_thresh")]] %||% Inf
        ),
        prune = c(
          lower = input[[paste0(trait, "_lower_prune")]] %||% -Inf,
          upper = input[[paste0(trait, "_upper_prune")]] %||% Inf
        )
      )
      
      current_penalty_scale[[trait]] = c(
        lower = input[[paste0(trait, "_lower_thres_penalty")]] %||% 0,
        upper = input[[paste0(trait, "_upper_thres_penalty")]] %||% 0
      )
    }
    list(bounds = current_trait_bounds, penalties = current_penalty_scale)
  })
  
  
  
  # Filtered population based on dynamic prune values
  
  pop_dta_filtered_for_ga = reactive({
    params = reactive_trait_params()
    current_pop_dta = pop_dta_initial_processed 
    
    for (trait in all_traits) {
      bounds = params$bounds[[trait]]
      if (is.finite(bounds$prune["lower"])) { # Check if a numeric value was entered
        current_pop_dta = current_pop_dta %>% filter(.data[[trait]] >= bounds$prune["lower"])
      }
      if (is.finite(bounds$prune["upper"])) { # Check if a numeric value was entered
        current_pop_dta = current_pop_dta %>% filter(.data[[trait]] <= bounds$prune["upper"])
      }
    }
    
    current_pop_dta = na.omit(current_pop_dta)
    current_pop_dta
  })
  
  # Penalty logic function
  penalty_logic_reactive = reactive({
    params = reactive_trait_params()
    function(values, bounds, trait_name){ 
      penalty_vec = rep(0, length(values))
      
      lower_thres = bounds$threshold["lower"] %||% -Inf
      upper_thres = bounds$threshold["upper"] %||% Inf
      lower_prune = bounds$prune["lower"] %||% -Inf 
      
      scale_lower = params$penalties[[trait_name]]["lower"] %||% 0
      scale_upper = params$penalties[[trait_name]]["upper"] %||% 0
      
      # Borderline violators get scaled penalty
      okayish_low = which(values >= lower_prune & values < lower_thres)
      okayish_high = which(values > upper_thres & values <= bounds$prune["upper"])
      
      penalty_vec[okayish_low]  = (lower_thres - values[okayish_low]) * scale_lower
      penalty_vec[okayish_high] = (values[okayish_high] - upper_thres) * scale_upper
      
      return(penalty_vec)
    }
  })
  
  ## Fitness Function 
  fitness_function_reactive = reactive({
    params = reactive_trait_params()
    current_pop_dta = pop_dta_filtered_for_ga() 
    current_penalty_logic = penalty_logic_reactive() 
    
    function(x){
      selected = which(x == 1)
      
      # Immediate penalty if no individuals selected or too many
      if (length(selected) == 0) return(-Inf) # No selection is not a valid solution
      
      
      # Check if selected indices are valid for the current_pop_dta
      if (max(selected) > nrow(current_pop_dta) || min(selected) < 1) {
        warning("Invalid indices generated. This should not happen with proper population initialization.")
        return(-Inf)
      }
      
      selected_dta = current_pop_dta[selected, ]
      
      # Penalties for secondary traits
      y2n_penalties = purrr::map_dbl(secondary_traits, function(trait){
        trait_values = selected_dta[[trait]]
        if(length(trait_values) == 0 || all(is.na(trait_values))) return(NA_real_) 
        
        bounds = params$bounds[[trait]] 
        
        penalties = current_penalty_logic(trait_values, bounds, trait)
        mean(penalties, na.rm=TRUE) 
      })
      total_y2n_penalties = sum(y2n_penalties, na.rm = TRUE)
      
      # Primary trait is yield
      norm_y1_colname = paste0('Norm_', input$primary_trait)
      norm_y1 = selected_dta[[norm_y1_colname]]
      if(length(norm_y1) == 0 || all(is.na(norm_y1))) return(-Inf)
      mean_y1 = mean(norm_y1, na.rm=TRUE)
      
      # Diversity Calculation
      all_parents = unique(c(selected_dta$parent1, selected_dta$parent2))
      all_parents = na.omit(all_parents)
      
      if (length(all_parents) == 0 || any(!(all_parents %in% rownames(coa)))) {
        return(-Inf) 
      }
      
      p = setNames(rep(0, length(all_parents)), all_parents)
      for(i in 1:nrow(selected_dta)){
        p[selected_dta$parent1[i]] <- p[selected_dta$parent1[i]] + 0.5
        p[selected_dta$parent2[i]] <- p[selected_dta$parent2[i]] + 0.5
      }
      
      p_vec = unlist(p)
      if (max(p_vec, na.rm = TRUE) == 0) { 
        p_vec_scaled = p_vec
      } else {
        p_vec_scaled = p_vec / max(p_vec, na.rm=TRUE)
      }
      
      
      tmp_pen = t(as.matrix(p_vec_scaled)) %*% coa[all_parents, all_parents]
      if (max(tmp_pen, na.rm = TRUE) == 0) { 
        tmp_pen_scaled = tmp_pen
      } else {
        tmp_pen_scaled = tmp_pen / max(tmp_pen, na.rm=TRUE)
      }
      
      div_penalty = as.numeric((tmp_pen_scaled %*% as.matrix(p_vec_scaled))/100) 
      
      div_score = 0.8 - div_penalty # should be 1 but div is high in this example
      div_bonus = ifelse(div_score >= input$div_target, 5, 0)
      
      # Parent overuse penalty
      parent_count = table(c(selected_dta$parent1, selected_dta$parent2))
      parent_count = parent_count[parent_count > 0]
      if (length(parent_count) > 0) {
        usage_penalty = sd(as.numeric(parent_count)) * overuse_weight 
        if (is.na(usage_penalty)) usage_penalty = 0 # Handle SD of 1 item
      } else {
        usage_penalty = 0
      }
      
      # Final Score calculation
      benefit = mean_y1 * yield_weight
      cost = total_y2n_penalties + div_penalty + usage_penalty
      bonus = div_bonus
      final_score = benefit - cost + bonus
      
      if (is.na(final_score) || !is.finite(final_score)) {
        return(-Inf) # Catch any non-finite 
      }
      
      return(final_score)
    }
  })
  
  # Custom Population function for GA to ensure n_sel is met
  custom_pop_reactive = reactive({
    current_n_sel = input$n_sel
    current_n_bits = nrow(pop_dta_filtered_for_ga())
    
    if (current_n_bits == 0) {
      stop("No individuals available after pruning to create initial population.")
    }
    if (current_n_sel > current_n_bits) {
      stop(paste0("Number to Select (", current_n_sel, ") cannot be greater than the number of available individuals (", current_n_bits, ") after pruning."))
    }
    if (current_n_sel <= 0) {
      stop("Number to Select (n_sel) must be positive.")
    }
    
    function(object){
      pop_mat = replicate(object@popSize, {
        x = rep(0, object@nBits)
        selected_idx = sample(1:object@nBits, current_n_sel, replace = FALSE)
        x[selected_idx] = 1
        x
      })
      t(pop_mat)
    }
  })
  
  
  
  
  # GA Model
  ga_results = eventReactive(input$submit, {

    # validate input
    req(input$n_sel, input$max_iter, input$div_target)
    
    current_fitness_function = fitness_function_reactive()
    current_custom_pop = custom_pop_reactive() 
    current_pop_dta_for_nBits = pop_dta_filtered_for_ga() 
    
    # Re-validate crucial parameters just before GA run
    if (nrow(current_pop_dta_for_nBits) == 0) {
      showNotification("No individuals remaining after pruning. Adjust prune bounds or data.", type = "error")
      return(NULL)
    }
    if (input$n_sel <= 0 || input$n_sel > nrow(current_pop_dta_for_nBits)) {
      showNotification(paste0("Number to Select (", input$n_sel, ") must be between 1 and the number of available individuals (", nrow(current_pop_dta_for_nBits), ") after pruning."), type = "error")
      return(NULL)
    }
    
   
  
    withProgress(message = 'Running GA Optimization', value = 0, {
      incProgress(0.1, detail = "Initializing GA...")
      ga_model = GA::ga(
        type='binary',
        fitness=current_fitness_function,
        nBits=nrow(current_pop_dta_for_nBits),
        popSize=pop_size,
        population=current_custom_pop,
        maxiter=input$max_iter,
        run=input$max_iter,
        keepBest=TRUE,
        monitor=F 
      )
      incProgress(1, detail = "GA Complete!") 
    })
    return(ga_model)
  })
  
  
  
  #  dataframe for plots - fitness, yield, diversity for solutions
  plot_data_df = reactive({
    req(ga_results())
    
    ga_mod = ga_results()
    current_pop_dta_for_metrics = pop_dta_filtered_for_ga()
    current_fitness_function = fitness_function_reactive()
    
    num_pop_individuals = nrow(ga_mod@population)
    if (num_pop_individuals == 0) return(NULL) # No population, no plot data
    
    sample_size = min(500, num_pop_individuals) # Sample up to 500 individuals for plots in case its a lot
    sampled_indices = sample(1:num_pop_individuals, sample_size, replace = FALSE)
    
    fitness_values = numeric(sample_size)
    yield_values = numeric(sample_size)
    diversity_values = numeric(sample_size)
    
    withProgress(message = 'Preparing plot data...', value = 0, {
      for (k in seq_along(sampled_indices)) {
        i = sampled_indices[k]
        incProgress(1/sample_size)
        
        x_chromo = ga_mod@population[i, ]
        selected_inds = which(x_chromo == 1)
        
        # Basic validation of selected_inds for safety
        if (length(selected_inds) == 0 || max(selected_inds) > nrow(current_pop_dta_for_metrics) || min(selected_inds) < 1) {
          fitness_values[k] = NA
          yield_values[k] = NA
          diversity_values[k] = NA
          next
        }
        
        current_selected_dta = current_pop_dta_for_metrics[selected_inds, ]
        
        current_yield_val = mean(current_selected_dta[[input$primary_trait]], na.rm = TRUE)
        
        current_diversity_val = {
          all_parents_curr = unique(c(current_selected_dta$parent1, current_selected_dta$parent2))
          all_parents_curr = na.omit(all_parents_curr)
          if (length(all_parents_curr) == 0 || any(!(all_parents_curr %in% rownames(coa)))) {
            NA_real_
          } else {
            p_curr = setNames(rep(0, length(all_parents_curr)), all_parents_curr)
            for(idx in 1:nrow(current_selected_dta)){
              p_curr[current_selected_dta$parent1[idx]] = p_curr[current_selected_dta$parent1[idx]] + 0.5
              p_curr[current_selected_dta$parent2[idx]] = p_curr[current_selected_dta$parent2[idx]] + 0.5
            }
            p_vec_curr = unlist(p_curr)
            p_vec_curr_scaled = p_vec_curr / max(p_vec_curr, na.rm=TRUE) 
            
            tmp_pen_curr = t(as.matrix(p_vec_curr_scaled)) %*% coa[all_parents_curr, all_parents_curr]
            tmp_pen_curr_scaled = tmp_pen_curr / max(tmp_pen_curr, na.rm=TRUE)
            
            # Handle case where tmp_pen_curr is all zeros -if coa is all zeros for selected parents
            if (any(is.nan(tmp_pen_curr_scaled))) tmp_pen_curr_scaled = 0
            
            div_penalty_curr = as.numeric((tmp_pen_curr_scaled %*% as.matrix(p_vec_curr_scaled))/100)
            round((1 - div_penalty_curr), 2)
          }
        }
        
        current_fitness_score = tryCatch(
          current_fitness_function(x_chromo),
          error = function(e) NA_real_ # Catch errors and set fitness to NA
        )
        
        fitness_values[k] = current_fitness_score
        yield_values[k] = current_yield_val
        diversity_values[k] = current_diversity_val
      }
    })
    
    df_plot = data.frame(
      fitness = fitness_values,
      yield = yield_values,
      diversity = diversity_values
    ) %>% na.omit() 
    
    return(df_plot)
  })
  
  
  # Calculate primary trait mean for each best solution in model
  all_best_solution_primary_trait_means = reactive({
    req(ga_results()) # Ensure the GA has run
    ga_mod = ga_results()
    current_pop_dta_for_metrics = pop_dta_filtered_for_ga() 
    
    # output type caused some issues. I have handled them here.
    tmp_list_of_chromosomes = list()
    if (is.list(ga_mod@bestSol) && length(ga_mod@bestSol) > 0) {
      # It's already a list, typically of 1-row matrices
      
      tmp_list_of_chromosomes = ga_mod@bestSol
    } else if (is.matrix(ga_mod@bestSol) && nrow(ga_mod@bestSol) > 0) {
      
      # If it's a matrix 
      for (r in 1:nrow(ga_mod@bestSol)) {
        tmp_list_of_chromosomes[[r]] <- matrix(ga_mod@bestSol[r, ], nrow = 1)
      }
    } else if (is.vector(ga_mod@bestSol) && length(ga_mod@bestSol) > 0) {
      # If it's a single vector
      tmp_list_of_chromosomes[[1]] <- matrix(ga_mod@bestSol, nrow = 1)
    } else {
      return("Unexpected structure for model or no solutions.")
    }
    
    if (length(tmp_list_of_chromosomes) == 0) {
      return("No best solutions found in ga_model@bestSol.")
    }
    
    mean_primary_traits_for_all_best_solutions = numeric(length(tmp_list_of_chromosomes))
    
    for (i in 1:length(tmp_list_of_chromosomes)) {
      if (!is.null(tmp_list_of_chromosomes[[i]]) && is.matrix(tmp_list_of_chromosomes[[i]]) && nrow(tmp_list_of_chromosomes[[i]]) >= 1) {
        best_solution_vector = tmp_list_of_chromosomes[[i]][1, ]
        
        selected_indices = which(best_solution_vector == 1)
        
        if (length(selected_indices) > 0) {
          # Validate indices against current_pop_dta_for_metrics
          if (max(selected_indices) > nrow(current_pop_dta_for_metrics) || min(selected_indices) < 1) {
            mean_primary_traits_for_all_best_solutions[i] <- NA_real_
            warning("Indices out of bounds for current_pop_dta_for_metrics in best solutions calculation.")
            next
          }
          selected_data <- current_pop_dta_for_metrics[selected_indices, ]
          mean_primary_traits_for_all_best_solutions[i] = mean(selected_data[[input$primary_trait]], na.rm = TRUE)
        } else {
          mean_primary_traits_for_all_best_solutions[i] = NA_real_ # No individuals selected
        }
      } else {
        mean_primary_traits_for_all_best_solutions[i] = NA_real_ # Invalid element in tmp_list
      }
    }
    return(mean_primary_traits_for_all_best_solutions)
  })
  
  # Best selected data
  best_selected_data = reactive({
    req(ga_results())
    ga_mod = ga_results()
    current_pop_dta_for_metrics = pop_dta_filtered_for_ga()
    
    best_chromosome = NULL
    if (is.list(ga_mod@bestSol) && length(ga_mod@bestSol) > 0) {
      best_chromosome <- ga_mod@bestSol[[1]][1, ]
    } else if (is.matrix(ga_mod@bestSol) && nrow(ga_mod@bestSol) > 0) {
      best_chromosome <- ga_mod@bestSol[1, ]
    } else if (is.vector(ga_mod@bestSol) && length(ga_mod@bestSol) > 0) {
      best_chromosome <- ga_mod@bestSol
    }
    
    if (is.null(best_chromosome) || all(is.na(best_chromosome))) {
      return(NULL)
    }
    
    selected_indices_best = which(best_chromosome == 1)
    if (length(selected_indices_best) == 0 ||
        max(selected_indices_best) > nrow(current_pop_dta_for_metrics) ||
        min(selected_indices_best) < 1) {
      return(NULL)
    }
    
    current_pop_dta_for_metrics[selected_indices_best, ]
  })
  
  
  
  
  best_solution_parent_counts = reactive({
    req(ga_results())
    ga_mod = ga_results()
    current_pop_dta_for_metrics = pop_dta_filtered_for_ga()
    
    best_chromosome = NULL
    if (is.list(ga_mod@bestSol) && length(ga_mod@bestSol) > 0) {
      best_chromosome <- ga_mod@bestSol[[1]][1, ]
    } else if (is.matrix(ga_mod@bestSol) && nrow(ga_mod@bestSol) > 0) {
      best_chromosome <- ga_mod@bestSol[1, ]
    } else if (is.vector(ga_mod@bestSol) && length(ga_mod@bestSol) > 0) {
      best_chromosome <- ga_mod@bestSol
    }
    
    if (is.null(best_chromosome) || all(is.na(best_chromosome))) {
      return(NULL)
    }
    
    selected_indices_best = which(best_chromosome == 1)
    if (length(selected_indices_best) == 0 ||
        max(selected_indices_best) > nrow(current_pop_dta_for_metrics) ||
        min(selected_indices_best) < 1) {
      return(NULL)
    }
    
    selected_dta_best = current_pop_dta_for_metrics[selected_indices_best, ]
    parent_count_best = table(c(selected_dta_best$parent1, selected_dta_best$parent2))
    parent_count_best = parent_count_best[parent_count_best > 0] # Remove counts of 0 if any
    return(as.data.frame(parent_count_best, stringsAsFactors = FALSE) %>%
             rename(Parent = Var1, Count = Freq))
  })
  
  
  # Diversity for selected pop
  selected_pop_individual_diversity_scores <- reactive({
    selected_dta <- best_selected_data()
    if (is.null(selected_dta) || nrow(selected_dta) == 0) return(NULL)
    
    scores = calculate_individual_cross_diversity(selected_dta, coa)
    data.frame(Diversity = scores) %>% filter(!is.na(Diversity))
  })
  
  
  
  
  # ------------------------------------------------------------------
  # Outputs
  # ------------------------------------------------------------------
  
  # Tab 2
  
  # Output for GA Diagnostics Log 
  output$log = renderPrint({
    req(ga_results())
    summary(ga_results())
  })
  
  # Output for GA Fitness Progression Plot
  output$fitnessPlot = renderPlot({
    req(ga_results())
    plot(ga_results(), main = "Fitness Progression", legend = FALSE)
    legend("bottomright", legend = c("Best", "Mean", "Median"),
           col = c("green", "blue", "lightgreen"),
           lty = c(1, 2, 1),
           pch = c(16, 1, NA),
           cex = 0.9,
           y.intersp = 0.7
    )
  })
  
  # Output for Scatter Plot (Fitness vs. Primary Trait)
  output$scatterPlot = renderPlot({
    df_plot = plot_data_df()
    req(df_plot)
    
    # Find the row with the maximum fitness
    max_fitness_row = df_plot[which.max(df_plot$fitness), ]
    
    ggplot(df_plot, aes(x = .data[["yield"]], y = fitness)) +
      geom_point(alpha = 0.6) + 
      geom_point(data = max_fitness_row, aes(x = .data[["yield"]], y = fitness, color = "Highest Fitness"),
                 size = 5, shape = 16) + 
      scale_color_manual(name = "Highlighted Points", values = c("Highest Fitness" = "red")) + 
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      labs(title = paste("Fitness vs. Mean", input$primary_trait),
           x = paste("Mean", input$primary_trait, "Value"),
           y = "Fitness Score") +
      theme_minimal() +
      theme(legend.position = "bottom") 
  })
  
  
  # Output for Diversity vs. Primary Trait (colored by Fitness)
  output$divYieldPlot = renderPlot({
    dta_for_plot = plot_data_df()
    req(dta_for_plot)
    ggplot(dta_for_plot, aes(x = .data[["yield"]], y = diversity, color = fitness)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_color_gradient(low = "yellow", high = "darkgreen") +
      labs(title = paste("Mean", input$primary_trait, "vs. Diversity (Colored by Fitness)"),
           x = paste("Mean", input$primary_trait, "Value"),
           y = "Diversity Score",
           color = "Fitness") +
      theme_minimal()
  })
  


  
  
  # Output for Selection Differential Viz
  output$selection_differential_plot = renderPlot({
    selected_data = best_selected_data()
    original_data_for_traits = pop_dta_raw
    req(selected_data, original_data_for_traits)
    
    traits_to_plot = c(input$primary_trait, secondary_traits)
    
    
    
    # --- Plotting function for individual trait density plots ---
    create_trait_density_plot = function(trait_name, original_df, selected_df, common_legend_name) {
      df_plot = bind_rows(
        original_df %>% dplyr::select(Value = all_of(trait_name)) %>% mutate(Group = "Original Population"),
        selected_df %>% dplyr::select(Value = all_of(trait_name)) %>% mutate(Group = "Selected Crosses")
      ) %>%
        filter(!is.na(Value))
      
      mean_original = mean(original_df[[trait_name]], na.rm = TRUE)
      mean_selected = mean(selected_df[[trait_name]], na.rm = TRUE)
      
      ggplot(df_plot, aes(x = Value, fill = Group)) +
        geom_density(alpha = 0.5, size = 1.2, aes(color = Group)) + 
        geom_vline(aes(xintercept = mean_original, color = "Original Mean"),
                   linetype = "dashed", size = 1) +
        geom_vline(aes(xintercept = mean_selected, color = "Selected Mean"),
                   linetype = "solid", size = 1) +
        scale_fill_manual(name = common_legend_name, 
                          values = c("Original Population" = "lightgrey",
                                     "Selected Crosses" = "skyblue")) +
        scale_color_manual(name = common_legend_name, 
                           values = c("Original Population" = "grey",
                                      "Selected Crosses" = "steelblue",
                                      "Original Mean" = "grey",
                                      "Selected Mean" = "steelblue")) +
        labs(
          title = paste("Distribution of", trait_name),
          x = trait_name,
          y = "Density"
        ) +
        theme_minimal() +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) 
    }
    
    plot_list = lapply(traits_to_plot, function(trait) {
      create_trait_density_plot(trait, original_data_for_traits, selected_data, "Population Group/Mean")
    })
    
    # --- Diversity Density Plot ---
    original_div_df_for_plot = original_pop_individual_diversity_scores()
    selected_div_df_for_plot = selected_pop_individual_diversity_scores()
    
    div_density_plot = NULL 
    if (!is.null(original_div_df_for_plot) && !is.null(selected_div_df_for_plot) &&
        nrow(original_div_df_for_plot) > 0 && nrow(selected_div_df_for_plot) > 0) {
      df_div_plot = bind_rows(
        original_div_df_for_plot %>% mutate(Group = "Original Population"),
        selected_div_df_for_plot %>% mutate(Group = "Selected Crosses")
      ) %>%
        filter(!is.na(Diversity))
      
      mean_original_div = mean(original_div_df_for_plot$Diversity, na.rm = TRUE)
      mean_selected_div = mean(selected_div_df_for_plot$Diversity, na.rm = TRUE)
      
      div_density_plot = ggplot(df_div_plot, aes(x = Diversity, fill = Group)) +
        geom_density(alpha = 0.5, size = 1.2, aes(color = Group)) +
        geom_vline(aes(xintercept = mean_original_div, color = "Original Mean"),
                   linetype = "dashed", size = 1) +
        geom_vline(aes(xintercept = mean_selected_div, color = "Selected Mean"),
                   linetype = "solid", size = 1) +
        geom_vline(aes(xintercept = input$div_target, color = "Diversity Target"),
                   linetype = "dotted", size = 1) +
        scale_fill_manual(name = "Population Group/Mean/Target", 
                          values = c("Original Population" = "lightgrey",
                                     "Selected Crosses" = "skyblue")) +
        scale_color_manual(name = "Population Group/Mean/Target", 
                           values = c("Original Population" = "grey", 
                                      "Selected Crosses" = "steelblue",   
                                      "Original Mean" = "grey",
                                      "Selected Mean" = "steelblue",
                                      "Diversity Target" = "forestgreen")) +
        labs(
          title = "Distribution of Individual Cross Diversity Scores",
          x = "Diversity Score (1 - COA(parent1, parent2))",
          y = "Density"
        ) +
        theme_minimal() +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) 
    } else {
      warning("Cannot plot diversity densities in Selection Differential: No valid individual diversity scores found for original or selected populations.")
    }
    
    # --- Parent Usage Histogram ---
    parent_counts_df = best_solution_parent_counts()
    parent_usage_plot = NULL
    if (!is.null(parent_counts_df) && nrow(parent_counts_df) > 0) {
      parent_counts_df = parent_counts_df %>%
        mutate(OverLimit = Count > input$parent_limit)
      
      parent_usage_plot = ggplot(parent_counts_df, aes(x = Count, fill = OverLimit)) +
        geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
        scale_fill_manual(name = "Parent Usage", 
                          values = c("TRUE" = "lightcoral", "FALSE" = "skyblue")) +
        labs(
          title = "Distribution of Parent Usage",
          x = "Times Parent Used",
          y = "Number of Parents",
          fill = paste0("Exceeds Limit (>", input$parent_limit, ")")
        ) +
        theme_minimal() +
        theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
        scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
    } else {
      warning("Cannot plot parent usage histogram: No valid parent count data.")
    }
    
    
    # --- Overused Parents Text Output ---
    
    overused_parents_text <- reactive({
      parent_counts_df = best_solution_parent_counts()
      if (is.null(parent_counts_df) || nrow(parent_counts_df) == 0) return("No parent usage data available.")
      
      overused_parents = parent_counts_df %>%
        filter(Count > input$parent_limit) %>%
        arrange(desc(Count))
      
      if (nrow(overused_parents) > 0) {
        text_content = paste0("Parents exceeding usage limit (>", input$parent_limit, "):\n",
                               paste(capture.output(print(overused_parents)), collapse = "\n"))
      } else {
        text_content = paste0("No parents exceeded the usage limit (max ", input$parent_limit, ").")
      }
      return(text_content)
    })
    
    
    
    overused_parents_grob = grid::grobTree(grid::textGrob(overused_parents_text(),
                                                           x = 0.05, y = 0.95, just = c("left", "top"),
                                                           gp = grid::gpar(fontsize = 10)))
    
    
    # --- Arrange Plots with Patchwork ---
    first_row_plots = wrap_plots(plot_list, ncol = 4, common.legend = TRUE, legend = "top") 
    
    second_row_plots = wrap_plots(
      div_density_plot,
      parent_usage_plot,
      overused_parents_grob,
      ncol = 3, 
      widths = c(1, 1, 0.7)
    )
    
    # Combine both rows
    first_row_plots / second_row_plots 
  })
  
  
  # Table for best solution
  selected_individuals_table_data = reactive({
    selected_data = best_selected_data() %>% dplyr::select(!starts_with('Norm_'))
    names(selected_data) = toupper(names(selected_data))
    if (!is.null(selected_data) && nrow(selected_data) > 0) {
      return(selected_data)
    }
    return(data.frame(Message = "No best solution found or invalid indices."))
  })
  
  
  output$summaryTable = DT::renderDataTable({
    selected_individuals_table_data()
    
  }, options = list( pageLength = 10, scrollX = TRUE),
  rownames = FALSE)
  
  
  output$downloadSummary = downloadHandler(
    filename = function() {
      paste("selected_individuals_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(selected_individuals_table_data(), file, row.names = FALSE)
    }
  )
  
  
  
}


shinyApp(ui = ui, server = server)