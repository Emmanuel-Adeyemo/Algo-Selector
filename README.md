# Algo-Selector
A shiny application that is designed to assist in selecting an optimal subset of crosses from a larger simulated crosses to create new breeding populations based on multiple genetic traits and diversity constraints. It leverages a Genetic Algorithm (GA) to maximize the primary trait while adhering to stipulated thresholds for the secondary traits, achieving the target genetic diversity among the selected group (if realistic/achievable), and limiting parental overuse. Instead of manually checking through thousands of potential crosses and trying to balance multiple traits, this app automates the complex selection process. It searches a vast solution space to find the best set of individuals based on the defined constraints, significantly reducing time, especially when working with a large list of potential crosses.

## Fitness Logic
Data Normalization: Numeric trait columns are automatically normalized for consistent fitness calculation.

Fitness Function: The heart of the GA - this is what the algorithm seeks to maximize. It calculates a final score for each potential solution (a set of selected crosses).

$finalscore = Mean(Primary Trait) - Secondary Trait Penalties - Diversity Penalty - Parent Overuse Penalty + Diversity Bonus$

Secondary Trait Penalties: Based on user definition. Values outside thresholds incur scaled penalties. Pruned values are automatically discarded.

Diversity Penalty: Calculated from the coa matrix for the selected parents. A diversity score lower than the diversity target results in a penalty.

Parent Overuse Penalty: Based on the standard deviation of how many times each parent is used in the selected group, discouraging reliance on a few parents. 

Diversity Bonus: A fixed bonus is awarded if the overall diversity score of the selected population meets or exceeds the diversity target.

## Running the app from your R console
To run the app locally, install the following R packages:

```r
install.packages(shiny)
install.packages(bslib)
install.packages(GA)
install.packages(dplyr)
install.packages(purrr) 
install.packages(janitor)
install.packages(ggplot2)
install.packages(GGally)
install.packages(patchwork)
```

Then run directly in your R console using:
```r
shiny::runGitHub("Algo-Selector", username = "emmanuel-adeyemo", launch.browser = TRUE)
```
