# Algo-Selector
A shiny application that is designed to assist in selecting an optimal subset of crosses from a larger population based on multiple genetic traits and diversity constraints. It leverages a Genetic Algorithm (GA) to maximize the primary trait while adhering to stipulated thresholds for the secondary traits, achieving the target genetic diversity among the selected group, and limiting parental overuse.

## Fitness Logic
Data Normalization: Numeric trait columns are automatically normalized for consistent fitness calculation.

Fitness Function: The heart of the GA - this is what the algorithm seeks to maximize. It calculates a final_score for each potential solution (a set of selected crosses).

$final_score = Mean(Primary Trait) - Secondary Trait Penalties - Diversity Penalty - Parent Overuse Penalty + Diversity Bonus$

Secondary Trait Penalties: Based on user definition. Values outside thresholds incur scaled penalties. Pruned values are automatically discarded.

Diversity Penalty: Calculated from the coa matrix for the selected parents. A diversity score lower than the diversity target results in a penalty.

Parent Overuse Penalty: Based on the standard deviation of how many times each parent is used in the selected group, discouraging reliance on a few parents. 

Diversity Bonus: A fixed bonus is awarded if the overall diversity score of the selected population meets or exceeds the diversity target.
