"11th_step_Ensemble_of_small_species"

This script enables the modelling of species with less than  occurrences after the 20% removal for outer validation.
This time we used the procedure proposed by ... et al. 
We compared a global maxent model to an ensemble of small model (build as a weighted average of bivariate Maxent models).
Then, using the best one, you can project species distribution model considering present (prediction) and future (projection) climate.
This script creates a projected species disctribution map (with presence/absence), per RCP and scenarios, for each species. 