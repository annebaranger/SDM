"7th_step_R_climate_variable_selection_per_species"

This script selects the climate variables that will be used for species modelling and projection.
The variables used are CHELSA_BIOCLIM_1979_2013.
Among them we removed BIO_2 and BIO_3 for their low interest in explaining our occurrences.
Then the script allow the removal of correlated variables that would lead to errors in modelling and prediction.
This removal prevent our model from overfitting with correlated predicted variables.
The removal is done through a stepwise procedure using the variance inflation factor.
The result of this script is a table with a line per species and, in each species line, you have a some columns in which the selected variables are written.
