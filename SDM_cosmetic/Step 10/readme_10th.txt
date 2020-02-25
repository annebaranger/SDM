"10th_step_workflow_big_species"

This script allows the modelling and projection of all species with more than 50 occurrences (after removal of 20% kept for outer validation of the model)
In this script we model species distribution with GLM, GBM and Maxent. Each model is implememnted 6 times for each species (3 repetitions with two pseudo-absence subsets).*
The standard models obtained are then combined to obtain one ensemble per type of model (1 ensemble for GLM, 1 for GBM and 1 for Maxent). 
These three ensemble models are made with the 6 repetitions for each model.
Then we implememted an ensemble model considering all types of model (GLM, GBM and Maxent). This time only models with TSS>0.7 where taken into account.
Finally we selected the best model between the 4 ensemble (1 per type and 1 with all types).
After the modelling procedure we projected the best model to present (prediction) and future (projection) climate data (using different RCP and scenarios).
This script creates a projected species distribution map (with presence/absence), per RCP and scenarios, for each species. 
