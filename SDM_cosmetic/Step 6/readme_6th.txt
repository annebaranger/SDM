"6th_step_species_ecoregions_and_nb_cropping"

This scripts enables the determination of all the ecoregions in which a species has occurrence records. All the neighbouring ecoregions are also added.
The ecoregions used are those from ... 2017
To determine the ecoregions in which the species has occurrence records, the script uses ArcGis.
First it selects all the ecoregions that "CONTAINS" occurrences records.
Second it selects all the ecoregions whose "'BOUNDARIES TOUCH" the first set of selected ecoregions.
Then, it creates one shapefile per species with all the selected ecoregions.
