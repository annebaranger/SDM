"4th_step_R_coordinate_cleaning_step_by_step"

This script enables the cleaning of the occurrence data obtained before.
It removes records that :
- have non-numeric and non available coordinates
- have lat >90, la <-90, lon > 180 or lon < -180
- show an erroneous conversion from a degree minute format to a decimal degree format
- have either zero longitude or latitude and a radius of 0.1 decimal degree around the point at zero longitude and zero latitude
- are within a 10km radius around country capitals
- are within a 1 kmradius around the geographic centroids of political countries and provinces
- are duplicated
- have equal longitude and latitude coordinates
- are within 0.5 degree radius around the GBIF headquarters in Copenhagen
- are assigned to the location of zoos, botanical gardens, herbaria, universities and museums
- are outside the referenced landmass
- are outliers in geographic space (data entry errors, imprecise geo-references)
Finally the result is a table of correctec occurrences per species.
