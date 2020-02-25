"9th_step_occurrences_aggregation"

This script realizes the aggregation of the corrected occurrences obtained earlier.
This way the occurrences data will have the same "resolution" than the climate data (aggregated earlier) that will be used for projection.
In order to have only one occurrences per raster cells the script used a climate raster (that will be used for projection),
and counted the number of occurrences per cell. If there is more than one, it removes the supplementary ones.
The result is a table of aggregated occurrences per species.

Cimate data are mising here cause it was to big for github but you can find it in BD SIG Chelsa Bioclim ...