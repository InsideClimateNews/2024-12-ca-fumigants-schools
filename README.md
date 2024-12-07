# Applications of the agricultural fumigants 1,3-dichloropropene and chloropicrin in California

Data, methodology and and [R](https://www.r-project.org/) code for the analysis underlying this [Inside Climate News](t/k) article, summarizing and mapping use of the agricultural fumigants 1,3-dichloropropene (1,3-D) and chloropicrin in California from 2018 to 2022, including estimates of applications within quarter-mile and one-mile buffer zones extended from the boundaries of school grounds.

Code for the analysis is in the script `fumigants.R` and its output is in the file `fumigants.RData`. The analysis also depended on geospatial buffers and intersections made using [QGIS](https://www.qgis.org/).

Data for applications of the pesticides is from the annual [Pesticide Use Reporting](https://files.cdpr.ca.gov/pub/outgoing/pur_archives/) data released by the California Department of Pesticide Regulation, which records applications at the level of sections in the [Public Land Survey System](https://www.usgs.gov/faqs/do-us-topos-and-national-map-have-a-layer-shows-public-land-survey-system-plss). (A typical PLSS section has an area of about one square mile.)

To consider these applications alongside socioeconomic and demographic variables, we aggregated the data to Census tracts, areas that are home to about 4,000 people. We then joined this to tract-level data from the CDC's 2022 [Social Vulnerability Index](https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html), plus additional variables from the 2018-2022 five-year [American Community Survey](https://www.census.gov/programs-surveys/acs), obtained using the [tidycensus](https://walker-data.com/tidycensus/) R package. Correlations between socioeconomic and demographic variables and the intensity of fumigant use, measured in pounds applied from 2018 to 2022 per unit area, were generally weak, see `correlation_matrix.png` in the `plots` folder. However, tracts with the greatest intensity of fumigant use had high proportions of residents born in Mexico or elsewhere in Central America, with limited proficiency in English, and children aged 17 or younger.

These analyses, plus summaries of fumigant use by county and year, considered applications over the calendar years 2018 to 2022.

For data on school boundaries, we used the [California School Campus Database](https://www.mapcollaborator.org/mapcollab_cscd/). For estimates of applications within buffer zones extended from these boundaries, and for the published interactive map showing fumigant applicantions at the level of PLSS sections and the locations of these schools, we considered applications over the California water years 2018 to 2022. (California water years run from October 1 of the previous calendar year to September 30.) This was because our estimates of applications depended on spatial intersections between buffer zones around the schools and active cropland, which is assessed annually by water year in the [Statewide Crop Mapping](https://data.cnra.ca.gov/dataset/statewide-crop-mapping) survey by the California Department of Water Resources.

For each year, PLSS section, and school buffer zone, we ran intersections between school buffer zones and surveyed active cropland. Applications within buffer zones were then estimated using the following formula, where `intersection_fraction` is the area of of intersection between the buffer zone and active cropland divided by the area of active cropland:

`pounds of fumigant applied in buffer zone = intersection fraction * pounds of fumigant applied in the PLSS section`

This makes the assumption that the fumigants would have been applied evenly over the cropland within any PLSS section and year, which will not literally be true, but is a reasonable approach given the data at our disposal. Because the intersection fractions tended to be larger for the one-mile buffer zones than for the quarter-mile buffer zones, see `intersection_fraction_histograms.png` in the `plots` folder, estimates for the larger buffer zones are likely to be more robust.

### Acknowledgements

Thanks to epidemiologist [Paul English](https://www.phi.org/experts/paul-english/), director of Tracking California, an environmental health initiative of the Public Health Institute in Oakland, for advice on the analysis of applications of the two fumigants near schools.

### Questions/Feedback

Email Peter Aldhous at [peter.aldhous\@insideclimatenews.org](mailto:peter.aldhous@insideclimatenews.org).
