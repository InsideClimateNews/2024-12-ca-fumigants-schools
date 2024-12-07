# set working directory to the folder containing this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load required packages
library(tidyverse)
library(sf)
library(tigris)
library(janitor)
library(tidycensus)
library(units)
library(corrplot)
library(leaflet)
library(htmlwidgets)

# don't use scientific notation for numbers
options(scipen = 999)

# cache TIGER shapefiles
options(tigris_use_cache = TRUE)

############################
# geodata

# load PLSS section-level map
# from California Department of Pesticide Regulation https://calpip.cdpr.ca.gov/plssFiles.cfm, cleaned up in QGIS
# this large file excluded from online version of Github repo 
ca_plss <- st_read("processed_data/plss/ca_plss.geojson") %>%
  clean_names() 

############################
# load and process Pesticide Use Reporting data 
# from California Department of Pesticide Regulation https://files.cdpr.ca.gov/pub/outgoing/pur_archives/
# extensive data files excluded from online version of Github repo 

# check chemical codes/names
chemical_files <- list.files(path = "data/pur", recursive = TRUE, full.names = TRUE, pattern = "chemical")
chemicals <- map_dfr(chemical_files, read_csv) %>%
  unique()

write_csv(chemicals, "processed_data/pur/chemical.csv", na = "")

# load PUR data from 2017-2022, filtering for 1,3-D and chloropicrin using their chem_code
pur_files <- list.files(path = "data/pur", recursive = TRUE, full.names = TRUE, pattern = "udc")
pur <- map_dfr(pur_files, ~ read_csv(.x, col_types = cols(.default = col_character())) %>%
                 filter(chem_code == "573" | chem_code == "136")) %>%
  unique() 

# some initial processing for correct data types, and to define calendar and California water years
pur <- pur %>%
  mutate(applic_dt = mdy(applic_dt),
         lbs_chm_used = as.double(lbs_chm_used),
         acre_treated = as.double(acre_treated),
         year = year(applic_dt),
         water_year = if_else(month(applic_dt) >= 10, year(applic_dt) + 1, year(applic_dt)))

# look at error flags
error_files <- list.files(path = "data/pur", recursive = TRUE, full.names = TRUE, pattern = "errors2")
errors <- map_dfr(error_files, ~ {
  year <- str_extract(.x, "\\d{4}") 
  read_csv(.x, col_types = cols(.default = col_character())) %>% 
    mutate(year = as.double(year))
})

pur_errors <- inner_join(pur, errors, by = c("use_no","year"))
write_csv(pur_errors, "processed_data/pur/pur_errors.csv", na = "") # exported for scrutiny

# remove any rows with error flags
pur <- pur %>%
  filter(is.na(error_flag))

# separate data frames for 1,3-D and chlorpicrin
pur_13d <- pur %>%
  filter(chem_code == "573")

pur_chloropicrin <- pur %>%
  filter(chem_code == "136")

############################
# 5-year fumigant use at PLSS section level for Mapbox interactive map

######
# 1,3-D

# filter for 2018-2022 California water years
pur_2018_2022_13d <- pur_13d %>%
  filter(water_year < 2023 & water_year > 2017) %>%
  mutate(comtrs = paste0(county_cd,base_ln_mer,township,tship_dir,range,range_dir,section)) 

# calculate total applications in pounds
pur_2018_2022_13d_total <- pur_2018_2022_13d %>%
  group_by(comtrs)%>%
  summarize(pounds_13d = sum(lbs_chm_used, na.rm = TRUE)) 

# join to sections map, calculate area of sections in acres and applications normalized by area
pur_2018_2022_13d_total_map <- inner_join(ca_plss, pur_2018_2022_13d_total, by = c("co_mtrs" = "comtrs")) %>%
   select(co_mtrs,pounds_13d) %>%
   mutate(section_total_acres = as.double(set_units(st_area(.), "acres")),
          pounds_13d_per_100_acres_total = pounds_13d/section_total_acres*100) %>%
   st_transform("EPSG:4326")

# quantiles, for breaks used on Mapbox map
quantile(pur_2018_2022_13d_total_map$pounds_13d_per_100_acres_total, probs = seq(0,1,0.2), na.rm = TRUE)
# 0%           20%           40%           60%           80%          100% 
#   0.1864698   690.9809856  1375.9936832  2391.0829394  4305.3990637 60805.6869800 

# export for Mapbox map
st_write(pur_2018_2022_13d_total_map, "processed_data/pesticide_applications/13d/water_yrs_2018_2022_13d.geojson")

######
# chloropicrin

# filter for 2017-2022 California water years
pur_2018_2022_chloropicrin <- pur_chloropicrin %>%
  filter(water_year < 2023 & water_year > 2017) %>%
  mutate(comtrs = paste0(county_cd,base_ln_mer,township,tship_dir,range,range_dir,section))

# calculate total applications in pounds
pur_2018_2022_chloropicrin_total <- pur_2018_2022_chloropicrin %>%
  group_by(comtrs)%>%
  summarize(pounds_chloropicrin = sum(lbs_chm_used, na.rm = TRUE)) 

# join to sections map, calculate area of sections in acres and applications normalized by area
pur_2018_2022_chloropicrin_total_map <- inner_join(ca_plss, pur_2018_2022_chloropicrin_total, by = c("co_mtrs" = "comtrs")) %>%
  select(co_mtrs,pounds_chloropicrin) %>%
  mutate(section_total_acres = as.double(set_units(st_area(.), "acres")),
         pounds_chloropicrin_per_100_acres_total = pounds_chloropicrin/section_total_acres*100) %>%
  st_transform("EPSG:4326")

# quantiles, for breaks used on Mapbox map
quantile(pur_2018_2022_chloropicrin_total_map$pounds_chloropicrin_per_100_acres_total, probs = seq(0,1,0.2), na.rm = TRUE)
# 0%            20%            40%            60%            80%           100% 
#   0.1129576    137.5960693    496.3009123   1408.2859818   6020.3266300 179219.5585856 

# export for Mapbox map
st_write(pur_2018_2022_chloropicrin_total_map, "processed_data/pesticide_applications/chloropicrin/water_yrs_2018_2022_chloropicrin.geojson")


############################
# process data for intersections in QGIS to allow estimates of applications near schools

######
# 1,3-D applications

by_water_yr_2018_2022_13d <- pur_2018_2022_13d %>%
  group_by(comtrs,water_year) %>%
  summarize(pounds_13d = sum(lbs_chm_used, na.rm = TRUE)) %>%
  arrange(water_year,comtrs,-pounds_13d)

by_water_yr_2018_2022_13d_map <- inner_join(ca_plss, by_water_yr_2018_2022_13d, by = c("co_mtrs" = "comtrs"))

for (y in unique(by_water_yr_2018_2022_13d$water_year)) {
  m <- by_water_yr_2018_2022_13d_map %>%
    filter(water_year == y)
  st_write(m, paste0("processed_data/pesticide_applications/13d/water_yr_",y,"_13d.geojson"))
}

######
# chloropicrin applications

by_water_yr_2018_2022_chloropicrin <- pur_2018_2022_chloropicrin %>%
  group_by(comtrs,water_year) %>%
  summarize(pounds_chloropicrin = sum(lbs_chm_used, na.rm = TRUE)) %>%
  arrange(-pounds_chloropicrin)

by_water_yr_2018_2022_chloropicrin_map <- inner_join(ca_plss, by_water_yr_2018_2022_chloropicrin, by = c("co_mtrs" = "comtrs"))

for (y in unique(by_water_yr_2018_2022_chloropicrin$water_year)) {
  m <- by_water_yr_2018_2022_chloropicrin_map %>%
    filter(water_year == y)
  st_write(m, paste0("processed_data/pesticide_applications/chloropicrin/water_yr_",y,"_chloropicrin.geojson"))
}

######
# crop survey data
# from California Department of Water Resources annual Statewide Crop Mapping surveys https://data.cnra.ca.gov/dataset/statewide-crop-mapping
# this survey uses water years, hence the choice of water years above
# these large files excluded from online version of Github repo 

crop_files <- list.files(path = "data/ca_crops", recursive = TRUE, full.names = TRUE, pattern = "\\.shp$")

for (f in crop_files) {
    c <- st_read(f) %>%
      clean_names() 
    c <- c %>%
      # filter for agricultural land on which the pesticides might be applied
      filter((class1 == "G" | class1 == "R" | class1 == "F" | class1 == "P" | class1 == "T" | class1 == "D" | class1 == "C" | class1 == "V" | class1 == "YP" | (class1 == "I" & subclass1 == "2"))
             | (class2 == "G" | class2 == "R" | class2 == "F" | class2 == "P" | class2 == "T" | class2 == "D" | class2 == "C" | class2 == "V" | class2 == "YP" | (class2 == "I" & subclass2 == "2"))
             | (class3 == "G" | class3 == "R" | class3 == "F" | class3 == "P" | class3 == "T" | class3 == "D" | class3 == "C" | class3 == "V" | class3 == "YP" | (class3 == "I" & subclass3 == "2"))
             | (class4 == "G" | class4 == "R" | class4 == "F" | class4 == "P" | class4 == "T" | class4 == "D" | class4 == "C" | class4 == "V" | class4 == "YP" | (class4 == "I" & subclass4 == "2")))
  y <- regmatches(f, gregexpr("\\d{4}", f))
  y <- sapply(y, function(x) if(length(x) > 0) x[1] else NA)
  st_write(c, paste0("processed_data/crops/crops",y,".gpkg"))
}
# these large files saved to the processed data are excluded from online version of Github repo

###################################################
# intersections performed in QGIS
###################################################

###################################################
# estimation of fumigant applications in buffer zones of 0.25 and 1 miles from school boundaries 

######
# 1,3-D
crops_13d_files <- list.files("processed_data/intersections/13d/crops_13d/", full.names = TRUE) 
crops_13d <- crops_13d_files %>%
  lapply(st_read) %>%
  bind_rows() %>%
  mutate(crops_acres = as.double(set_units(st_area(.), "acres"))) %>%
  select(co_mtrs,water_year,crops_acres, pounds_13d) %>%
  st_drop_geometry()

crops_13d %>% get_dupes() # no dupes
 
# 0.25 mile buffer
schools_0.25_13d_files <- list.files("processed_data/intersections/13d/buffer_0.25", full.names = TRUE)
schools_0.25_13d <- schools_0.25_13d_files %>%
  lapply(st_read) %>%
  bind_rows() %>%
  clean_names() %>%
  mutate(intersection_acres = as.double(set_units(st_area(.), "acres")))  %>%
  select(2:13,stacked,stack_cnt,grades_offered,grades_served,co_mtrs,water_year,pounds_13d,intersection_acres) %>%
  st_drop_geometry()

schools_0.25_13d_join <- inner_join(schools_0.25_13d,crops_13d, by = c("water_year","co_mtrs")) %>%
  mutate(intersection_fraction = intersection_acres/crops_acres,
         pounds_13d_intersection = intersection_fraction*pounds_13d.x)
# all intersection fractions are between 0 and 1

schools_0.25_13d_2018_2022 <- schools_0.25_13d_join %>%
  group_by(across(c(1:16))) %>%
  summarize(pounds_13d_intersection = sum(pounds_13d_intersection)) %>%
  arrange(-pounds_13d_intersection)

write_csv(schools_0.25_13d_2018_2022,"processed_data/pesticie_applications_schools/schools_0.25_13d_2018_2022.csv",na = "")

# 1 mile buffer
schools_1_13d_files <- list.files("processed_data/intersections/13d/buffer_1", full.names = TRUE)
schools_1_13d <- schools_1_13d_files %>%
  lapply(st_read) %>%
  bind_rows() %>%
  clean_names() %>%
  mutate(intersection_acres = as.double(set_units(st_area(.), "acres")))  %>%
  select(2:13,stacked,stack_cnt,grades_offered,grades_served,co_mtrs,water_year,pounds_13d,intersection_acres) %>%
  st_drop_geometry() %>%
  unique() # this handles dupes in schools data for Summit Charter Academy 54718370109009

schools_1_13d_join <- inner_join(schools_1_13d,crops_13d, by = c("water_year","co_mtrs")) %>%
  mutate(intersection_fraction = intersection_acres/crops_acres,
         pounds_13d_intersection = intersection_fraction*pounds_13d.x)
# all intersection fractions are between 0 and 1

schools_1_13d_2018_2022 <- schools_1_13d_join %>%
  group_by(across(c(1:16))) %>%
  summarize(pounds_13d_intersection = sum(pounds_13d_intersection)) %>%
  arrange(-pounds_13d_intersection)

write_csv(schools_1_13d_2018_2022,"processed_data/pesticide_applications_schools/schools_1_13d_2018_2022.csv",na = "")

######
# chloropicrin
crops_chloropicrin_files <- list.files("processed_data/intersections/chloropicrin/crops_chloropicrin/", full.names = TRUE) 
crops_chloropicrin <- crops_chloropicrin_files %>%
  lapply(st_read) %>%
  bind_rows() %>%
  mutate(crops_acres = as.double(set_units(st_area(.), "acres"))) %>%
  select(co_mtrs,water_year,crops_acres, pounds_chloropicrin) %>%
  st_drop_geometry()

# 0.25 mile buffer
schools_0.25_chloropicrin_files <- list.files("processed_data/intersections/chloropicrin/buffer_0.25", full.names = TRUE)
schools_0.25_chloropicrin <- schools_0.25_chloropicrin_files %>%
  lapply(st_read) %>%
  bind_rows() %>%
  clean_names() %>%
  mutate(intersection_acres = as.double(set_units(st_area(.), "acres")))  %>%
  select(2:13,stacked,stack_cnt,grades_offered,grades_served,co_mtrs,water_year,pounds_chloropicrin,intersection_acres) %>%
  st_drop_geometry()

schools_0.25_chloropicrin_join <- inner_join(schools_0.25_chloropicrin,crops_chloropicrin, by = c("water_year","co_mtrs")) %>%
  mutate(intersection_fraction = intersection_acres/crops_acres,
         pounds_chloropicrin_intersection = intersection_fraction*pounds_chloropicrin.x)
# all intersection fractions are between 0 and 1

schools_0.25_chloropicrin_2018_2022 <- schools_0.25_chloropicrin_join %>%
  group_by(across(c(1:16))) %>%
  summarize(pounds_chloropicrin_intersection = sum(pounds_chloropicrin_intersection)) %>%
  arrange(-pounds_chloropicrin_intersection)

write_csv(schools_0.25_chloropicrin_2018_2022,"processed_data/pesticide_applications_schools/schools_0.25_chloropicrin_2018_2022.csv",na = "")

# 1 mile buffer
schools_1_chloropicrin_files <- list.files("processed_data/intersections/chloropicrin/buffer_1", full.names = TRUE)
schools_1_chloropicrin <- schools_1_chloropicrin_files %>%
  lapply(st_read) %>%
  bind_rows() %>%
  clean_names() %>%
  mutate(intersection_acres = as.double(set_units(st_area(.), "acres")))  %>%
  select(2:13,stacked,stack_cnt,grades_offered,grades_served,co_mtrs,water_year,pounds_chloropicrin,intersection_acres) %>%
  st_drop_geometry()

schools_1_chloropicrin %>%
  get_dupes() # good, no dupes

schools_1_chloropicrin_join <- inner_join(schools_1_chloropicrin,crops_chloropicrin, by = c("water_year","co_mtrs")) %>%
  mutate(intersection_fraction = intersection_acres/crops_acres,
         pounds_chloropicrin_intersection = intersection_fraction*pounds_chloropicrin.x)
# all intersection fractions are between 0 and 1

schools_1_chloropicrin_join %>%
  filter(pounds_chloropicrin.x != pounds_chloropicrin.y) # gives zero rows, as it should

write_csv(schools_1_chloropicrin_join,"processed_data/pesticide_applications_schools/schools_1_chloropicrin_year.csv", na = "")

schools_1_chloropicrin_2018_2022 <- schools_1_chloropicrin_join %>%
  group_by(across(c(1:16))) %>%
  summarize(pounds_chloropicrin_intersection = sum(pounds_chloropicrin_intersection)) %>%
  arrange(-pounds_chloropicrin_intersection)

write_csv(schools_1_chloropicrin_2018_2022,"processed_data/pesticide_applications_schools/schools_1_chloropicrin_2018_2022.csv",na = "")

#################
# processing the above for Datawrapper tables on fumigant applications near schools

# school enrollment data for 2021-2022, from California Department of Education https://data-cdegis.opendata.arcgis.com/datasets/712403d542894040a3ec01281cc2ebaf_0/explore
enrollment <- read_csv("data/schools/SchoolSites2122.csv") %>%
  clean_names() %>%
  select(cds_code,enroll_total)

schools_13d_2018_2022 <- full_join(
  schools_0.25_13d_2018_2022 %>%
    rename(pounds_0.25 = pounds_13d_intersection),
  schools_1_13d_2018_2022 %>%
    rename(pounds_1 = pounds_13d_intersection)
) %>%
  left_join(enrollment, by = "cds_code") %>%
  ungroup() %>%
  select(cds_code,school,district,enroll_total,pounds_0.25,pounds_1) %>%
  arrange(-pounds_1)

write_csv(schools_13d_2018_2022, "processed_data/pesticide_applications_schools/schools_13d_2018_2022.csv", na = "")

######
# chloropicrin

schools_chloropicrin_2018_2022 <- full_join(
  schools_0.25_chloropicrin_2018_2022 %>%
    rename(pounds_0.25 = pounds_chloropicrin_intersection),
  schools_1_chloropicrin_2018_2022 %>%
    rename(pounds_1 = pounds_chloropicrin_intersection)
) %>%
  left_join(enrollment, by = "cds_code") %>%
  ungroup() %>%
  select(cds_code,school,district,enroll_total,pounds_0.25,pounds_1) %>%
  arrange(-pounds_1)

write_csv(schools_chloropicrin_2018_2022, "processed_data/pesticide_applications_schools/schools_chloropicrin_2018_2022.csv", na = "")

############################
# locations of schools in the above data, for Mapbox interactive map
# from California School Campus Database https://www.mapcollaborator.org/mapcollab_cscd/
# this data excluded from online version of Github repo

st_layers("data/schools/CSCD_2021.gdb")

schools_centroids <- st_read("data/schools/CSCD_2021.gdb", layer = "School_Centroids") %>%
  clean_names() %>%
  select(cds_code) %>%
  st_transform("EPSG:4326") %>%
  group_by(cds_code) %>%
  slice_head(n = 1) # handles duplicate points for Summit Charter Academy 54718370109009

schools_13d_2018_2022_sf <- inner_join(schools_centroids,schools_13d_2018_2022, by = "cds_code") 
st_write(schools_13d_2018_2022_sf,"processed_data/pesticide_applications_schools/schools_13d_2018_2022.geojson")

schools_chloropicrin_2018_2022_sf <- inner_join(schools_centroids,schools_chloropicrin_2018_2022, by = "cds_code")
st_write(schools_chloropicrin_2018_2022_sf,"processed_data/pesticide_applications_schools/schools_chloropicrin_2018_2022.geojson")

############################
# analysis at Census tract level, allowing pesticide exposures to be considered in socioeconomic context

# load CDC Social Vulnerability Index data (based on 2018-2022 American Community Survey)
# data from https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html
st_layers("data/svi/SVI2022_CALIFORNIA_tract.gdb")

svi <- st_read("data/svi/SVI2022_CALIFORNIA_tract.gdb", layer = "SVI2022_CALIFORNIA_tract") %>%
  clean_names() %>%
  select(county,fips,area_sqmi,
         housing_units = e_hu,
         households = e_hh,
         svi = rpl_themes,
         socioeconomic_svi = rpl_theme1,
         household_svi = rpl_theme2,
         racial_svi = rpl_theme3,
         housing_transport_svi = rpl_theme4) %>%
  mutate(across(where(is.numeric) & !Shape, ~ replace(., . == -999, NA)))

# load 2028-2022 American Community Survey data using tidycensus

v22 <- load_variables(2022, "acs5", cache = TRUE)

variables <- c(
  population = "B01001_001", # denominator for age/sex and for born in Central America
  poverty = "B17001_002",
  population_assessed_poverty = "B17001_001", # denominator for poverty
  born_central_america = "B05006_154",
  population_5_and_older = "B16005_001", # denominator for limited English
  median_household_income = "B19013_001",
  limited_english_1 = "B16005_007",
  limited_english_2 = "B16005_008",
  limited_english_3 = "B16005_012",
  limited_english_4 = "B16005_013",
  limited_english_5 = "B16005_017",
  limited_english_6 = "B16005_018",
  limited_english_7 = "B16005_022",
  limited_english_8 = "B16005_023",
  limited_english_9 = "B16005_029",
  limited_english_10 = "B16005_030",
  limited_english_11 = "B16005_034",
  limited_english_12 = "B16005_035",
  limited_english_13 = "B16005_039",
  limited_english_14 = "B16005_040",
  limited_english_15 = "B16005_044",
  limited_english_16 = "B16005_045",
  male_under_5_yrs = "B01001_003",
  male_5_9_yrs = "B01001_004",
  male_10_14_yrs = "B01001_005",
  male_15_17_yrs = "B01001_006",
  female_under_5_yrs = "B01001_027",
  female_5_9_yrs = "B01001_028",
  female_10_14_yrs = "B01001_029",
  female_15_17_yrs = "B01001_030",
  latino = "B03003_003"
)

acs_tracts <- get_acs(
  geography = "tract",
  state = "CA",
  variables = variables,
  year = 2022,
  survey = "acs5"
) %>%
clean_names()

acs_state <- get_acs(
  geography = "state",
  state = "CA",
  variables = variables,
  year = 2022,
  survey = "acs5"
) %>%
  clean_names()

acs_tracts <- acs_tracts %>%
  select(-moe) %>%
  group_by(geoid) %>%
  pivot_wider(names_from = variable, values_from = estimate)

acs_tracts <- acs_tracts %>%
  ungroup() %>%
  mutate(age_17_and_below = rowSums(select(., contains("yrs")), na.rm = TRUE),
         limited_english = rowSums(select(., contains("limited")), na.rm = TRUE),
         pc_age_17_and_below = round(age_17_and_below/population*100,1),
         pc_poverty = round(poverty/population_assessed_poverty*100,1),
         pc_latino = round(latino/population*100,1),
         pc_limited_english = round(limited_english/population_5_and_older*100,1),
         pc_poverty = round(poverty/population_assessed_poverty*100,1),
         pc_born_central_america = round(born_central_america/population*100,1))

acs_tracts <- acs_tracts %>%
  select(geoid,name,population, median_household_income, contains("pc_"))

acs_state <- acs_state %>%
  select(-moe) %>%
  group_by(geoid) %>%
  pivot_wider(names_from = variable, values_from = estimate)

acs_state <- acs_state %>%
  ungroup() %>%
  mutate(age_17_and_below = rowSums(select(., contains("yrs")), na.rm = TRUE),
         pc_poverty = round(poverty/population_assessed_poverty*100,1),
         limited_english = rowSums(select(., contains("limited")), na.rm = TRUE),
         pc_age_17_and_below = round(age_17_and_below/population*100,1),
         pc_poverty = round(poverty/population_assessed_poverty*100,1),
         pc_latino = round(latino/population*100,1),
         pc_limited_english = round(limited_english/population_5_and_older*100,1),
         pc_born_central_america = round(born_central_america/population*100,1))

acs_state <- acs_state %>%
  select(geoid,name,population, median_household_income, contains("pc_"))

write_csv(acs_state, "processed_data/acs_state.csv", na = "")

# comparison of tract level to state level data
acs_tracts_relative <- sweep(acs_tracts[4:9], 2, as.numeric(acs_state[4:9]), "/") # simple ratio

acs_tracts_relative_long <- acs_tracts_relative %>%
  pivot_longer(cols = 1:6, names_to = "variable", values_to = "value")

colnames(acs_tracts_relative) <- paste0(colnames(acs_tracts_relative), "_relative")

acs_tracts <- bind_cols(acs_tracts,acs_tracts_relative)

# combine SVI and ACS data
svi_acs <- inner_join(svi,acs_tracts, by = c("fips" = "geoid"))

# filter pesticide applicaton data for 2018-2022 calendar years
cal_yrs_2018_2022_13d <- pur_13d %>%
  filter(year < 2023 & year > 2017) %>%
  mutate(comtrs = paste0(county_cd,base_ln_mer,township,tship_dir,range,range_dir,section)) %>%
  group_by(comtrs) %>%
  summarize(acres_treated = sum(acre_treated, na.rm = TRUE),
            pounds_13d = sum(lbs_chm_used, na.rm = TRUE))

cal_yrs_2018_2022_chloropicrin <- pur_chloropicrin %>%
  filter(year < 2023 & year > 2017) %>%
  mutate(comtrs = paste0(county_cd,base_ln_mer,township,tship_dir,range,range_dir,section)) %>%
  group_by(comtrs) %>%
  summarize(acres_treated = sum(acre_treated, na.rm = TRUE),
            pounds_chloropicrin = sum(lbs_chm_used, na.rm = TRUE))

# aggregate pesticide data to census tracts

sf_use_s2(FALSE)

# load data for plss/tract intersection, processed in QGIS
ca_plss_tracts <- st_read("processed_data/plss/ca_plss_tracts.geojson") %>%
  clean_names()

# calculate areas
ca_plss_tracts <- ca_plss_tracts %>%
  mutate(area_intersect = as.double(set_units(st_area(.), "acres")))

# calculate areas and drop geometry
ca_plss_areas <- ca_plss %>%
  mutate(area = as.double(set_units(st_area(.), "acres"))) %>%
  st_drop_geometry()

# calculate proportion of each intersected area in each section
ca_plss_tracts <- ca_plss_tracts %>%
  inner_join(ca_plss_areas, by = "co_mtrs") %>%
  mutate(area_proportion = area_intersect / area)

# cdpr county codes
cdpr_county_codes <- read_csv("data/county_codes/cdpr_county_codes.csv") %>%
  mutate(county_code = sprintf("%02d", county_code))

# join to pesticide application data and aggregate by tract
ca_tracts_2018_2022_13d <- ca_plss_tracts %>%
  inner_join(cal_yrs_2018_2022_13d, by = c("co_mtrs" = "comtrs")) %>%
  mutate(weighted_pounds_13d = pounds_13d * area_proportion) %>%
  st_drop_geometry() %>%
  group_by(geoid) %>%
  summarize(pounds_13d_2018_2022 = sum(weighted_pounds_13d,na.rm = TRUE))

# join to chloropicrin 2018-2022 and aggregate by tract
ca_tracts_2018_2022_chloropicrin <- ca_plss_tracts %>%
  inner_join(cal_yrs_2018_2022_chloropicrin, by = c("co_mtrs" = "comtrs")) %>%
  mutate(weighted_pounds_chloropicrin = pounds_chloropicrin * area_proportion) %>%
  st_drop_geometry() %>%
  group_by(geoid) %>%
  summarize(pounds_chloropicrin_2018_2022 = sum(weighted_pounds_chloropicrin,na.rm = TRUE))

# join pesticide to svi/acs data

ca_tracts <- tracts(state = "CA", cb = TRUE) %>%
  st_transform("EPSG:3310")

ca_tracts_land_area_acres <- ca_tracts %>%
  clean_names() %>%
  select(geoid, aland) %>%
  mutate(land_area_acres = aland * 0.000247105) %>%
  st_drop_geometry() %>%
  select(-aland)

ca_tracts_svi_acs_pesticides <- svi_acs %>%
  left_join(ca_tracts_land_area_acres, by = c("fips" = "geoid")) %>%
  left_join(ca_tracts_2018_2022_13d, by = c("fips" = "geoid")) %>%
  left_join(ca_tracts_2018_2022_chloropicrin, by = c("fips" = "geoid")) %>%
  mutate(pounds_13d_per_100_acres_2018_2022 = pounds_13d_2018_2022/land_area_acres*100,
         pounds_chloropicrin_per_100_acres_2018_2022 = pounds_chloropicrin_2018_2022/land_area_acres*100,
         pounds_13d_per_1000_people_2018_2022 = pounds_13d_2018_2022/population*1000,
         pounds_chloropicrin_per_1000_people_2018_2022 = pounds_chloropicrin_2018_2022/population*1000) %>%
  rename(geometry = Shape) %>%
  select(county,
         fips,
         name,
         area_sqmi,
         land_area_acres,
         population,
         housing_units,
         households,
         svi,
         socioeconomic_svi,
         household_svi,
         racial_svi,
         housing_transport_svi,
         median_household_income,
         pc_poverty,
         pc_age_17_and_below,
         pc_latino,
         pc_limited_english,
         pc_born_central_america,
         median_household_income_relative,
         pc_age_17_and_below_relative,
         pc_poverty_relative,
         pc_latino_relative,
         pc_limited_english_relative,
         pc_born_central_america_relative,
         pounds_13d_2018_2022,
         pounds_chloropicrin_2018_2022,
         pounds_13d_per_100_acres_2018_2022,
         pounds_chloropicrin_per_100_acres_2018_2022,
         pounds_13d_per_1000_people_2018_2022,
         pounds_chloropicrin_per_1000_people_2018_2022,
         geometry)

write_csv(ca_tracts_svi_acs_pesticides %>% st_drop_geometry(), "processed_data/socioeconomic/ca_tracts_svi_acs_pesticides.csv", na = "")
write_csv(ca_tracts_svi_acs_pesticides %>% st_drop_geometry(), "processed_data/socioeconomic/ca_tracts_svi_acs_pesticides.csv", na = "")

st_write(ca_tracts_svi_acs_pesticides, "processed_data/socioeconomic/ca_tracts_svi_acs_pesticides.geojson")

# correlation matrix and plot
cor_matrix <- ca_tracts_svi_acs_pesticides %>%
  st_drop_geometry() %>%
  select(9:19,contains("per")) %>%
  cor(use = "complete.obs") # use "complete.obs" to handle NA values

corrplot(cor_matrix,
         method = "color",
         type = "lower",
         addCoef.col = 'gray',
         tl.col = "black",
         tl.cex = 0.7,
         number.cex = 0.4)





