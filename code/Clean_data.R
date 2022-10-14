#--------------------------------------------------------------------------------------------
# Intrahabitat complexity drives the distribution of fish trait groups on coral reefs
# By Noelle Helder
# Updated 4/5/2021
# 
# This script outputs cleaned, analysis-ready data frames (.csv) from raw fish survey data and organizes an intra-reef 
# complexity metrics csv. 

#--------------------------------------------------------------------------------------------
# 1. Clean, organize, and merge survey data frames (fish, cluster, habitat); create survey data frame that matches site-grid-survey # with
#     observer data, restoration status (y/n); update language to match manuscript (grid = plot; site = reefscape; CC = CR2, etc.)

# 2. Calculate fish community summaries per plot: abundance (inds/25 m2), biomass (g/m2), 
#     diversity (Shannon), richness (# of species)

# 3. Calculate functional group level summaries per plot: abundance (inds/25 m2)

# 4. Organize data frame of final plot-level reef complexity predictors: VRM_1cm, VRM_4cm (VRM_4cm-VRM_1cm), 
#    absolute profile curvature, height_diff
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# Summary data set descriptions
#
# fish -  fish observations per plot from FL Keys 2019
#         @survey a unique ID for each site, visit, and transect 
#         @grid 25m^2 unit of observation per reefscape. Value 1-96 for NDR1, NDR2, GR1, GR2 and 1-100 for CR1, CR2, PR1, PR2
#         @common_name common name of fish species
#         @life_phase observed life phase (A=adult, J=juvenile, I=intermediate)
#         @total_length total length in cm
#         @transit transit or NA if individual was moving through the plot rather than remaining within

# species - fish species data derived from Fish Base to match with common_name from fish data 

# survey - metadata for each unique survey 

# coords - relative x,y coordinates for plots 1-100 within each reefscape. 0,0 originates at plot 1

# complex - reef complexity data per plot measured from DEMs at 1cm resolution using script from Fukunaga et al. 2021

# cluster.load - results from reef fish trait-based clustering results

# transect - metadata to align transect info with plot and reefscape
#--------------------------------------------------------------------------------------------

# Load Packages: 
library(devtools)
library(tidyverse)
library(vegan)
library(cowplot)
library(viridis)
library(plyr)
library(lubridate)
library(corrplot)

# Load data 
fish.load <- read_csv(here::here('./data/raw/fish_surveys.csv'))
cluster.load <- read_csv(here::here('./data/raw/functional_groups.csv'))
complex.load <- read_csv(here::here('./data/raw/intrareef_complexity_metrics.csv'))

# survey metadata
species <- read_csv(here::here('./data/raw/species_list.csv'))
survey <- read_csv(here::here('./data/raw/survey_data.csv'))
coords <- read_csv(here::here('./data/raw/plot_coordinates.csv'))
transect <- read_csv(here::here('./data/raw/transect_id.csv'))


#####################################################################################################
# 1. DATA CLEANING + PREP
# Combine info from multiple .csv files to creat 1 survey  metadata file associated with each plot
#####################################################################################################
# Drop index column
species <- species %>% dplyr::select(-'...1')

# rename var to match with survey df
transect$trans <- transect$transect
transect <- transect %>% dplyr::select(-'transect')

# Calculate survey time
survey$time.st.pos <- as.POSIXct(survey$time_start) # convert start time
survey$time.end.pos <- as.POSIXct(survey$time_end) # convert end time
survey$time <- (survey$time.end.pos-survey$time.st.pos) # diff


# add variable for restoration status at the site level (Restored vs. Control) and Protection status 
# (Protected, Not Protected)
site <- as.factor(c("CR", "GR", "NDR", "PR", "CC", "GC", "NDRC", "PC"))
protection <- c("Y", "Y", "N", "N", "Y", "Y", "N", "N")
restored <- c("Y", "Y", "Y", "Y", "N", "N", "N", "N")
site.df <- data.frame(site, protection, restored)
coords <- coords %>% left_join(site.df)
coords$grid <- as.integer(coords$grid) 

# Create DF includes a row for every grid that was surveyed during the field season. 
# Final community aggregated df should match w/ 4308 rows (with 0s) + 13 vars. 
metadata <- (survey) %>% dplyr::select(survey, site, trans, site_survey, surveyor, current, time) %>% 
  separate(site_survey, c("abbrev", "visit")) %>% 
  left_join(transect) %>% 
  dplyr::select(-'site') %>% 
  dplyr::rename(site=abbrev) %>% 
  dplyr::left_join(coords) 

# confirm total # observations = 4308


m <- metadata
m$grid[m$grid %in% "1"] <- "01"
m$grid[m$grid %in% "2"] <- "02"
m$grid[m$grid %in% "3"] <- "03"
m$grid[m$grid %in% "4"] <- "04"
m$grid[m$grid %in% "5"] <- "05"
m$grid[m$grid %in% "6"] <- "06"
m$grid[m$grid %in% "7"] <- "07"
m$grid[m$grid %in% "8"] <- "08"
m$grid[m$grid %in% "9"] <- "09"
metadata <- m

# export metadata file
#write.csv(metadata, "metadata.csv")  


# CLEAN + PREP FISH DATA: 
# change fish data plot numbers (this helps later on). 
fish <- fish.load
m <- fish
m$grid[m$grid %in% "1"] <- "01"
m$grid[m$grid %in% "2"] <- "02"
m$grid[m$grid %in% "3"] <- "03"
m$grid[m$grid %in% "4"] <- "04"
m$grid[m$grid %in% "5"] <- "05"
m$grid[m$grid %in% "6"] <- "06"
m$grid[m$grid %in% "7"] <- "07"
m$grid[m$grid %in% "8"] <- "08"
m$grid[m$grid %in% "9"] <- "09"
fish <- m

# check
# table(fish$grid)

# change coords dataset grid numbers. 
m <- coords
m$grid[m$grid %in% "1"] <- "01"
m$grid[m$grid %in% "2"] <- "02"
m$grid[m$grid %in% "3"] <- "03"
m$grid[m$grid %in% "4"] <- "04"
m$grid[m$grid %in% "5"] <- "05"
m$grid[m$grid %in% "6"] <- "06"
m$grid[m$grid %in% "7"] <- "07"
m$grid[m$grid %in% "8"] <- "08"
m$grid[m$grid %in% "9"] <- "09"
coords <- m


# Fix data entry issues to match with species list

# NAs to correct for: grouper sp., pompano sp., Porgy sp. 
fish[fish == "Grouper sp."] <- "Black Grouper" #unknown grouper? use data for black grouper
fish[fish == "Pompano sp."] <- "Florida Pompano" #unknown species coded as Pompano. 
fish[fish == "Porgy sp."] <- "Saucereye Porgy" #unknown porgy? use data for Saucereye porgy
fish[fish == "Chub sp."] <- "Brassy chub" #unknown chub? use data for brassy chub. 
fish[fish == "Slippery dick"] <- "Slippery Dick" #rename as needed to align with species list.


# Exclude animals that are not part of analysis - only interested in resident fish species >15cm
filter_fish <- fish %>% 
  dplyr::filter(common_name != "Reef shark") %>% 
  dplyr::filter(common_name != "Nurse Shark") %>% 
  dplyr::filter(common_name != "Tarpon") %>% 
  dplyr::filter(common_name != "Southern Stingray") %>% 
  dplyr::filter(common_name != "Yellow Stingray") %>% 
  dplyr::filter(common_name != "Permit") %>% 
  dplyr::filter(common_name != "Florida Pompano") %>% 
  dplyr::filter(common_name != "Brown Chromis") %>% 
  dplyr::filter(common_name != "Blue Chromis") %>% 
  dplyr::filter(common_name != "Yellowtail Damselfish") %>% 
  dplyr::filter(common_name != "Yellowtail Hamlet") %>% 
  dplyr::filter(common_name != "Blue Hamlet") %>% 
  dplyr:: filter(common_name != "Black Hamlet") %>% 
  dplyr::filter(common_name != "Barred Hamlet") %>% 
  dplyr::filter(common_name != "Butter Hamlet") %>% 
  dplyr::filter(common_name != "Mackerel Scad") %>% # only seen once
  dplyr::filter(common_name != "Blue spotted cornetfish") %>%  # only seen once + not included in sp. matrix
  dplyr::filter(common_name != "Yellow Jack")

# Check common names for spelling, weird occurrences, etc. 
unique(filter_fish$common_name) 
length(unique(filter_fish$common_name)) #82 unique species observed that meet criteria. 

# Link species name with species info from FishBase spreadsheet 
dat <- merge(filter_fish, species, by.x = "common_name", 
             by.y = "common_name", all.x = TRUE, all.y = FALSE)
length(unique(dat$family))

# Fix incorrect a and b values for Queen Parrotfish (as found during herbivory analysis)
dat$a <- ifelse(dat$common_name %in% 'Queen Parrotfish', 0.0144,dat$a)
dat$b <- ifelse(dat$common_name %in% 'Queen Parrotfish', 3.05, dat$b)

# Join all dataframes
dat <- dat %>% 
  left_join(survey) %>% 
  left_join(coords) %>% 
  separate(site_survey, c("abbrev", "visit")) %>% 
  
  
  # Filter small fish (<15cm): leaving in transit for now to compare
  dplyr::filter(total_length >= 15) %>% 
  # test <- dat %>% dplyr::select(common_name, trophic_group) %>% unique()
  dplyr::filter(is.na(transit)) #remove any species in transit. 
                                # This also includes large resting
                                #school of schoolmasters at PC.

# Select variables to use in data vis/analysis
names(dat)
dat <- dat %>% dplyr::select(survey, abbrev, grid, visit, common_name, total_length, life_phase,
                             a, b) # 29,680 observations


# Clean and prepare habitat metrics csv 

# Recode grids to match other dataframes
m <- complex
m$grid[m$grid %in% "1"] <- "01"
m$grid[m$grid %in% "2"] <- "02"
m$grid[m$grid %in% "3"] <- "03"
m$grid[m$grid %in% "4"] <- "04"
m$grid[m$grid %in% "5"] <- "05"
m$grid[m$grid %in% "6"] <- "06"
m$grid[m$grid %in% "7"] <- "07"
m$grid[m$grid %in% "8"] <- "08"
m$grid[m$grid %in% "9"] <- "09"
complex <- m

complex$grid <- as.factor(complex$grid)
levels(complex$grid)

# Separate by resolution (not needed, but will keep in for now)
res1 <- complex %>% filter(fac == 1)
res2 <- complex %>% filter(fac == 2)
res4 <- complex %>% filter(fac == 4)
# res8 <- complex %>% filter(fac == 8)
# res16 <- complex %>% filter(fac == 16)
# res32 <- complex %>% filter(fac == 32)

# mean TL = 21.06 cm 
summary(dat$total_length) # 21.06

###################################################################################################
# INTRA-REEF COMPLEXITY CSV (habitat predictors)
###################################################################################################
# pull out just 1 cm resolution data
predictors <- res1 %>% dplyr::select(site, grid, mean_vrm, height_diff, mean_curv, 
                                        mean_pro_curv, mean_plan_curv)

# Rename mean_VRM to vrm_1cm to contrast with vrm_4cm (from res4 df)
colnames(predictors)[3] <- "vrm_1cm"

# Add VRM 4cm to analysis dataframe
temp <- res4 %>% dplyr::select(site, grid, mean_vrm) # subset res4 df
colnames(temp)[3] <- "vrm_4cm" # rename
predictors <- predictors %>% left_join(temp) # join

# Calculate VRM deviation (4cm) from Fukunaga et al. 2021
# 4cm VRM and 1cm VRM are highly corrlated: p=0.93
corrplot(cor(predictors[3:8], use = "complete.obs"), method="number", type = "upper",
         diag=FALSE,number.cex= 7/ncol(predictors), title = "Final Metrics", mar=c(0,0,2,0))

# Calculate difference VRM deviation following Fukunaga et al. 2020: now it is uncorrelated. 
# This has been shown to capture complexity from mounding and boulder type corals, etc: 
predictors$vrm_dev <- (predictors$vrm_4cm-predictors$vrm_1cm)

# Check that vrm_dev and vrm_1cm are no longer correlated
corrplot(cor(predictors[3:9], use = "complete.obs"), method="number", type = "upper",
         diag=FALSE,number.cex= 7/ncol(predictors), title = "Final Metrics", mar=c(0,0,2,0))

# drop extra variable 
predictors <- predictors %>% dplyr::select(-'vrm_4cm')

# Export final complexity data
# write.csv(predictors, "intrareef_complexity_clean.csv")


#####################################################################################################
# FISH COMMUNITY METRICS ()
# 1. Abundance 2. Density 3. Diversity/Richness  4. Biomass 
#####################################################################################################
dat$site <- dat$abbrev

#1) Total Abundance --------------------------------------
abund <- dat %>% 
  dplyr::group_by(site, visit, grid) %>%
  dplyr::summarise(abundance=n()) %>% 
  dplyr::ungroup()

#Create a unique identifying column 
#abund_mean$uniqueid <- paste(abund_mean$site, abund_mean$grid, sep="_")
abund$uniqueid <- paste(abund$site, abund$grid, abund$visit, sep="_")

#2) Community Density ------------------------------------------
area <- res1 %>% dplyr::select(site, grid, area) 

abund.area <- abund %>% left_join(area)

# Missing area values?. Swap out NAs for 25m2
abund.area$area <- ifelse(abund.area$area %in% 
                            NA, 25,abund.area$area)

# Calculate density
density <- abund.area %>% 
  dplyr::mutate(density=(abundance/area))


#3) Richness and Diversity  --------------------------------------- 
# Calculate abundance for each species to create species x abundance matrix
community <- dat %>% 
  dplyr::group_by(site, grid, common_name, visit) %>% 
  dplyr::summarize(abundance=n()) %>% 
  ungroup() #%>% 
#dplyr::group_by(site, grid, common_name)# %>% 
# dplyr::summarize(abundance.mean=mean(abundance, na.rm=TRUE))

# create unique identifier col.
community$uniqueid <- paste(community$site,
                            community$grid, community$visit, sep="_")
community.wide <- community %>% 
  spread(key=common_name, value=abundance) #change data from LONG to WIDE

community.wide[is.na(community.wide)] <- 0

# Richness: # of species
SpeciesRichness <- ddply(community.wide, ~uniqueid, function(x) {
  data.frame(Richness=sum(x[-3]>0))
})
SpeciesRichness

# Diversity 
#Calculate Shannon-Wiener Index (H') 
# It is an information index and is the most commonly used diversity index in ecology.
# Technically, the Shannon-Wiener Index (when applied to ecology) 
# quantifies the uncertainty associated with predicting the identity of a new taxa
# given number of taxa and evenness in abundances of individuals within each taxa 
# Values of H' can range from 0 to 5, although they typically range from 1.5 to 3.5
# The Shannon-Wiener Index assumes that the sample for site was collected randomly.

df <- community.wide[,4:86]

ShannonW <- ddply(df, ~uniqueid, function(x) {
  data.frame(shannon=diversity(x[-1], index="shannon"))
})
ShannonW

# Join with richness
df <- SpeciesRichness %>% left_join(ShannonW)

# Join with other metrics dataframe 
fish.data <- df %>% full_join(density)


#4) Community Biomass (g/m2) --------------------------------- 
bio <-  dat %>% 
  dplyr::mutate(biomass=(a*(total_length^b))) %>% 
  dplyr::filter(!is.na(biomass)) %>% 
  dplyr::group_by(site, visit,grid) %>% 
  dplyr::summarize(biomass=sum(biomass)) %>% 
  dplyr::left_join(area)

# Missing area values for NDRC, CC? Swap out NAs for 25m2 (assumed grid size)
bio$area <- ifelse(bio$area %in% NA, 25,bio$area)

biomass <- bio%>% 
  dplyr::mutate(biomass=biomass/area) %>%
  dplyr::ungroup()
#dplyr::group_by(site, grid) #%>% 
#dplyr::summarize(biodensity.mean=mean(biodensity, na.rm=TRUE), biodensity.sd=sd(biodensity))

fish.data <- fish.data %>% left_join(biomass) %>% 
  dplyr::select(c(uniqueid, site, grid, visit, abundance, density, Richness, 
                  shannon, biomass)) #organize; 4308 rows

# Add 0s: 
total.df <- metadata %>% left_join(fish.data) %>%
  replace_na(list(abundance = 0, biomass = 0,
                  density = 0, Richness = 0, shannon = 0)) %>% 
  dplyr::select(c(site, grid, visit, abundance, density, biomass, Richness, shannon))


# Export analysis-ready community level data
# write.csv(total.df, "fish_community_summaries_clean.csv", row.names=FALSE)

###################################################################################################
# FUNCTIONAL GROUP RESPONSE METRIC
# 1. Abundance 2. Density (inds/m2) 3. Biomass (g/m2) 
###################################################################################################
# Calculate biomass and abundance by functional group (k=9: 
cluster <- cluster.load %>% dplyr::select(common_name, fish.clust.num9)
cluster$comb.clust <- cluster$fish.clust.num9

df <- dat %>% full_join(cluster) %>% filter(common_name != 'Glassy Sweeper') %>% # not included in traits
  filter(common_name != 'Yellow Jack') # not included in analysis

# ABUNDANCE + DENSITY 
abund.cl <- df %>% 
  dplyr::group_by(site, grid, visit, fish.clust.num9) %>% 
  dplyr::summarize(abundance=n()) %>% 
  dplyr::left_join(area)
abund.cl$area <- ifelse(abund.cl$area %in% NA, 25,abund.cl$area)
abund.cl <- abund.cl %>% 
  dplyr::mutate(density=abundance/area)



# BIOMASS: 
bio.cl <-  df %>% 
  dplyr::mutate(biomass=(a*(total_length^b))) %>% 
  dplyr::filter(!is.na(biomass)) %>% 
  dplyr::group_by(site, grid, visit, fish.clust.num9) %>% 
  dplyr::summarize(biomass=sum(biomass)) %>% 
  dplyr::left_join(area)
bio.cl$area <- ifelse(bio.cl$area %in% NA, 25,bio.cl$area)
bio.cl <- bio.cl %>% 
  dplyr::mutate(biomass=biomass/area)

cluster.df <- abund.cl %>% left_join(bio.cl)

# 
func.group.df <- cluster.df %>% 
  dplyr::select(-'area') %>% 
  ungroup() %>% # must ungroup for complete to work. 
  tidyr::complete(nesting(site, grid), fish.clust.num9, visit, fill=list(biomass=0, 
                                                                         abundance=0,
                                                                         density=0)) 
func.group.df$uniqueid <- paste(func.group.df$site, func.group.df$visit,  sep="_") # add unique ID col 
# to drop site-visit pairs that
# do not exist. 
func.group.df<- func.group.df %>% # 
  dplyr::filter(uniqueid != "PC_5") %>% # PC only had 4 visits (drop 5 & 6)
  dplyr::filter(uniqueid != "PC_6") %>% 
  dplyr::filter(uniqueid != "PR_6") %>% # PR only had 5 visits
  dplyr::filter(uniqueid != "GC_6") # GC only had 5 visits; total = 43,080 (=4308*9 groups)

# levels(as.factor(func.group.df$uniqueid))


# Export functional group abundance data
# write.csv(func.group.df, "functional_group_summaries_clean.csv", row.names=FALSE)