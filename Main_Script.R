# Build aster models to test for fitness effects of dispersal distance and maternal performance for the three focal Lasthenia species. 
# Will need to bring in fitness data for both experimental individuals and maternal plants for all three years of the experiment 

# libraries
library(aster)
library(tidyverse)

# Bring in 2018 experimental & maternal data #
# experimental
data_2018 <- read.csv("data/Jepson_2018_Traits.csv")
data_2018$year <- 2018

#maternal
maternal_2018 <- read_csv("data/Jepson_2017_Maternal.csv")
maternal_2018$year <- 2018
###########################################################

# Bring in 2019 experimental data 
#experimental
data_2019 <- read.csv("data/Jepson_2019_Traits.csv")
data_2019$year <- 2019

#maternal
maternal_2019 <- read_csv("data/JP_2018_Maternal.csv")
maternal_2019$year <- 2019
###########################################################

# Bring in 2022 experimental data 
#experimental
data_2022 <- source()
data_2022$year <- 2022

#maternal 
maternal_2022 <- read_csv("data/JP_2020_Maternal.csv")
maternal_2022$year <- 2022
###########################################################

# Change dataframes so they talk nicely to each other any are prepared for required aster model structure
data_2018 <- data_2018[ ,c("group","mayID", "sp", "trans", "dist", "trt","loc","year", "status", "vs")]
colnames(data_2018) <-c("group","id", "sp", "tran", "dist", "trt","loc","year", "status", "vs")
data_2018 <-
    data_2018 %>% 
    filter(!is.na(status)) %>% 
    unite(col = "mom", remove = FALSE, sp, group, id) %>% 
    unite(col = "indiv",remove = FALSE, group, sp, id, tran, dist, trt, loc, year)
data_2018

data_2019 <- data_2019[ ,c("group","id", "sp", "tran", "dist", "trt","loc","year", "status", "vs")]
data_2019 <-
    data_2019 %>% 
    filter(!is.na(status)) %>% 
    unite(col = "mom", remove = FALSE, sp, group, id) %>% 
    unite(col = "indiv",remove = FALSE, group, sp, id, tran, dist, trt, loc, year)
data_2019

# trait dataframe for 2022 only includes processed samples not all location with seeds --> need to find dataframe with status column to get an accurate estimate on how many individuals didnt germinate
data_2022 <- data_2022[ ,c("group","id", "sp", "tran", "dist", "trt","loc","year", "status", "vs")]
data_2019 <-
    data_2019 %>% 
    filter(!is.na(status)) %>% 
    unite(col = "mom", remove = FALSE, sp, group, id) %>% 
    unite(col = "indiv",remove = FALSE, group, sp, id, tran, dist, trt, loc, year)
data_2019

