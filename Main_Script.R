# Build aster models to test for fitness effects of dispersal distance and maternal performance for the three focal Lasthenia species. 
# Will need to bring in fitness data for both experimental individuals and maternal plants for all three years of the experiment 

# libraries
library(aster)
library(tidyverse)

# Bring in 2018 experimental & maternal data #
# experimental
data_2018 <- read_csv("data/Jepson_2018_Traits.csv")
data_2018$year <- 2018

#maternal
maternal_2018 <- read_csv("data/Jepson_2017_Maternal.csv")
maternal_2018$year <- 2018

###########################################################
# Bring in 2019 experimental data 
#experimental
data_2019 <- read_csv("data/Jepson_2019_Traits.csv")
data_2019$year <- 2019

#maternal
maternal_2019 <- read_csv("data/JP_2018_Maternal.csv")
maternal_2019$year <- 2019

###########################################################
# Bring in 2022 experimental data 
#experimental
source("scripts/Combining_2022SurveySample.R") # creates prelim data_2022 dataframe

#maternal 
maternal_2022 <- read_csv("data/JP_2020_Maternal.csv")
maternal_2022$year <- 2022

###########################################################
# Modifydataframes so they talk nicely to each other
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

data_2022 <- data_2022[ ,c("group","id", "sp", "tran", "dist", "trt","loc","year", "status", "vs")]
data_2022 <-
    data_2022 %>% 
    filter(!is.na(status)) %>% 
    unite(col = "mom", remove = FALSE, sp, group, id) %>% 
    unite(col = "indiv",remove = FALSE, group, sp, id, tran, dist, trt, loc, year)
data_2022

# combine all data
all_data <- rbind(data_2018, data_2019, data_2022)
all_data


#####################################################################
# make vs an integer
all_data$vs <- as.integer(all_data$vs)

# reshape for aster - need it to be a long-form where there is a row for each event - for survival and for the number of seeds 
re_data<-
 all_data %>%
 pivot_longer(11:12, names_to = "vars", values_to = "resp") %>% 
  arrange(vars)

# change vars (status, vs) to a factor
re_data$vars <- as.factor(re_data$vars)
re_data

# add needed columns/vectors
re_data$root <- 1
pred <- c(0, 1)
fam <- c(1, 2); sapply(fam.default(), as.character)[fam]

# make a fitness vector
re_data$fit <- as.numeric(re_data$vars == "vs")

#Make a numeric 'id' row
nindvs_2019 <- 1:nrow(unique(re_all_field_data[re_all_field_data$year == 2019, "indiv"]))
nindvs_2018 <-1:nrow(unique(re_all_field_data[re_all_field_data$year == 2018, "indiv"]))
nindvs_2018 <- max(nindvs_2019) + nindvs_2018
re_all_field_data$id <- c(nindvs_2019 , nindvs_2018 , nindvs_2019, nindvs_2018)
re_all_field_data$id  <- as.integer(re_all_field_data$id)
