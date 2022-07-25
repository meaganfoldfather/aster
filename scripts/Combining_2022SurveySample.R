#Combine 2022 survey and sample data

# pull in sample data
sample_2022 <- read_csv("data/Jepson_2022_Samples.csv")
sample_2022

# pull in survey data 
survey_2022 <- read_csv("data/Jepson_2022_Survey.csv")
survey_2022 <-
survey_2022 %>% 
  pivot_longer(cols = L1:L4, names_to = "locs")

# mark status as 1 if germinated (B,F,V,S) and NA if missing (M,G)
survey_2022[survey_2022$value == 0, "status"] <- 0
survey_2022[survey_2022$value %in% c("V","B","F","f"), "status"] <- 1
survey_2022[survey_2022$value %in% c("M","G"), "status"] <- NA
survey_2022[survey_2022$sp == "glab", "sp"] <- "gla"
unique(survey_2022$status)

#combine survey and sample
colnames(survey_2022) <- c("group", "sp", "id", "mtag", "tran", "dist","trt","tag","side","Notes","loc", "value","status")
data_2022 <- left_join(survey_2022, sample_2022, by = c("group", "sp","id","tran","dist","trt","loc"))

# add year columns
data_2022$year <- 2022

data_2022


