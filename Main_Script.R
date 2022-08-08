# Build aster models to test for fitness effects of dispersal distance and maternal performance for the three focal Lasthenia species. 
# Will need to bring in fitness data for both experimental individuals and maternal plants for all three years of the experiment 

# libraries
library(aster)
library(tidyverse)
library(cowplot)
library(emmeans)

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
maternal_2019 <- read_csv("data/Jepson_2018_Maternal.csv")
maternal_2019$year <- 2019

###########################################################
# Bring in 2022 experimental data 
#experimental
source("scripts/Combining_2022SurveySample.R") # creates prelim data_2022 dataframe

#maternal 
maternal_2022 <- read_csv("data/Jepson_2020_Maternal.csv")
maternal_2022$year <- 2022

###########################################################
# Modify dataframes so they talk nicely to each other
data_2018 <- data_2018[ ,c("group","mayID", "sp", "trans", "dist", "trt","loc","year", "status", "vs")]
colnames(data_2018) <-c("group","id", "sp", "tran", "dist", "trt","loc","year", "status", "vs")
data_2018[data_2018$sp == "glab", "sp"] <- "gla"
data_2018 <-
    data_2018 %>% 
    filter(!is.na(status)) %>% 
    unite(col = "mom", remove = FALSE, sp, group, id) %>% 
    unite(col = "indiv",remove = FALSE, group, sp, id, tran, dist, trt, loc, year)
data_2018

data_2019 <- data_2019[ ,c("group","id", "sp", "tran", "dist", "trt","loc","year", "status", "vs")]
data_2019[data_2019$sp == "glab", "sp"] <- "gla"
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

# make number of viable seeds a zero if status is zero
all_data[all_data$status == 0,"vs"] <- 0
all_data <- 
all_data %>% 
    filter(!is.na(vs))

# reshape for aster - need it to be a long-form where there is a row for each event - for survival and for the number of seeds 
re_data<-
 all_data %>%
 pivot_longer(11:12, names_to = "vars", values_to = "resp") %>% 
  arrange(vars)

# change vars (status, vs) to a factor
re_data$vars <- as.factor(re_data$vars)
re_data

# change year to a factor 
re_data$year <- as.factor(re_data$year)
re_data

# add needed columns/vectors
re_data$root <- 1
pred <- c(0, 1)
fam <- c(1, 2); sapply(fam.default(), as.character)[fam]

# make a fitness vector
re_data$fit <- as.numeric(re_data$vars == "vs")

#Make a numeric 'id' row
nindvs_2018 <- 1:nrow(unique(re_data[re_data$year == 2018, "indiv"])) # 4221
nindvs_2019 <- 1:nrow(unique(re_data[re_data$year == 2019, "indiv"])) # 5734
nindvs_2022 <- 1:nrow(unique(re_data[re_data$year == 2022, "indiv"])) #3485
nindvs_2019 <- max(nindvs_2018) + nindvs_2019
nindvs_2022 <- max(nindvs_2019) + nindvs_2022
re_data$id <- c(nindvs_2018 , nindvs_2019,  nindvs_2022, nindvs_2018 , nindvs_2019,  nindvs_2022)
re_data$id  <- as.integer(re_data$id)

###########################################
# Aster models looking for an effect of distance, seperated by species
# Lasthenia californica model
aster_cal <- reaster(resp ~ vars + fit:(year+trt*dist),
                   list(group = ~ 0 + fit:group,
                        mom = ~ 0 + fit:mom),
                   pred, fam, vars, id, root, data = re_data[re_data$sp == "cal",])
summary(aster_cal, show.graph=TRUE)
# distance has a non-significant negative effect, NR has a significant positive effect, years 2019 and 2022 have significantly lower fitness (if include interaction between trt and dist. neither dist nor dist*trt is significant)

# Lasthenia fremontii model
aster_fre <- reaster(resp ~ vars + fit:(year + trt * dist),
                   list(group = ~ 0 + fit:group,
                        mom = ~ 0 + fit:mom),
                   pred, fam, vars, id, root, data = re_data[re_data$sp == "fre",])
summary(aster_fre, show.graph=TRUE)
# distance has a non-significant positive effect, NR has a significant positive effect, years 2019 and 2022 have significantly lower fitness (if include interaction between trt and dist. neither dist nor dist*trt is significant)

# Lasthenia glaberrima model
aster_gla <- reaster(resp ~ vars + fit:(year + trt * dist),
                   list(group = ~ 0 + fit:group,
                        mom = ~ 0 + fit:mom),
                   pred, fam, vars, id, root, data = re_data[re_data$sp == "gla",])
summary(aster_gla, show.graph=TRUE)
# distance has a non-significant negative effect, NR has a significant positive effect, no effect of year (if include interaction between trt and dist. the positive effect of NR is reduced by dispersing farther distances)

##################### 
# Modify maternal data so years match up 
maternal_2018
maternal_2018[maternal_2018$sp == "californica", "sp"] <- "cal"
maternal_2018[maternal_2018$sp == "fremontii", "sp"] <- "fre"
maternal_2018[maternal_2018$sp == "glaberrima", "sp"] <- "gla"
maternal_2018 <-
maternal_2018 %>% 
  filter(!is.na(group)) %>% 
  unite(col = "mom", sp, group, matID)
maternal_2018 <- maternal_2018[,c(1,3,7)]
colnames(maternal_2018) <- c("mom", "vw", "year")
maternal_2018


maternal_2019
maternal_2019[maternal_2019$sp == "CAL", "sp"] <- "cal"
maternal_2019[maternal_2019$sp == "FRE", "sp"] <- "fre"
maternal_2019[maternal_2019$sp == "GLA", "sp"] <- "gla"
maternal_2019<- 
maternal_2019 %>%
   unite(col = "mom", sp, group, mom)
maternal_2019 <- maternal_2019[,c(1,15,2)]
maternal_2019

plot(density(as.numeric(na.omit(maternal_2019$vw))))

maternal_2022
maternal_2022[maternal_2022$sp == "CAL", "sp"] <- "cal"
maternal_2022[maternal_2022$sp == "FRE", "sp"] <- "fre"
maternal_2022[maternal_2022$sp == "GLA", "sp"] <- "gla"
maternal_2022<- 
maternal_2022 %>%
   unite(col = "mom", sp, group, mom)
maternal_2022 <- maternal_2022[,c(1,15,2)]
maternal_2022

maternal <- rbind(maternal_2018,maternal_2019, maternal_2022)
maternal$sp <- substr(maternal$mom, start = 1, stop = 3)
maternal$year <- as.factor(maternal$year)
maternal %>% 
  ggplot(aes(year, vw, group = year))+
  geom_boxplot()+
  facet_wrap(.~sp)+
  theme_bw()

maternal <-
maternal %>% 
  group_by(year, sp) %>% 
  mutate(year_vw = round(mean(vw, na.rm = T),3))

# join fitness and maternal data together
fitness_maternal <- left_join(re_data, maternal, by = c("mom", "year", "sp"))
fitness_maternal

fitness_maternal <-
fitness_maternal %>% 
  filter(!is.na(vw))

#########################################################
# Aster models including maternal vw, seperated by species
# Lasthenia californica model
aster_cal_vw <- reaster(resp ~ vars + fit:(year + trt * vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])
summary(aster_cal_vw, show.graph=TRUE)
# positive effect of NR, fitness lower in 2019, 2022, positive effect of vw, no additive or interactive effect of dist (same if used average year vw) (if 3-way trt*vw*dist interaction included nothing beyond year, trt significant)

aster_cal_null <- reaster(resp ~ vars + fit:(trt + year + vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])
summary(aster_cal_null, show.graph=TRUE)
anova(aster_cal_null, aster_cal_vw) # interaction not supported with trt or vw with dist

# Lasthenia fremontii model
aster_fre_vw <- reaster(resp ~ vars + fit:(year + trt * vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])
summary(aster_fre_vw, show.graph=TRUE)
# positive effect of NR, fitness lower in 2019, 2022, significant interaction between vw and dist --> If your mom is big and you are far away, you have higher fitness--> if your mom was small and you are close by, you have higher fitness (this effect of dist is gone if use year averages of maternal vw) (if 3-way trt*vw*dist interaction included everything is significant! with both maternal vw as well)

aster_fre_null <- reaster(resp ~ vars + fit:(year + trt + vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])
summary(aster_fre_null, show.graph=TRUE)
anova(aster_fre_null, aster_fre_vw) # interaction supported

# Lasthenia glaberrima model
aster_gla_vw <- reaster(resp ~ vars + fit:(trt * vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])
summary(aster_gla_vw, show.graph=TRUE)
# positive effect of NR, no effect of year, positive effect of vw, non-significant positive effect of dist --> the positive effect of a big mom is washed out farther away; no direct effect of distance (this effect of dist is gone if use year averages of maternal vw)(3-way interaction model will not converge, unless year is removed) --> then just a positive effect of distance, confusing

aster_gla_null <- reaster(resp ~ vars + fit:(trt + vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])
summary(aster_gla_null, show.graph=TRUE)
anova(aster_gla_null, aster_gla_vw) # interaction supported

##############################################################
# Try to make predictions from the models - bases on no RE models - estimates are very different....
aster_cal_noRE <- aster(resp ~ vars + fit:(year + trt* vw * dist),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])

#aster_cal_noRE <- aster(resp ~ vars + fit:(year + vw * dist),
                   #pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal" & #fitness_maternal$trt == "C",])

summary(aster_cal_noRE) 
preds <- predict(aster_cal_noRE,varvar = vars, idvar = id, root = root,
    newdata = fitness_maternal[fitness_maternal$sp == "cal",], se.fit = TRUE)
# new_data_cal <- as.data.frame(expand.grid(vars = unique(fitness_maternal$vars), trt = unique(fitness_maternal$trt), year = unique(fitness_maternal$year), dist = unique(fitness_maternal$dist), vw = seq(min(fitness_maternal[fitness_maternal$sp == "cal", "vw"]), max(fitness_maternal[fitness_maternal$sp == "cal", "vw"]), by = .001)))
# new_data_cal$fit <- as.numeric(new_data_cal$vars == "vs")
#preds <- predict(aster_cal_noRE, varvar = vars, idvar = id, root = root,
#    newdata = new_data_cal, se.fit = TRUE)
#preds # how to best interpret the $fit values?

cal_data <- fitness_maternal[fitness_maternal$sp == "cal",]
cal_data$preds <- preds$fit 
cal_data$pred.se <- preds$se.fit 
cal_data

plot_cal <-
cal_data %>% 
  ggplot(aes(vw, preds, col = as.factor(dist))) + 
  geom_smooth(method = "lm", se = F) + 
  theme_classic()+
  scale_color_viridis_d()+
  xlab("Maternal Vegetation Weight")+
  ylab("Predicted Fitness")+
  ggtitle("L. cal")+
  facet_grid(.~trt)
plot_cal

#######
aster_fre_noRE <- aster(resp ~ vars + fit:(year + trt * vw * dist),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])

# aster_fre_noRE <- aster(resp ~ vars + fit:(year + vw * dist),
#                    pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre" & fitness_maternal$trt == "C",])

summary(aster_fre_noRE) 
preds <- predict(aster_fre_noRE,varvar = vars, idvar = id, root = root,
    newdata = fitness_maternal[fitness_maternal$sp == "fre",], se.fit = TRUE)

fre_data <- fitness_maternal[fitness_maternal$sp == "fre",]
fre_data$preds <- preds$fit 
fre_data$pred.se <- preds$se.fit 
fre_data

plot_fre<-
  fre_data %>% 
  ggplot(aes(vw, preds, col = as.factor(dist))) + 
  geom_smooth(method = "lm", se = F) + 
  theme_classic()+
  scale_color_viridis_d()+
  xlab("Maternal Vegetation Weight")+
  ylab("Predicted Fitness")+
  ggtitle("L. fre")+
  facet_grid(.~trt)
plot_fre


###########
aster_gla_noRE <- aster(resp ~ vars + fit:(trt * vw * dist),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])

# aster_gla_noRE <- aster(resp ~ vars + fit:(vw * dist),
#                    pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla" & fitness_maternal$trt == "C",])

summary(aster_gla_noRE) # issue with this model if year included - removed
preds <- predict(aster_gla_noRE,varvar = vars, idvar = id, root = root,
    newdata = fitness_maternal[fitness_maternal$sp == "gla",], se.fit = TRUE)

gla_data <- fitness_maternal[fitness_maternal$sp == "gla",]
gla_data$preds <- preds$fit 
gla_data$pred.se <- preds$se.fit 
gla_data

plot_gla <-
gla_data %>% 
  ggplot(aes(vw, preds, col = as.factor(dist))) + 
 geom_smooth(method = "lm", se = F) + 
  theme_classic()+
  scale_color_viridis_d()+
  xlab("Maternal Vegetation Weight")+
  ylab("Predicted Fitness")+
  ggtitle("L. gla")+
  facet_grid(.~trt)
plot_gla

plot_grid(plot_cal, plot_fre, plot_gla, nrow = 1)

###############
#Make Figures that have dispersal distance on the x-axis model fits split by large and small maternal size
# cal
cal_data
plot(density(cal_data$vw))
quantile(cal_data$vw)# median is 0.02

cal_data$size <- "Little"
cal_data[cal_data$vw > 0.02, "size"] <- "Big"

cal_data %>% 
  ggplot(aes(dist, preds, col = as.factor(size), fill = as.factor(size))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T, lty = "dashed") + 
  theme_classic()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle("L. cal")+
  facet_grid(.~trt)

cal_data %>% 
  ggplot(aes(dist, preds, col = as.factor(trt), fill = as.factor(trt))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = F, lty = "dashed") + 
  theme_classic()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle("L. cal")+
  facet_grid(.~size)

# fre
fre_data
plot(density(fre_data$vw))
quantile(fre_data$vw)# median is 0.023

fre_data$size <- "Little"
fre_data[fre_data$vw > 0.023, "size"] <- "Big"

fre_data %>% 
  ggplot(aes(dist, preds, col = as.factor(size), fill = as.factor(size))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T) + 
  theme_classic()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle("L. fre")+
  facet_grid(.~trt)

fre_data %>% 
  ggplot(aes(dist, preds, col = as.factor(trt), fill = as.factor(trt))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T) + 
  theme_classic()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle("L. fre")+
  facet_grid(.~size)
  
# gla
gla_data
plot(density(gla_data$vw))
quantile(gla_data$vw)# median is 0.029

gla_data$size <- "Little"
gla_data[gla_data$vw > 0.029, "size"] <- "Big"

gla_data %>% 
  ggplot(aes(dist, preds, col = as.factor(size),  fill = as.factor(size))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T,) + 
  theme_classic()+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle("L. gla")
  #facet_grid(.~trt)
  #ylim(0,1)
