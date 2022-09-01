# Build aster models to test for fitness effects of dispersal distance and maternal performance for the three focal Lasthenia species. 
# Will need to bring in fitness data for both experimental individuals and maternal plants for all three years of the experiment 

# libraries
library(aster)
library(aster2)
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
maternal_2018_sub <- maternal_2018[,c(1,3,7)]
colnames(maternal_2018_sub) <- c("mom", "vw", "year")
maternal_2018_sub

maternal_2019
maternal_2019[maternal_2019$sp == "CAL", "sp"] <- "cal"
maternal_2019[maternal_2019$sp == "FRE", "sp"] <- "fre"
maternal_2019[maternal_2019$sp == "GLA", "sp"] <- "gla"
maternal_2019<- 
maternal_2019 %>%
   unite(col = "mom", sp, group, mom)
maternal_2019_sub <- maternal_2019[,c(1,15,2)]
maternal_2019_sub

plot(density(as.numeric(na.omit(maternal_2019_sub$vw))))

maternal_2022
maternal_2022[maternal_2022$sp == "CAL", "sp"] <- "cal"
maternal_2022[maternal_2022$sp == "FRE", "sp"] <- "fre"
maternal_2022[maternal_2022$sp == "GLA", "sp"] <- "gla"
maternal_2022<- 
maternal_2022 %>%
   unite(col = "mom", sp, group, mom)
maternal_2022_sub <- maternal_2022[,c(1,15,2)]
maternal_2022_sub

maternal <- rbind(maternal_2018_sub,maternal_2019_sub, maternal_2022_sub)
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

aster_cal_null_trt <- reaster(resp ~ vars + fit:(year + trt + vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])
summary(aster_cal_null_trt, show.graph=TRUE)

aster_cal_null_vw <- reaster(resp ~ vars + fit:(year + vw + trt * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])
summary(aster_cal_null_vw, show.graph=TRUE)

aster_cal_null_both <- reaster(resp ~ vars + fit:(year + vw*trt + trt * dist + vw*dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])
summary(aster_cal_null_both, show.graph=TRUE)


anova(aster_cal_null_trt, aster_cal_vw) # interaction not supported 
anova(aster_cal_null_vw, aster_cal_vw) # interaction not supported 
anova(aster_cal_null_vw, aster_cal_vw) # interaction not supported 

# updated_cal_model <- reaster(resp ~ vars + fit:(year + trt + vw + dist),
#                    list(group = ~ 0 + fit:group),
#                    pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "cal",])
# summary(updated_cal_model, show.graph=TRUE)

aout <- aster_cal_vw$obj
# aout <- updated_cal_model$obj
# create new data
preds <- predict.aster(aout, model.type = "conditional", is.always.parameter = F, parm.type = "mean.value", newdata = fitness_maternal[fitness_maternal$sp == "cal",])

cal_data <- fitness_maternal[fitness_maternal$sp == "cal",]
cal_data$preds <- preds

cal_data
plot(density(cal_data$vw))
quantile(cal_data$vw)# median is 0.02

#cal_data$size <- "Small Maternal Plant"
#cal_data[cal_data$vw > 0.02, "size"] <- "Large Maternal Plant"

cal_data$size <- "Average Maternal Plant"
cal_data[cal_data$vw > 0.034, "size"] <- "Large Maternal Plant"
cal_data[cal_data$vw < 0.013, "size"] <- "Small Maternal Plant"
cal_data$size <- as.factor(cal_data$size)
cal_data$size <- factor(cal_data$size, levels = c( "Small Maternal Plant","Average Maternal Plant", "Large Maternal Plant"))

cal_fitness_plot<-cal_data %>% 
  ggplot(aes(dist, preds, col = as.factor(size), fill = as.factor(size))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T, alpha = .3) + 
  theme_classic()+
  scale_color_viridis_d(option = "E")+
  scale_fill_viridis_d(option = "E")+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle(expression(italic("Lasthenia californica")))+
  facet_grid(.~trt)+
  theme(legend.title=element_blank(), text = element_text(size = 16))

# new_data_cal <- as.data.frame(expand.grid(vars = unique(fitness_maternal$vars), trt = unique(fitness_maternal$trt), year = unique(fitness_maternal$year), dist = unique(fitness_maternal$dist), vw = seq(min(fitness_maternal[fitness_maternal$sp == "cal", "vw"]), max(fitness_maternal[fitness_maternal$sp == "cal", "vw"]), by = .001)))
# new_data_cal$fit <- as.numeric(new_data_cal$vars == "vs")
# cal_data <- fitness_maternal[fitness_maternal$sp == "cal",]
# cal_data$preds <- preds$fit 
# cal_data$pred.se <- preds$se.fit 
# cal_data

# Lasthenia fremontii model
aster_fre_vw <- reaster(resp ~ vars + fit:(year + trt * vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])
summary(aster_fre_vw, show.graph=TRUE)
# positive effect of NR, fitness lower in 2019, 2022, significant interaction between vw and dist --> If your mom is big and you are far away, you have higher fitness--> if your mom was small and you are close by, you have higher fitness (this effect of dist is gone if use year averages of maternal vw) (if 3-way trt*vw*dist interaction included everything is significant! with both maternal vw as well)

aster_fre_null_trt <- reaster(resp ~ vars + fit:(year + trt + vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])
summary(aster_fre_null_trt, show.graph=TRUE)

aster_fre_null_vw <- reaster(resp ~ vars + fit:(year + vw + trt * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])
summary(aster_fre_null_vw, show.graph=TRUE)

aster_fre_null_both <- reaster(resp ~ vars + fit:(year + vw*trt + trt * dist + vw*dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "fre",])
summary(aster_fre_null_both, show.graph=TRUE)

anova(aster_fre_null_trt, aster_fre_vw) # interaction supported
anova(aster_fre_null_vw, aster_fre_vw) # interaction supported
anova(aster_fre_null_both, aster_fre_vw) # interaction supported

aout <- aster_fre_vw$obj
preds <- predict.aster(aout, model.type = "conditional", is.always.parameter = F, parm.type = "mean.value", newdata = fitness_maternal[fitness_maternal$sp == "fre",])

fre_data <- fitness_maternal[fitness_maternal$sp == "fre",]
fre_data$preds <- preds

# fre
fre_data
plot(density(fre_data$vw))
quantile(fre_data$vw)# median is 0.023

fre_data$size <- "Average Maternal Plant"
fre_data[fre_data$vw > 0.039, "size"] <- "Large Maternal Plant"
fre_data[fre_data$vw < 0.012, "size"] <- "Small Maternal Plant"
fre_data$size <- as.factor(fre_data$size)
fre_data$size <- factor(fre_data$size, levels = c( "Small Maternal Plant","Average Maternal Plant", "Large Maternal Plant"))

#fre_data$size <- "Small Maternal Plant"
#fre_data[fre_data$vw > 0.023, "size"] <- "Large Maternal Plant"

fre_fitness_plot <- fre_data %>% 
  ggplot(aes(dist, preds, col = as.factor(size), fill = as.factor(size))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T, alpha = 0.3) + 
  theme_classic()+
  scale_color_viridis_d(option = "E")+
  scale_fill_viridis_d(option = "E")+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle(expression(italic("Lasthenia fremontii")))+
  facet_grid(.~trt)+
  theme(legend.title=element_blank(), text = element_text(size = 16))

# Lasthenia glaberrima model
aster_gla_vw <- reaster(resp ~ vars + fit:(vw * trt * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])
summary(aster_gla_vw, show.graph=TRUE)
# positive effect of NR, no effect of year, positive effect of vw, non-significant positive effect of dist --> the positive effect of a big mom is washed out farther away; no direct effect of distance (this effect of dist is gone if use year averages of maternal vw)(3-way interaction model will not converge, unless year is removed) --> then just a positive effect of distance, confusing - always bad to disperse far away, except when it may not matter if very large in control plots

aster_gla_null_trt <- reaster(resp ~ vars + fit:(trt + vw * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])
summary(aster_gla_null_trt, show.graph=TRUE)

aster_gla_null_vw <- reaster(resp ~ vars + fit:(vw + trt * dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])
summary(aster_gla_null_vw, show.graph=TRUE)

aster_gla_null_both <- reaster(resp ~ vars + fit:(vw*trt + trt * dist + vw*dist),
                   list(group = ~ 0 + fit:group),
                   pred, fam, vars, id, root, data = fitness_maternal[fitness_maternal$sp == "gla",])
summary(aster_gla_null_both, show.graph=TRUE)

anova(aster_gla_null_trt, aster_gla_vw) # interaction supported
anova(aster_gla_null_vw, aster_gla_vw) # interaction supported
anova(aster_gla_null_both, aster_gla_vw) # 3-way interaction not supported

aout <- aster_gla_vw$obj
preds <- predict.aster(aout, model.type = "conditional", is.always.parameter = F, parm.type = "mean.value", newdata = fitness_maternal[fitness_maternal$sp == "gla",])

gla_data <- fitness_maternal[fitness_maternal$sp == "gla",]
gla_data$preds <- preds

# gla
gla_data
plot(density(gla_data$vw))
quantile(gla_data$vw)# median is 0.029

gla_data$size <- "Average Maternal Plant"
gla_data[gla_data$vw > 0.039, "size"] <- "Large Maternal Plant"
gla_data[gla_data$vw < 0.012, "size"] <- "Small Maternal Plant"
gla_data$size <- as.factor(gla_data$size)
gla_data$size <- factor(gla_data$size, levels = c( "Small Maternal Plant","Average Maternal Plant", "Large Maternal Plant"))

#gla_data$size <- "Small Maternal Plant"
#gla_data[gla_data$vw > 0.029, "size"] <- "Large Maternal Plant"

gla_fitness_plot <- gla_data %>% 
  ggplot(aes(dist, preds, col = as.factor(size),  fill = as.factor(size))) + 
  #geom_point()+
  geom_smooth(method = "lm", se = T, alpha = 0.3) + 
  theme_classic()+
  scale_color_viridis_d(option = "E")+
  scale_fill_viridis_d(option = "E")+
  xlab("Distance")+
  ylab("Predicted Fitness")+
  ggtitle(expression(italic("Lasthenia glaberrima")))+
  facet_grid(.~trt)+
  theme(legend.title=element_blank(), text = element_text(size = 16))

library(cowplot)
all <-plot_grid(cal_fitness_plot, fre_fitness_plot, gla_fitness_plot)

all_fitness_plot <- plot_grid(
  cal_fitness_plot + theme(legend.position="none"),
  fre_fitness_plot + theme(legend.position="none"),
  gla_fitness_plot + theme(legend.position="none"),
  get_legend(cal_fitness_plot + theme(legend.text=element_text(size=20))),
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 2
)
ggsave2(filename = "figures/fitness_preds_allSp.jpeg", all_fitness_plot)

