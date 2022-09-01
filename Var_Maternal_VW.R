# Look at variability in maternal veg through time
source("Main_Script.R")

library(lme4)

# dataframe with all maternal veg data
maternal

# Temporal variability in suitability relative the overall spatial variation 
# Is this the same for seed count and veg weight? Means different things for different species
# Negative/positive correlation in habitat values between years? 

# add in 2020 maternal collections to this analyses?

wide_maternal <- 
 maternal %>% 
pivot_wider(id_cols = mom, names_from = year, values_from = vw)

wide_maternal
wide_maternal$sp <- substr(wide_maternal$mom, 1, 3)
colnames(wide_maternal) <- c("mom","vw18", "vw19", "vw22", "sp")

wide_maternal %>% 
  ggplot(aes(vw18, vw19))+
  geom_point()+
  facet_wrap(.~sp)

wide_maternal %>% 
  ggplot(aes(vw18, vw22))+
  geom_point()+
  facet_wrap(.~sp)

wide_maternal %>% 
  ggplot(aes(vw19, vw22))+
  geom_point()+
  facet_wrap(.~sp)

# look at CV by species
maternal_analysis <- maternal
maternal_analysis$sp <- substr(maternal_analysis$mom, 1, 3)

maternal_analysis %>% 
  group_by(sp, mom) %>% 
  summarise(cv = sd(vw)) %>% 
  ggplot(aes(sp, cv))+
  geom_boxplot(outlier.size=0)+
  theme_classic()+
  ylim(0,0.1)

maternal_analysis %>% 
  ggplot(aes(year, vw, group = mom, color = sp))+ 
  geom_line()+
  facet_wrap(.~sp)

# Spatial variabilty 
#Within year variability for all moms within a site, seed set, veg weight? Are those two correlated?
# Make a site column
maternal_analysis$site <- substr(maternal_analysis$mom, 5, 5)
maternal_analysis

maternal_analysis %>% 
  group_by(year, sp, site) %>% 
  summarise(cv = sd(vw, na.rm =T)) %>% 
  ggplot(aes(sp, cv))+
  geom_boxplot(outlier.size=0)+
  theme_classic()+
  ylim(0,0.1)+
  facet_wrap(.~year)

# Temporal variabilty
# Same mom through time

# Hydrologic variability 
# soil moisture data; need some soil conversions for comparisons across species
# CV didn't appear different between species (from Courtney)

# Look at relationships between maternal vegetation weight and seed set and avg. seed weight --> need to go back to original year maternal datasets, only applicable for 2019 and 2022 maternal datasets
maternal_twoyears <- rbind(maternal_2019, maternal_2022)  
maternal_twoyears$sp <- substr(maternal_twoyears$mom, start = 1, stop = 3)
maternal_twoyears$year <- as.factor(maternal_twoyears$year)
maternal_twoyears$site <- substr(maternal_twoyears$mom, 5, 5)
maternal_twoyears

# Veg-weight ~ Seed Count
maternal_twoyears %>% 
  ggplot(aes(vw, vs))+
  facet_wrap(.~sp)+
  theme_classic()+
  geom_point()+
  geom_smooth(method = "lm", col = "black")
# Strong positive relationship, supported by the sp-specific models

summary(glmer(vs ~ vw + year + (1|site), family = "poisson", data = maternal_twoyears[maternal_twoyears$sp == "cal",]))
summary(glmer(vs ~ vw + year+ (1|site), family = "poisson", data = maternal_twoyears[maternal_twoyears$sp == "fre",]))
summary(glmer(vs ~ vw + year + (1|site), family = "poisson", data = maternal_twoyears[maternal_twoyears$sp == "gla",]))

# Veg-weight ~ Average seed weight (infl weight/counts of both viable and viable seeds)
maternal_twoyears$seed_w <- maternal_twoyears$iw/(maternal_twoyears$vs + maternal_twoyears$is)
maternal_twoyears %>% 
  ggplot(aes(vw, seed_w))+
  facet_wrap(.~sp)+
  theme_classic()+
  geom_point()+
  geom_smooth(method = "lm", col = "black")
# Positive relationship

summary(lmer(seed_w ~ vw + year + (1|site), data = maternal_twoyears[maternal_twoyears$sp == "cal",]))
summary(lmer(vs ~ vw + year+ (1|site), data = maternal_twoyears[maternal_twoyears$sp == "fre",]))
summary(lmer(vs ~ vw + year + (1|site), data = maternal_twoyears[maternal_twoyears$sp == "gla",]))

# Veg-weight ~ Veg height
maternal_twoyears %>% 
  ggplot(aes(vw, ih))+
  facet_wrap(.~sp)+
  theme_classic()+
  geom_point()+
  geom_smooth(method = "lm", col = "black")
# Strong positive relationship,

summary(lmer(ih ~ vw + year + (1|site), data = maternal_twoyears[maternal_twoyears$sp == "cal",]))
summary(lmer(ih ~ vw + year+ (1|site), data = maternal_twoyears[maternal_twoyears$sp == "fre",]))
summary(lmer(ih ~ vw + year + (1|site), data = maternal_twoyears[maternal_twoyears$sp == "gla",]))
