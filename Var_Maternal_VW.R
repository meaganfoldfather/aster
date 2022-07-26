# Look at variability in maternal veg through time
source("Main_Script.R")

# dataframe with all maternal veg data
maternal

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


  

