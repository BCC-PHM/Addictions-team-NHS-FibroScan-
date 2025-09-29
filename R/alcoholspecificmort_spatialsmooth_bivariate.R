setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/addiction team/Jad_liverbus")



library(biscale)
library(cowplot)
library(ggspatial)

library("tidyverse")
library("dplyr")
library("sf")
library("sp")
library("spdep")
library("INLAspacetime")
library(INLA)
library(fmesher)
library(readr)
library(data.table)

alcohol_specific_mortality <- read_excel("alcohol_specific_mortality.xlsx")


specificmort = alcohol_specific_mortality %>% 
  select(DEC_AGEC, REG_DATE, WARD_OF_RESIDENCE_CODE) %>%
  mutate(
    AdmissionDate = as.Date(REG_DATE, format = "%m-%d-%Y"),
    AgeGroup = case_when(
      DEC_AGEC == 0 ~ "UNDER 1",
      DEC_AGEC >= 1 & DEC_AGEC <= 4 ~ "1-4",
      DEC_AGEC >= 5 & DEC_AGEC <= 9 ~ "5-9",
      DEC_AGEC >= 10 & DEC_AGEC <= 14 ~ "10-14",
      DEC_AGEC >= 15 & DEC_AGEC <= 19 ~ "15-19",
      DEC_AGEC >= 20 & DEC_AGEC <= 24 ~ "20-24",
      DEC_AGEC >= 25 & DEC_AGEC <= 29 ~ "25-29",
      DEC_AGEC >= 30 & DEC_AGEC <= 34 ~ "30-34",
      DEC_AGEC >= 35 & DEC_AGEC <= 39 ~ "35-39",
      DEC_AGEC >= 40 & DEC_AGEC <= 44 ~ "40-44",
      DEC_AGEC >= 45 & DEC_AGEC <= 49 ~ "45-49",
      DEC_AGEC >= 50 & DEC_AGEC <= 54 ~ "50-54",
      DEC_AGEC >= 55 & DEC_AGEC <= 59 ~ "55-59",
      DEC_AGEC >= 60 & DEC_AGEC <= 64 ~ "60-64",
      DEC_AGEC >= 65 & DEC_AGEC <= 69 ~ "65-69",
      DEC_AGEC >= 70 & DEC_AGEC <= 74 ~ "70-74",
      DEC_AGEC >= 75 & DEC_AGEC <= 79 ~ "75-79",
      DEC_AGEC >= 80 & DEC_AGEC <= 84 ~ "80-84",
      DEC_AGEC >= 85 & DEC_AGEC <= 89 ~ "85-89",
      DEC_AGEC >= 90 ~ "90+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(AdmissionDate >= "2024-04-01" & AdmissionDate <="2025-3-31" ) %>% 
  group_by(WARD_OF_RESIDENCE_CODE, AgeGroup) %>%
  summarise(count = n(), .groups = 'drop') 




#just use the latest estimate 
pop_estimate <- read_excel("~/R projects/addiction team/BDAP/pop_estimate_by_ward_2023_2018.xlsx")

pop_estimate = pop_estimate %>% 
  filter(Year == "2023") %>% 
  select(-Year)



#merge population estimate to mortaluty data 


mergedward = merge(
  pop_estimate,
  specificmort ,
  all.x = TRUE,
  by.x = c("Ward", "AgeGroup"),
  by.y = c("WARD_OF_RESIDENCE_CODE", "AgeGroup")
)


############################################################################

#1) Build adjacency graph with poly2nb

ward_map <- st_read("Birmingham Ward map.json")
ward_map <- st_set_crs(ward_map, 4326)

#Create new IDs in order of data
ward_map$new_id = 1:nrow(ward_map)


#Create the graph file to indicate who are and are not neighbors
#Reproject to British National Grid (EPSG:27700)
ward_map <- st_transform(ward_map, 27700)  # meters, avoids angular distortion
ward_map <- st_make_valid(ward_map)
which(!st_is_valid(ward_map))


# Export to INLA graph
mcnty_nb = poly2nb(ward_map, row.names = ward_map$new_id, queen = F)

nb2INLA("ward.adj", mcnty_nb )
g = inla.read.graph("ward.adj")




# create index for ward
ward_levels = ward_map %>%
  st_drop_geometry() %>%
  pull(new_id)  

lookup <- ward_map %>%
  st_drop_geometry() %>%
  transmute(Ward = as.character(WD21CD),
            new_id = as.character(new_id))

mergedward <- mergedward %>%
  mutate(Ward = as.character(Ward)) %>%
  left_join(lookup, by = "Ward")


# Build the INLA indices for age group
age_levels = c(
  "UNDER 1","1-4","5-9","10-14","15-19","20-24","25-29","30-34",
  "35-39","40-44","45-49","50-54","55-59","60-64","65-69",
  "70-74","75-79","80-84","85-89","90+"
)



inla_dat  = mergedward %>%
  mutate(
    AgeGroup = as.numeric(factor(AgeGroup, levels = age_levels)),
    new_id   = as.numeric(factor(new_id, levels = ward_levels)),  # align to map order
    Ward_index = as.integer(new_id),                  # 1..N matching graph nodes
    Age_index  = as.integer(AgeGroup),
    count      = ifelse(is.na(count), 0L, as.integer(count)),
    local_pop  = ifelse(is.na(local_pop), 0L, as.integer(local_pop))
  )




# make sure local pop is numeric
inla_dat = inla_dat %>% mutate(local_pop = as.numeric(local_pop))




#########################################################################
#paper that i follow
# https://ij-healthgeographics.biomedcentral.com/articles/10.1186/s12942-020-00251-z

#define the fucntion
formula2<- count ~ 1 +     #intercept
  f(Age_index, model = "ar1", constr = TRUE,  #assumption of first order autoregressive age dependence structure
    hyper = list(
      theta1 = list(prior = "pc.prec", param = c(1, 0.01)),  # standard deviation is probably small (keeps the effect smooth),
      theta2  = list(prior = "pc.cor1", param = c(0, 0.9))    # The correlation is likely to be strong
    ))+      
  f(Ward_index,                                 
    model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01)),
      phi  = list(prior = "pc",      param = c(0.5, 2/3))
    ),
    group = Age_index,
    control.group = list(model = "ar1")
  )


#run the function using poisson
fit2 = inla(
  formula2, 
  family = "poisson",
  offset = log(local_pop),
  data = inla_dat,
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  verbose = TRUE
)


#####################################
# look at the summary!
summary(fit2) 
head(fit2$summary.fitted.values)



#bind them to the smoothed rate estimates you just calculated
inla_result = inla_dat %>%
  cbind(fit2$summary.fitted.values[, c("mean","0.025quant","0.975quant")]) %>%
  rename(smoothedcount= mean, lower025 = `0.025quant`, upper975 = `0.975quant`) %>% 
  mutate(age_standardised_count = (smoothedcount/local_pop)* Standard_Population ) %>% 
  group_by(Ward ) %>% 
  summarise(
    local_event = sum(count),
    local_pop = sum(local_pop),
    Standard_Population = 100000,
    age_standardised_count = sum(age_standardised_count),
    .groups = 'drop'
  )




########################################################

ward_map_inla <- st_read("Birmingham Ward map.json")

ward_map_inla <- st_set_crs(ward_map_inla , 4326)

#########################################################
#join inla result
ward_map_inla = ward_map_inla %>% 
  left_join(inla_result, by = c("WD21CD" = "Ward"))


#load imd data and process
imd_ward<- read_csv("imd-income-deprivation-score-percentage-birmingham-wards.csv")


imd_ward = imd_ward %>% 
  filter(`Period Label` == "2019") %>% 
  rename(IMD = Value)

colnames(imd_ward) = gsub(" ", "_", colnames(imd_ward))

#join imd data 
warddatainla = ward_map_inla  %>% 
  left_join(imd_ward, by = c("WD21CD" = "Area_Code")) 



#apply bivaraiate function
smootheddata = bi_class(warddatainla , x = age_standardised_count, y = IMD, style = "quantile", dim = 4)


#################################################################
#plot the map 
smoothedplot = ggplot() +
  geom_sf(data = smootheddata, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  geom_sf_label(data = smootheddata %>% 
                  filter(bi_class == "4-4"), aes(label = WD21NM), size = 2.8) + 
  bi_scale_fill(pal = "GrPink2", dim = 4) +
  labs(
    title = "Smoothed alcohol-specific mortality in 24/25 \nand IMD score",
    caption = "Contains OS data ©️ Crown copyright and database right 2025, Source: \nOffice for National Statistics licensed under the Open Government Licence v.3.0"
  ) +
  bi_theme(base_size = 14)+
  theme( plot.caption.position = "plot",
         plot.caption = element_text(hjust = 0, size = 10),
         plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x  = element_blank(),
         axis.text.y  = element_blank(),
         axis.ticks   = element_blank())+
  ggspatial::annotation_north_arrow(
    location = "br", which_north = "true",
    height = grid::unit(3.2, "cm"),     
    width  = grid::unit(3.2, "cm"),     
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("black", "white"),
      line_col = "grey20",
      text_family = "ArcherPro Book"
    )
  )

#create label
labels1  = bi_class_breaks(warddatainla, x = age_standardised_count, y = IMD, style = "quantile", 
                           dim = 4, dig_lab = 4, split = FALSE)

#label squeezed together so i add\n
labels1$bi_x = sub("-", "\n–", labels1$bi_x) 
#Creating Legends

legend <- bi_legend(pal = "GrPink2",
                    dim = 4,
                    xlab = "Higher mortality rate \nper 100,000",
                    ylab = "Higher deprivation %",
                    breaks = labels1,
                    size = 11)



# combine map with legend
alcoholspecilmort_sARs = ggdraw() +
  draw_plot(smoothedplot, 0, 0, 1, 1) +
  draw_plot(legend, -0.01, .6, 0.35, 0.35)

#can only save in view w=792 h=829
ggsave("alcoholspecilmort_sARs.png", alcoholspecilmort_sARs, dpi=600)



###########################################################################################












