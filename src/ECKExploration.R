library(tidyverse)
library(lubridate)
library(vegan)

load("species-abundance_2022-04-07.Rdata")

#1. Repeat visits assessment----
locs <- species.long %>% 
  dplyr::select(location, recording_date, recording_time, latitude, longitude) %>% 
  unique() %>% 
  mutate(date = ymd_hms(paste0(recording_date, recording_time)),
         doy = yday(date),
         hour = hour(date)) %>% 
  group_by(location, latitude, longitude) %>% 
  mutate(n = n()) %>% 
  ungroup()

table(locs$n)
ggplot(locs) +
  geom_hex(aes(x=doy, y=hour))

#2. Looking at data for amphibs of interest----
am <- species.long %>% 
  dplyr::filter(species_code %in% c("CATO", "WETO", "PLSP", "GPTO")) %>% 
  dplyr::select(location, recording_date, recording_time, latitude, longitude, species_code) %>%
  group_by(species_code, location, latitude, longitude) %>% 
  mutate(n.am = n()) %>% 
  ungroup() %>% 
  right_join(locs) %>% 
  mutate(species_code = ifelse(is.na(species_code), "none", species_code))

table(am$n.am, am$species_code)
ggplot(am) +
  geom_jitter(aes(x=doy, y=hour)) +
  facet_wrap(~species_code)

#3. Quick and dirty NMDS----
#Sample one survey per year
#Remove unknown species
set.seed(1234)
sp <- species.wide %>% 
  group_by(Location, Year) %>%
  sample_n(1) %>% 
  ungroup() %>% 
  dplyr::select(ALFL:TUSW, UPSA:YRWA) %>% 
  dplyr::select(-NONE)

#Remove rare species
sp.n <- colSums(sp)
sp <- sp[,which(sp.n > 5)]

#Remove rows==0
sp.0 <- rowSums(sp)
sp <- sp[which(sp.0 > 0),]

#Run (this takes ~ 10 hours)
nmds2 <- metaMDS(sp, k=2)
saveRDS(nmds2, "nmds2.rds")

#Wrangle for plotting
sp2 <- data.frame(nmds2$species) %>% 
  mutate(sp = rownames(nmds2$species),
         amphib = ifelse(sp %in% c("CATO", "WETO", "BCFR", "WOFR"), 1, 0)) %>% 
  dplyr::filter(!sp %in% c("NONE")) %>% 
  dplyr::filter(MDS1>0.001,
                MDS2>-.0005)

#Plot
ggplot(sp2) +
  geom_text(aes(x=MDS1, y=MDS2, label=sp, colour=factor(amphib)))
