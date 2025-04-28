# ---- Юманова Е; начало! ----
library(tidyverse)
library(dplyr)

curl::curl_download(
    "https://api.gbif.org/v1/occurrence/download/request/0002236-250402121839773.zip", 
    "0002236-250402121839773.zip")
d <- readr::read_delim(unzip("0002236-250402121839773.zip", "occurrence.txt")) %>% 
    # readr::read_delim("dwca-visim_spiders-v1.7/occurrence.txt") %>%
  select(family, genus, specificEpithet, scientificName, 
         eventDate, day, month, year, parentEventID, habitat,locality, individualCount)%>%
  rename(line = parentEventID) %>% 
  filter(line != !is.na(line)) %>%  # оставляем данные только многолетних площадок
  filter(!line %in% c("Line 1", "Line 2"))

#Находим число особей каждого вида в разные года на разных площадках
d1 <- d %>% 
  group_by(scientificName, year, line) %>% 
  summarise(totalCount = sum(individualCount), .groups = "drop")
#Находим число видов за разные года на разных площадках
d2 <- d %>% 
  group_by(year, line) %>% 
  summarise(speciesCount=n(), .groups = "drop")





# Измайлов Е.; Различные расчеты ------------------------------------------
library(abdiv)
library(vegan)
library(iNEXT)

#тут подсчет менхинника и шеннона (просто скопировано с доработками)
Indexes <-  d1 %>% 
  group_by(year, line, scientificName) %>%
  summarise(individualCount = sum(totalCount), .groups = "drop_last") %>%
  summarise(
    N = sum(individualCount), 
    S = n(), 
    iMn = S / sqrt(N),
    iH = diversity(
      individualCount, 
      index = "shannon"))

#Тут перевод таблицы в формат для iNEXT
d1_wide <- d1 %>% 
  filter(totalCount > 0) %>% 
  unite("ID", year, line, sep = "_") %>% 
  pivot_wider(names_from = ID, values_from = totalCount, values_fill = 0) %>% 
  column_to_rownames("scientificName")

#iNEXT
res <- iNEXT(d1_wide, q = 0, se = FALSE, 
             datatype = "abundance", 
             size = c(400), nboot = 0)
res <- res$iNextEst$size_based %>% 
  select(-Order.q, -qD.LCL:-SC.UCL) %>% 
  mutate_if(is.numeric, ~round(.x, 1)) %>% 
  as_tibble() %>% 
  separate(Assemblage, into = c("year", "line"), sep = "_")

# ---- Юманова Е; продолжение! ----
res1 <- res %>% 
  group_by(year, line) %>% 
  mutate(totalCount = m[Method == "Observed"],
         speciesCount = qD[Method == "Observed"],
         species400 = qD[m == 400]
         ) %>% 
  slice(1) %>% 
  select(-m, -Method, -qD) %>% 
  ungroup()

#Дополнение таблицы экологией видов
v <- rio::import(
  "https://raw.githubusercontent.com/ANSozontov/Revda_2024/main/clean_data12.xlsx", 
  which = "Arachnida.Traits")%>% 
  select(taxa, lifeform, storey__) %>% 
  rename(storey=storey__) %>% 
  mutate(taxa = str_remove(taxa, "-[fm]$")) %>% 
  distinct()


v1 <- d1
v1$scientificName <- str_replace(v1$scientificName, "\\s*\\(.*\\)", "")
v1 <- v1 %>% 
  rename(taxa = scientificName)

v2 <- left_join(v1, v, by = c("taxa"))
