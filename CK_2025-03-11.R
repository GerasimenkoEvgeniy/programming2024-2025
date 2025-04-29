library(tidyverse)
library(ggplot2)
# Visim <- readr::read_delim("occurrence.txt")
Visim <- readr::read_delim(unzip("dwca-visim_spiders-v1.7.zip", "occurrence.txt"))
V <- Visim %>% 
    select(year, i = individualCount) %>%
    filter(!is.na(year)) %>% 
    mutate(i = case_when(is.na(i) ~ 1, TRUE ~ i)) %>% 
    group_by(year) %>%
    summarise(i_total = sum(i))

ggplot(V, aes(year,i_total))+
    geom_line()

Visim %>% 
    # mutate(
    #     habitat = str_replace(habitat, "\\(.*?\\)", ""), 
    #     habitat = str_replace(habitat, "with.*+", ""),
    #     habitat = str_squish(habitat)) %>% 
    group_by(habitat, scientificName) %>% 
    summarise() %>% 
    summarise(nsp = n()) %>% 
    arrange(desc(nsp))





    
    
    
