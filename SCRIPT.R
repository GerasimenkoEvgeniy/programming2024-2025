# libraries ---------------------------------------------------------------
library(tmap)#пакет для работы с картой tmap
library(sf)#для работы с sf объектами
library(leaflet)#для интерактивной карты
library(raster)#для интерактивной карты
library(tidyverse)
library(vegan)
library(iNEXT)
library(ggspatial)

# ---- Юманова Екатерина ----

curl::curl_download(
    "https://api.gbif.org/v1/occurrence/download/request/0002236-250402121839773.zip", 
    "0002236-250402121839773.zip")
d <- readr::read_delim(unzip("0002236-250402121839773.zip", "occurrence.txt")) %>% 
  select(family, genus, specificEpithet, scientificName, decimalLongitude,decimalLatitude,
         eventDate, day, month, year, parentEventID, habitat,locality, individualCount)%>%
  rename(line = parentEventID) %>% 
  filter(line != !is.na(line), str_detect(line, "PZP"))
# Выудили данные только по многолетним (более 4 лет) площадкам

#Находим число особей каждого вида в разные года на разных площадках
d1 <- d %>% 
  group_by(scientificName, year, line) %>% 
  summarise(totalCount = sum(individualCount), .groups = "drop")
#Находим число видов за разные года на разных площадках
d2 <- d %>% 
  filter(!line %in% c("Line 1", "Line 2","Line 3", "Line 4")) %>% 
  group_by(line, year, scientificName) %>%
  slice(1) %>% 
  summarise(.groups = "drop_last") %>% 
  summarise(speciesCount=n(), .groups = "drop_last") %>% 
  mutate(max = max(speciesCount), min = min(speciesCount), mean = mean(speciesCount)) %>% 
  ungroup()
# базовая статистика
d3 <- d2 %>% 
  select(-year, -speciesCount) %>% 
  group_by(line) %>% 
  slice(1) %>% 
  ungroup()

# Измайлов Е.; Различные расчеты ------------------------------------------
# подсчет менхинника и шеннона 
Indexes <-  d1 %>% 
  group_by(year, line, scientificName) %>%
  summarise(individualCount = sum(totalCount), .groups = "drop_last") %>%
  summarise(
    N = sum(individualCount), 
    S = n(), 
    iMn = S / sqrt(N),
    iH = diversity(
      individualCount, 
      index = "shannon"), 
    .groups = "drop")

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

#Дополнение таблиц экологией видов
v <- rio::import(
  "https://raw.githubusercontent.com/ANSozontov/Revda_2024/main/clean_data12.xlsx", 
  which = "Arachnida.Traits")%>% 
  select(taxa, lifeform, storey__) %>% 
  rename(storey1 = storey__, lifeform1 = lifeform) %>% 
  mutate(taxa = str_remove(taxa, "-[fm]$")) %>% 
  distinct()

sp <- rio::import(
    "https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/species_traits.xlsx", 
    which = "species")%>%
  rename(taxa = Species, storey2 = storey, lifeform2 = lifeform) %>% 
  select(taxa, storey2, lifeform2)

v1 <- d1
v1$scientificName <- str_replace(v1$scientificName, "\\s*\\(.*\\)", "")
v1 <- v1 %>% 
  rename(taxa = scientificName)

v2 <- left_join(v1, v, by = c("taxa"))
v3 <- left_join(v2, sp, by = c("taxa")) %>% 
  mutate(lifeform = coalesce(lifeform1, lifeform2),
         storey = coalesce(storey1, storey2)
        ) %>% 
  select(-lifeform1, -lifeform2, -storey1, -storey2)
# сколько всего видов обнаружено за всё время на площадках
v3 %>%
  group_by(line) %>%
  summarise(unique_taxa = n_distinct(taxa)) 
# процент видов, имеющих экологическую характеристику
v3 %>%
    summarise(
        total_rows = n(),
        complete_rows = sum(complete.cases(.)),
        percent_complete = (complete_rows / total_rows) * 100
    )

# #Карта Висима -----------------------------------------------------------
#Считывание файлов со слоями и предварительная обработка данных
granitsy<-read_sf('https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/visim_mini.gpkg', layer='ural_oopt')
rusla<-read_sf('https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/visim_mini.gpkg', layer='reki_rusla')
pritoki<-read_sf('https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/visim_mini.gpkg', layer='reki_pritoki')
iso_poly<-st_read('https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/visim_mini.gpkg', layer='iso_poly')#тут рельеф, взят из карты 5 занятия
occur<-d%>%   
dplyr::select(scientificName,individualCount,line,decimalLongitude,decimalLatitude)%>%#выбираем нужные колонки
  as.data.frame()
occur[occur == ""] <- "none"  # Заменяем пустые строки на слово none
occure<-filter(occur,!line %in% c("Line 1", "Line 2","Line 3", "Line 4")) %>% filter(line!='none') 

#Создание точек площадок
platform<-filter(occure,line!='none')%>%
  rename(lon=decimalLongitude,lat=decimalLatitude,species=scientificName) %>% 
  st_as_sf(coords=c('lon','lat'),crs=4326)#создаём точки для площадок

#Таблица для точек площадок и регулирования их размера
platformp<-platform%>%
  group_by(line)%>%
  summarise(unique_species = n_distinct(species))#подсчёт количества уникальных видов для каждой площадки
platformp$index <- row.names(platformp)

#Создание точки для заповедника на карте России и рамки
data('World')
world=st_union(World)
a<-granitsy%>%
  st_bbox()
x<-mean(a[c(1,3)])
y<-mean(a[c(2,4)])
cord<-data.frame(col.names = c('lon'=x,'lat'=y))
point <- st_point(c(y, x))
point_sf <- st_sf(geometry = st_sfc(point), crs = 4326)
ramka<-st_bbox(c(xmin = 10, xmax = 150, ymax = 85, ymin = 10), crs = st_crs(4326))#рамка(обрезка)карты(1.создание рамки, как тут 2.указываем в тм_шейп в ббоксе)

#Создание карты с точкой
tmap_options(check.and.fix = TRUE)
t<-tm_shape(world,bbox=ramka)+#использование объекта
  tm_borders(lwd = 0.5) +
  tm_fill(col='white')+
  tm_layout(bg.color = 'skyblue')+
  tm_shape(point_sf)+
  tm_dots(col='red',size = 1)
tmap_save(t,'точка.png',dpi=300,width = 6, height = 3)

#Карта заповедника в ggplot
ggplot() +
  geom_sf(data = iso_poly, aes(fill = ELEV, color = ELEV)) +
  scale_color_gradient2(midpoint = 500, low = "lightblue", mid = "lightyellow", high = "pink") +
  geom_sf(data = granitsy, color = "black", fill = NA, linewidth = 0.5) +
  scale_fill_gradient2(midpoint = 500, low = "lightblue", mid = "lightyellow", high = "pink") +
  geom_sf(data = rusla, color = "blue", linewidth = 0.5) +
  geom_sf(data = pritoki, color = "blue", linewidth = 0.3) +
  theme_bw() +
  geom_sf(data = platformp, color = "black", fill = "red", size = platformp$unique_species / 30, shape = 21) +  
  geom_sf_text(data = platformp,aes(label = index),nudge_x = 0.001,nudge_y = 0.004,color = "black",size = 4,check_overlap = FALSE) +
  labs(title = "Висимский заповедник") +
  theme(legend.position = 'none', axis.title = element_blank())+
  annotation_north_arrow(location = "tl",style = north_arrow_nautical() ) 
ggsave('картав.png',dpi=300,width = 6, height = 4)

# Интерактивная карта -----------------------------------------------------
r <- raster("https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/visim_elevation.tiff")
terraincolors <- colorNumeric(
  palette = terrain.colors(10),
  domain = values(r),
  na.color = "transparent")
const<-read.delim('occurrence.txt')%>%
  dplyr::select(scientificName,individualCount,parentEventID,decimalLongitude,decimalLatitude)%>%#выбираем нужные колонки
  rename(species=scientificName,count=individualCount,lon=decimalLongitude,lat=decimalLatitude,line = parentEventID)%>%#переименовываем
  mutate(comment = case_when(
    line == "PZP-02" ~ "PZP-02.<br>Заложена в 1985.\nПапоротниково-злаковый пихтово-ельниковый лес.\nПочвы горно-лесные бурые, средне- или тяжелосуглинистые.\nПожар в 2010 году.",
    line == "PZP-07" ~ "PZP-07.<br>Заложена в 1991.\nВейниково-злаковый березовый лес.\nПочвы бурые горно-лесные, средне- или тяжелосуглинистые.",
    line == "PZP-19" ~ "PZP-19.<br>Заложена в 1996.\nПапоротниково-злаковый лес.\nПочвы горно-лесные бурые, средне- или тяжелосуглинистые.",
    line == "PZP-20" ~ "PZP-20.<br>Заложена в 1998.\nПочвы здесь бурые горно-лесные, средне- и тяжелосуглинистые.\nВ 1998 году прошел сильный низовой пожар.\nПозднее в 2010 году произошел слабый низовой пожар.\nСостав растительности сильно менялся во времени."))%>%
  filter(comment != !is.na(comment)) %>%  # оставляем данные только многолетних площадок
  distinct(line, comment, lon, lat) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
leaflet() %>% 
  addRasterImage(r, colors = terraincolors) %>% 
  addPolygons(data = granitsy, color="black", opacity = 1, weight = 2, fillOpacity = 0) %>%
  addPolylines(data=pritoki, color="blue", opacity = 1, weight = 1) %>% 
  addPolylines(data=st_zm(rusla), color="blue", opacity = 1, weight = 3) %>% 
  addCircleMarkers(data=const, popup = ~ comment, label = ~line, radius = 5, weight = 1, color="black", fillColor = "red", opacity = 1, fillOpacity = 1)

# Корепанова М. Графики ---------------------------------------------------
#График расчётов
fires<-rio::import('https://raw.githubusercontent.com/GerasimenkoEvgeniy/programming2024-2025/main/Forest_fires.xlsx')
res1_long<-res1%>%
  pivot_longer(cols = c(totalCount, speciesCount, species400),names_to = "variable",values_to = "count")
res1_long %>%
  mutate(y = count, year = as.numeric(year)) %>% 
  filter(variable != 'totalCount') %>%
  ggplot( aes(x = year, y = y, color = variable)) + 
  geom_vline(
    data = fires,
    mapping = aes(xintercept = year),
    color = "darkred",
    linewidth = 1) +
  geom_line()+
  scale_color_manual(values = c(
    "speciesCount" = "green", "species400" = "orange"))+
  facet_wrap(~line,nrow=2,ncol=2)+
  scale_x_continuous(breaks = seq(1985, 2018, by = 5)) +
  theme_bw() + 
  labs(y='Число видов', x= 'Год')+
  theme(axis.text.x = element_text(angle=90), 
        legend.position = 'none')
ggsave('график_расчётов.png',dpi=300, width = 6, height = 6)

#График числа особей
res1_long %>%
  mutate(y = count, year = as.numeric(year)) %>% 
  filter(variable == 'totalCount') %>%
  ggplot( aes(x = year, y = y)) + 
  geom_vline(
    data = fires,
    mapping = aes(xintercept = year),
    color = "darkred",
    linewidth = 1) +
  geom_line(col = 'blue')+
  facet_wrap(~line,nrow=2,ncol=2)+
  scale_x_continuous(breaks = seq(1985, 2018, by = 5)) +
  theme_bw() + 
  labs(y='Число особей', x= 'Год')+
  theme(axis.text.x = element_text(angle=90), 
        legend.position = 'none')
ggsave('график_числа_особей.png',dpi=300, width = 6, height = 6)

#График по Менхинику и Шеннону
indexes<-Indexes%>%
  dplyr::select(-N,-S) %>% 
  pivot_longer(cols = c(iMn,iH),names_to = "index",values_to = "count")
indexes %>%
  mutate(y = count, year = as.numeric(year)) %>% 
  ggplot( aes(x = year, y = y, color = index)) + 
  geom_vline(
    data = fires,
    mapping = aes(xintercept = year),
    color = "darkred",
    linewidth = 1) +
  geom_line()+
  scale_color_manual(values = c('iH'='blue','iMn'='magenta'))+
  facet_wrap(~line,nrow=2,ncol=2)+
  scale_x_continuous(breaks = seq(1985, 2018, by = 5)) +
  theme_bw() + 
  labs(y='Значение', x= 'Год')+
  theme(axis.text.x = element_text(angle=90),legend.position = 'none')
ggsave('график_индексов.png',dpi=300, width = 6, height = 6)

#Графики по экологии
v3 <- v3 %>% 
  filter(!is.na(storey)) %>% 
  filter(!is.na(lifeform))
v3.1<-group_by(v3,year,storey,line) %>% summarise(sum=sum(totalCount))
v3.2<-group_by(v3,year,lifeform,line) %>% summarise(sum=sum(totalCount))
#storey
ggplot()+
  geom_bar(data=v3.1, aes(x=year, y=sum, fill=storey), stat = "identity",position = 'fill', width = 0.9)+
  geom_vline(data = fires,mapping = aes(xintercept = year),color = "darkred",linewidth = 1) +
  scale_fill_hue(labels = c("Герпетобионты", "Хортобионты","Стратобионты"))+
  labs(y='Доля особей',x="Год",fill='Ярус')+
  scale_x_continuous(breaks = seq(1985, 2018, by = 5))+
  theme_bw()+
  facet_wrap(~line)
ggsave("storey.png",dpi=300,width = 9, height = 6)
#lifeform
ggplot()+
  geom_bar(data=v3.2, aes(x=year, y=sum,fill=lifeform), stat = "identity",position = 'fill', width = 0.9)+
  geom_vline(data = fires,mapping = aes(xintercept = year),color = "darkred",linewidth = 1) +
  scale_fill_hue(labels = c("Охотники", "Тенетники","Охотники-тенетники"))+
  labs(y='Доля особей',x="Год",fill='Жизненная форма')+
  scale_x_continuous(breaks = seq(1985, 2018, by = 5))+
  theme_bw()+
  facet_wrap(~line)
ggsave("lifeform.png",dpi=300,width = 9, height = 6)