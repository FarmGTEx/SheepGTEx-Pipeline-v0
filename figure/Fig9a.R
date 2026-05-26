library(ggplot2)
library(maps)
library(maptools)
library(data.table)
library(openxlsx)

#setwd("D:/西农")
mydata <- read.xlsx("C:\\Users\\Administrator\\Desktop\\sample.year.xlsx", sheet = "sample.site")
visit.x <- mydata$Longitude
visit.y <- mydata$Latitude

mydata$breed <- factor(mydata$Breed)
mydata$normalize_number <- sqrt(mydata$size) / sqrt(max(mydata$size))

# 设置地图数据的范围, 经度需要横跨360度
mapWorld <- map_data('world', wrap = c(-180, 180), ylim = c(-60, 60))

mp1 <- ggplot() +
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group),
               fill = "#D4D4D4",
               color = "#FFFAFA",
               size = 0.5
  ) +
  geom_point(data = mydata, aes(x = visit.x, y = visit.y), fill = mydata$color,
             shape = mydata$shape, color = mydata$color, alpha=0.9,size = mydata$normalize_number * 2
  ) +
  theme(panel.background = element_rect(fill = '#FFFAFA'),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "None",
        axis.text = element_blank()
  ) + labs(x = "", y = "") +
  # 设置地图范围
  coord_cartesian(xlim = c(-5, 150), ylim = c(25,60)) 

print(mp1)
