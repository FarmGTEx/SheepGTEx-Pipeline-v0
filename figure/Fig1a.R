library(ggplot2)
library(openxlsx)
library(maps)

setwd("D:/博士/pic(PDF&TIFF)/breed_map")
mydata <- read.xlsx("breed_location_20241113.xlsx", sheet = "location")
mydata$lognumber <- log(mydata$numbers)
mydata$size_mapped <- mydata$lognumber
summary(mydata$size_mapped)

# 检查数据中的 Species 列的唯一值
print(unique(mydata$Species))

# 设置地图数据的范围, 经度需要横跨360度
mapWorld <- map_data('world', wrap = c(-180, 180), ylim = c(-90, 90))
mydata <- na.omit(mydata)

# 定义颜色映射
colors <- c(
  "RNA-seq"="#1D5EAD",
  "paired RNA-seq and WGS"="#F1AA00")

# 定义填充颜色映射（更浅的颜色）
fill_colors <- c(
  "RNA-seq"="#A3B9D9",
  "paired RNA-seq and WGS"="#F3D9A5"
)

# 分别筛选出 Sheep 和 Goat 的数据
sheep_data <- mydata[mydata$Species == "Sheep", ]
goat_data <- mydata[mydata$Species == "Goat", ]

# 绘制 Sheep 的地图
mp_sheep <- ggplot() +
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group),
               fill = "#ffffff",
               color = "#ffffff",
               size = 0.5
  ) +
  geom_point(data = sheep_data[sheep_data$Type == "Public", ], aes(x = Longitude, y = Latitude, 
                                                                   fill = classification, shape = Type,
                                                                   color = classification, size = lognumber), 
             alpha=0.8,stroke = 1.5  # 修改点的外观
  ) +
  geom_point(data = sheep_data[sheep_data$Type == "This study", ], aes(x = Longitude, y = Latitude, 
                                                                       fill = classification, shape = Type,
                                                                       color = classification, size = lognumber), 
             alpha=0.6,stroke = 1.5  # 修改点的外观
  ) + 
  scale_radius(range = c(1, 13),
               breaks = log(c(50, 100, 500, 1000, 2000)),
               labels = c("50", "100", "500", "1000", "2000")
  ) + 
  scale_shape_manual(values = c("This study" = 24, "Public" = 21)) +  # 使用不同的形状
  scale_color_manual(values = colors) +  # 使用更深的颜色作为边框
  scale_fill_manual(values = fill_colors) +  # 使用更浅的颜色作为填充
  theme(panel.background = element_rect(alpha('#457EAC',0.3)),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        axis.text = element_blank(),
        legend.text = element_text(size = 14)
  ) + labs(x = "", y = "") +
  coord_cartesian(xlim = c(-120, 174), ylim = c(-43, 66)) + 
  theme(legend.position = c(0.99, 0.65),
        legend.justification = c(0.95, 0.5), 
        legend.background = element_rect(fill = rgb(1, 1, 1, alpha = 0.001), colour = NA),
        legend.key = element_rect(fill = "transparent"),
        legend.spacing.y = unit(0.1, 'cm')) +
  guides(
    shape = guide_legend(title = NULL, override.aes = list(size = 4, stroke = 1.5)),
    size = guide_legend(title = NULL),
    color = guide_legend(title = NULL, override.aes = list(size = 4, stroke = 1.5)),
    fill = guide_legend(title = NULL, override.aes = list(size = 4, stroke = 1.5))
  ) 

# 打印和保存图形
print(mp_sheep)
ggsave('breed_map_sheep_11.13.pdf', plot = mp_sheep, width = 15, height = 7.5)

#print(mp_goat)
#ggsave('breed_map_goat_11.7.pdf', plot = mp_goat, width = 15, height = 7.5)