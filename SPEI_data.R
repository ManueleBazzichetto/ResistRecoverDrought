##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
##Extract SPEI data

library(sf)
library(raster)
library(elevatr)
library(ggplot2)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggspatial)
library(rcartocolor)

##----------------------------------------------------------------Compute region-specific centroids

#import kml of vegetation plots
Kml_obs.dir <- list.files(path = "~/Documents/ClimateExtremes/KML_obs/", pattern = ".kml", full.names = T)

Kml_obs <- lapply(Kml_obs.dir, st_read)

names(Kml_obs) <- c("A", "H", "S")

str(Kml_obs$A)

Kml_obs <- lapply(Kml_obs, function(.) {
  .$MyDescr <- ifelse(.$Description == "VIP - Forest" | .$Description == "EP - Forest", "Forest", "Grassland")
  .
})

plot(st_as_sfc(st_bbox(do.call(rbind, Kml_obs))))

invisible(lapply(Kml_obs, function(.) {
  plot(st_geometry(.[.$MyDescr == "Forest", ]), col = "blue", add = T)
  plot(st_geometry(.[.$MyDescr == "Grassland", ]), col = "green", add = T)
}))

#https://gis.stackexchange.com/questions/37219/recommended-coordinate-reference-system-crs-for-germany
lapply(Kml_obs, function(.) {
  Pl_obs <- st_transform(., crs = 25832)
  Cntrd <- st_centroid(st_union(Pl_obs))
  Cntrd.geo <- st_transform(Cntrd, crs = 4326)
})

#centroids for grasslands only! -> these will be used to extract SPEI values 
Kml_obs.gr <- lapply(Kml_obs, function(.) {
  df.gr <- .[which(.$MyDescr == "Grassland"), ]
  return(df.gr)
}) 

sapply(Kml_obs.gr, function(.) unique(.$MyDescr))

Centroids.gr <- lapply(Kml_obs.gr, function(.) {
  Pl_obs <- st_transform(., crs = 25832)
  Cntrd <- st_centroid(st_union(Pl_obs))
  Cntrd.geo <- st_transform(Cntrd, crs = 4326)
})

plot(st_as_sfc(st_bbox(do.call(rbind, Kml_obs.gr))))
plot(Centroids.gr$A, add = T, col = "blue")
plot(Centroids.gr$H, add = T, col = "purple")
plot(Centroids.gr$S, add = T, col = "red")

#https://gis.stackexchange.com/questions/251641/how-to-combine-sfc-objects-from-r-package-sf
Centroids.gr <- c(Centroids.gr$A, Centroids.gr$H, Centroids.gr$S)

plot(st_as_sfc(st_bbox(do.call(rbind, Kml_obs.gr))))
plot(Centroids.gr, add = T)

##----------------------------------------------------------------Map of the study area

#get DEM for the 3 regions
Bbox_ALB <- st_as_sf(data.frame(ID = seq_len(4), x = c(9.2, 9.2, 9.6, 9.6), y = c(48.35, 48.525, 48.35, 48.525)), coords = c("x", "y"))
Bbox_HAI <- st_as_sf(data.frame(ID = seq_len(4), x = c(10.12, 10.12, 10.8, 10.8), y = c(50.92, 51.36, 50.92, 51.36)), coords = c("x", "y"))
Bbox_SCH <- st_as_sf(data.frame(ID = seq_len(4), x = c(13.55, 13.55, 14.09, 14.09), y = c(52.83, 53.18, 52.83, 53.18)), coords = c("x", "y"))

st_crs(Bbox_ALB) <- st_crs(Bbox_HAI) <- st_crs(Bbox_SCH) <- 4326

DEM_ALB <- elevatr::get_elev_raster(locations = Bbox_ALB, z = 10, clip = "locations")
DEM_ALB <- as.data.frame(DEM_ALB, na.rm = T, xy = T)

DEM_HAI <- elevatr::get_elev_raster(locations = Bbox_HAI, z = 10, clip = "locations")
DEM_HAI <- as.data.frame(DEM_HAI, na.rm = T, xy = T)

DEM_SCH <- elevatr::get_elev_raster(locations = Bbox_SCH, z = 10, clip = "locations")
DEM_SCH <- as.data.frame(DEM_SCH, na.rm = T, xy = T)

#rename elevation column
colnames(DEM_ALB)[3] <- colnames(DEM_HAI)[3] <- colnames(DEM_SCH)[3] <- "Elevation"

#derive df of centroids of the 3 regions
Centr_df <- st_sf(Region = c("South-West", "Center", "North-East"), Centroids.gr)

Germany_border <- st_union(ne_states(country = "germany", returnclass = "sf"))

Germany_map <- ggplot(data = Germany_border) +
  geom_sf(fill = NA, col = "black", size = 1) +
  geom_sf(data = Centr_df, mapping = aes(col = Region), size = 6) +
  annotate(geom = "text", label = "South-West", x = 9.6, y = 48.75, col = "black", fontface = "bold", size = 9) +
  annotate(geom = "text", label = "Center", x = 10.4, y = 51.6, col = "black", fontface = "bold", size = 9) +
  annotate(geom = "text", label = "North-East", x = 13, y = 53.4, col = "black", fontface = "bold", size = 9) +
  scale_color_manual(values = c("South-West" = "#009E73", "Center" = "#CC79A7", "North-East" = "#0072B2")) +
  annotation_north_arrow(which_north = "true", location = "br", height = unit(4, "cm"), width = unit(2.2, "cm")) +
  theme_void() +
  theme(legend.position = "none")

ggsave(plot = Germany_map, filename = "~/Documents/ClimateExtremes/CleanCode/PNG_figs/Germany_map.png", device = "png", width = 10, height = 10, dpi = 300, bg = "white")
#ggsave(plot = Germany_map, filename = "~/Documents/ClimateExtremes/CleanCode/PNG_figs/Germany_map.svg", device = "svg", width = 10, height = 10, dpi = 300, bg = "white")

ALB_map <- ggplot() +
  geom_raster(data = DEM_ALB, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_carto_c(palette = "Earth", name = "Altitude m asl") +
  geom_sf(data = Kml_obs.gr$A, col = "black", size = 3) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.1, "cm"), height = unit(0.5, "cm"), text_cex = 1) +
  ylab(NULL) + xlab(NULL) +
  ggtitle("South-West") +
  theme_pubr() +
  theme(legend.position = "right", legend.title = element_text(size = 16), title = element_text(size = 18),
        axis.text.x.bottom = element_text(size = 10), axis.text.y.left = element_text(size = 10), legend.text = element_text(size = 14), axis.ticks = element_blank(),
        panel.border = element_rect(color = "#009E73", fill = NA, size = 1.5), axis.line = element_blank())

ggsave(plot = ALB_map, filename = "~/Documents/ClimateExtremes/CleanCode/PNG_figs/ALB_map.png", device = "png", width = 9, height = 5, dpi = 300, bg = "white")

HAI_map <- ggplot() +
  geom_raster(data = DEM_HAI, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_carto_c(palette = "Earth", name = "Altitude m asl") +
  geom_sf(data = Kml_obs.gr$H, col = "black", size = 3) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.1, "cm"), height = unit(0.5, "cm"), text_cex = 1) +
  ylab(NULL) + xlab(NULL) +
  ggtitle("Center") +
  theme_pubr() +
  theme(legend.position = "right", legend.title = element_text(size = 16), title = element_text(size = 18),
        axis.text.x.bottom = element_text(size = 10), axis.text.y.left = element_text(size = 10), legend.text = element_text(size = 14), axis.ticks = element_blank(),
        panel.border = element_rect(color = "#CC79A7", fill = NA, size = 1.5), axis.line = element_blank())

ggsave(plot = HAI_map, filename = "~/Documents/ClimateExtremes/CleanCode/PNG_figs/HAI_map.png", device = "png", width = 9, height = 9, dpi = 300, bg = "white")

SCH_map <- ggplot() +
  geom_raster(data = DEM_SCH, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_carto_c(palette = "Earth", name = "Altitude m asl") +
  geom_sf(data = Kml_obs.gr$S, col = "black", size = 3) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.1, "cm"), height = unit(0.5, "cm"), text_cex = 1) +
  ylab(NULL) + xlab(NULL) +
  ggtitle("North-East") +
  theme_pubr() +
  theme(legend.position = "right", legend.title = element_text(size = 16), title = element_text(size = 18),
        axis.text.x.bottom = element_text(size = 10), axis.text.y.left = element_text(size = 10), legend.text = element_text(size = 14), axis.ticks = element_blank(),
        panel.border = element_rect(color = "#0072B2", fill = NA, size = 1.5), axis.line = element_blank())

ggsave(plot = SCH_map, filename = "~/Documents/ClimateExtremes/CleanCode/PNG_figs/SCH_map.png", device = "png", width = 9, height = 9, dpi = 300, bg = "white")

##----------------------------------------------------------------Extract value of SPEI 3/12/24

#SPEI data can be downloaded from here: https://digital.csic.es/handle/10261/202305

#import directories of SPEI netCDF
Dir_SPEI <- list.files(path = "~/Documents/ClimateExtremes/SPEI_nc", pattern = ".nc", full.names = T)

SPEI_nc <- lapply(Dir_SPEI, brick)

names(SPEI_nc) <- paste("SPEI", c(3, 12, 24), sep = "_")

#set date interval to get data (everything between April and June)
Low_time <- as.Date(paste(seq(2008, 2018), "04-01", sep = "-"), format = "%Y-%m-%d")
Up_time <- as.Date(paste(seq(2008, 2018), "06-30", sep = "-"), format = "%Y-%m-%d")

#get indices of layers with the data for the desired date interval
SPEI_subs <- lapply(SPEI_nc, function(brk, T1 = Low_time, T2 = Up_time) {
  datetime <- brk@z[[1]]
  subs_lyr <- Map(function(t1, t2) which(datetime < t2 & datetime > t1), t1 = T1, t2 = T2)
  subs_lyr <- unlist(subs_lyr)
})

#check all indices are the same
for(i in seq_len(length(SPEI_subs)-1)) {
  print(all.equal(SPEI_subs[[i]], SPEI_subs[[i+1]]))
}

#as all are the same, only take one vector of indices
SPEI_subs <- SPEI_subs[[1]]

#extract layers
SPEI_nc <- lapply(SPEI_nc, '[[', SPEI_subs)

#get SPEI value from centroids
SPEI_nc.val <- lapply(names(SPEI_nc), function(nm, cntrd = Centroids.gr) {
  df <- extract(SPEI_nc[[nm]], st_coordinates(cntrd))
  df <- t(df)
  Col_nm <- rownames(df)
  Nrow <- nrow(df)
  dim(df) <- NULL
  df_ind <- data.frame(Explo = rep(c("A", "H", "S"), each = Nrow), Datetime = rep(Col_nm, 3), df, Ind_nm = nm)
  colnames(df_ind)[3] <- "Val"
  df_ind$Year <- substr(df_ind$Datetime, start = 2, stop = 5)
  return(df_ind)
}
)

names(SPEI_nc.val) <- names(SPEI_nc)

#get SPEI values only for May (25/10/2022)
SPEI_nc.may <- do.call(rbind, lapply(SPEI_nc.val, function(df) {
  df <- df[grep(pattern = ".05.", x = df$Datetime), ]
  return(df)
}))

row.names(SPEI_nc.may) <- seq_len(nrow(SPEI_nc.may))

#check if data are chronologically ordered
identical(SPEI_nc.may$Year, as.character(rep(seq(2008, 2018), 9)))

#change Ind_nm colname in Type and get rid of _ in Ind_nm
colnames(SPEI_nc.may)[4] <- "Type"

SPEI_nc.may$Type <- sub(pattern = "_", replacement = "", x = SPEI_nc.may$Type)

#create separate data.frame for plotting adding 2019 with no value to make it match with biomass time-series
SPEI_nc.may.plot <- SPEI_nc.may

SPEI_nc.may.plot <- rbind(SPEI_nc.may.plot,
                          data.frame(Explo = rep(c("A", "H", "S"), 3), Datetime = NA, Val = NA,
                                     Type = rep(unique(SPEI_nc.may.plot$Type), each = 3), Year = "2019"))


SPEI_nc.may.plot$Type <- factor(SPEI_nc.may.plot$Type, levels = c("SPEI3", "SPEI12", "SPEI24"))

#plot trend of all SPEI time-scales - here exclude 2019 as this plot goes to appendix
AllSPEI_plot.May <- ggplot(SPEI_nc.may.plot[!SPEI_nc.may.plot$Year %in% c("2008", "2019"), ], aes(x = Year, y = Val, col = Type, group = Type)) +
  geom_line(aes(lty = Type), lwd = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c('SPEI3' = "#1E88E5", 'SPEI12' = "#FFC107", 'SPEI24' = "#004D40")) +
  scale_linetype_manual(values = c('SPEI3' = "longdash", 'SPEI12' = "solid", 'SPEI24' = "dotted")) +
  facet_wrap(~ Explo, labeller = labeller(Explo = c(A = "South-West", H = "Center", S = "North-East"))) +
  ylab("SPEI May") + xlab(NULL) +
  theme_pubr() + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
                       strip.background = element_blank(), legend.position = "top")

ggsave(plot = AllSPEI_plot.May, "~/Documents/ClimateExtremes/CleanCode/SPEI_app.jpeg", device = "jpeg", dpi = 300, width = 20, height = 18, units = "cm")

#palette for colours to associate with drought intensity
Drought_pal <- c(ExtremeDry = "burlywood1", ModerateDry = "lightblue1", Normal = "#38ADC3", ModerateWet = "#206EAF", ExtremeWet = "#00116F")

#spei plot for main manuscript figure
MainSPEI_plot <- ggplot(SPEI_nc.may.plot[SPEI_nc.may.plot$Type %in% c("SPEI3", "SPEI12"), ], aes(x = Year, y = Val, group = Type)) +
  annotate(geom = "rect", ymax = -1.28, ymin = -Inf, xmin = -Inf, xmax = Inf, fill = Drought_pal[["ExtremeDry"]], alpha = 0.7) +
  annotate(geom = "rect", ymax = -0.67, ymin = -1.28, xmin = -Inf, xmax = Inf, fill = Drought_pal[["ModerateDry"]], alpha = 0.7) +
  annotate(geom = "rect", ymax = 0.67, ymin = -0.67, xmin = -Inf, xmax = Inf, fill = Drought_pal[["Normal"]], alpha = 0.7) +
  annotate(geom = "rect", ymax = 1.28, ymin = 0.67, xmin = -Inf, xmax = Inf, fill = Drought_pal[["ModerateWet"]], alpha = 0.7) +
  annotate(geom = "rect", ymax = Inf, ymin = 1.28, xmin = -Inf, xmax = Inf, fill = Drought_pal[["ExtremeWet"]], alpha = 0.7) +
  geom_line(aes(lty = Type), lwd = .8) +
  scale_linetype_manual(values = c("SPEI3" = 3, "SPEI12" = 1), labels = c('SPEI3' = "SPEI-3", 'SPEI12' = "SPEI-12")) +
  geom_point(aes(pch = Type), cex = 3) +
  scale_shape_manual(values = c("SPEI3" = 16, "SPEI12" = 18), labels = c('SPEI3' = "SPEI-3", 'SPEI12' = "SPEI-12")) +
  facet_wrap(~ Explo, labeller = labeller(Explo = c(A = "South-West", H = "Center", S = "North-East"))) +
  ylab("SPEI May") + xlab(NULL) +
  theme_pubr() + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
                       legend.position = "bottom", strip.background = element_blank(), strip.text = element_blank())


##----------------------------------------------------------------Classify value of SPEI 3/12/24

SPEI_cat <- SPEI_nc.may[c("Explo", "Val", "Type", "Year")]

SPEI_cat$Explo <- with(SPEI_cat, ifelse(Explo == "A", "ALB",
                                            ifelse(Explo == "H", "HAI", "SCH")))

colnames(SPEI_cat)[2] <- "SPEI"

SPEI_cat.ls <- lapply(unique(SPEI_cat$Type), function(typ) {
  df <- SPEI_cat[SPEI_cat$Type == typ, ]
  df$SPEI_cat <- with(df, ifelse(SPEI < -1.281552 | SPEI > 1.281552, "Extreme",
                             ifelse((SPEI >= -1.281552 & SPEI <= -0.6744898) | (SPEI <= 1.281552 & SPEI >= 0.6744898), "Moderate", "Normal")))
  df$SPEI_final_cat <- with(df, ifelse(SPEI_cat == "Extreme" & SPEI < 0, paste0(SPEI_cat, "Dry"),
                                       ifelse(SPEI_cat == "Moderate" & SPEI < 0, paste0(SPEI_cat, "Dry"),
                                              ifelse(SPEI_cat == "Extreme" & SPEI > 0, paste0(SPEI_cat, "Wet"),
                                                     ifelse(SPEI_cat == "Moderate" & SPEI > 0, paste0(SPEI_cat, "Wet"), SPEI_cat)))))
  return(df)
  })


names(SPEI_cat.ls) <- unique(SPEI_cat$Type)

#number of droughts (and wet events)
lapply(SPEI_cat.ls, function(.) with(., table(Explo, SPEI_final_cat)))


#select years to test resistance
Resist_yrs <- lapply(SPEI_cat.ls, function(.) {
  df <- .[.$SPEI_final_cat %in% c("ModerateDry", "ExtremeDry"), c("Explo", "Year")]
  return(df)
  })

#select years to test recovery
Recovery_yrs <- lapply(SPEI_cat.ls, Find_recovery_yrs)


##------------------------------------------------------Association between SPEI and median biomass

all(sapply(lapply(SPEI_cat.ls, '[[', "Explo"), function(i) i == rep(c("ALB", "HAI", "SCH"), each = length(seq(2008, 2018)))))

all(sapply(lapply(SPEI_cat.ls, '[[', "Year"), function(i) i == rep(seq(2008, 2018), 3)))

Cor_SPEI_biomass <- do.call(cbind, lapply(names(SPEI_cat.ls), function(i) {
  df <- SPEI_cat.ls[[i]][SPEI_cat.ls[[i]]$Year != "2008", "SPEI", drop = F]
  colnames(df) <- i
  if(i == "SPEI3") {
    df <- data.frame(df, SPEI_cat.ls[[i]][SPEI_cat.ls[[i]]$Year != "2008", c("Explo", "Year")])
    df <- df[c(2, 3, 1)]
    df$Year <- as.integer(df$Year)
    }
  return(df)
  }))

Cor_SPEI_biomass <- dplyr::left_join(Cor_SPEI_biomass, biomass.med.q[c("Explo", "year", "Med")],
                                     by = c("Explo" = "Explo", "Year" = "year"))



ggplot(Cor_SPEI_biomass, aes(x = Med, y = SPEI3)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_wrap(~ Explo)

ggplot(Cor_SPEI_biomass, aes(x = Med, y = SPEI12)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_wrap(~ Explo)

ggplot(Cor_SPEI_biomass, aes(x = Med, y = SPEI24)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_wrap(~ Explo)


