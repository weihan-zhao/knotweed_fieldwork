# The organized codes for uploading to the repository 'knotweed_fieldwork_clustering'

# *****************************************************************************
# Content in this document:

# - loading packages

# - loading data

# - Tidying data for clustering and visualization

# - 0. the situation for the ploidy information for 150 populations

# - 1. Hierarchical clustering for all 150 populations
# for Figure 3 and Supplementary Figure 4, Supplementary Table 7.

### Tidying the data
### Selecting the suitable clustering method
### Hierarchical clustering
### The dendrogram for 150 populations

# - 2. Linear Discriminants Analysis (LDA) for all 150 populations

# - 3. Radial plots for all 150 populations
### Overlap Radial plots together

# - 4. Mapping the clustering results for all 150 populations

# - 5. Hierarchical clustering for the populations include 8x individuals
# for Supplementary Figure 5, 6, and 7

### Tidying the data - Filtering out those not 8x individuals
### Selecting the suitable clustering method
### Hierarchical clustering
### The dendrogram for 126 populations

# - 6. The paired dendrograms, trait9_dend and P8X_dend

# - 7. Linear Discriminants Analysis (LDA) for 126 populations

# - 8. Radial plots for 126 populations
### Overlap Radial plots together for 126 populations

# - 9. Mapping the clustering results for 8x clustering results
# *****************************************************************************

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# loading packages

pacman::p_load(
  dplyr,tibble,readr,stringr,psych,MASS,
  cluster,FactoMineR,factoextra,corrplot,dendextend,
  ggord,ggsci,ggplot2,cowplot,ggrepel,
  car,phia,emmeans,fmsb,
  sf,sp,terra,spData,stars
  ) 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# loading data ####

cat("\014")
rm(list = ls())
options(na.action = "na.fail")
gc()

# Loading .RData
load("loading_data.RData")

# data - this table is just a organized version of Ramona's original table
# nc0 - the shp file from gadm, country level map
# nc1 - the shp file from gadm, only for America and China, with povince(state) level polygens
# gadm: https://gadm.org/

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Tidying data for clustering and visualization

names(data)

# [1] "Range", "Individual_Code", "PopulationID", "Latitude", "Ploidy",

# The traits are:
# [7] "Biomass",           
# [8] "SLA",        
# [9] "Chlorophyll_leaf",
# [10] "Toughness_leaf",  
# [11] "leaf_CN",
# [12] "tannin",
# [13] "alkaloid",        
# [14] "lignin",
# [15] "flavonoids"

# chekcing the ploidy level for three ranges:
table(data$Range,data$Ploidy)
#     4X  6X  8X
# CN  18  10 215
# EU   0  23 223
# US  14  75 138

sum(data %>% 
      filter(Range == "CN") %>% 
      select(Ploidy) %>% 
      is.na()) # 1 NA for CN

sum(data %>% 
      filter(Range == "EU") %>% 
      select(Ploidy) %>% 
      is.na()) # 4 NA for EU

sum(data %>% 
      filter(Range == "US") %>% 
      select(Ploidy) %>% 
      is.na()) # 23 NA for US

count(data,Ploidy)
#   Ploidy   n
# 1     4X  32
# 2     6X 108
# 3     8X 576
# 4   <NA>  28

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 0. the situation for the ploidy information for 150 populations ####

# AtLeastOne8X - whether the individuals within a population have at least One 8X. True for the population that has at least one 8x, FALSE for the population without 8X
# ProportionOf8X - the number of 8X individuals within a population. '5', '4',' 3',' 2'，'1', '0'

test_pop <-  
  data %>% 
  group_by(PopulationID) %>% 
  summarise(AtLeastOne8X = any(Ploidy == "8X", na.rm = TRUE),
            ProportionOf8X = sum(Ploidy == "8X", na.rm = TRUE))

count(test_pop, AtLeastOne8X)
#   AtLeastOne8X     n
# 1 FALSE           24
# 2 TRUE           126

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 1. Hierarchical clustering for all 150 populations ####

# Figure 3 and Supplementary Figure 4, Supplementary Table 7.

# = = = = = = = = = = = = = = = = = = = = 
# Tidying the data

trait9 <- data[,c(1,2,4,7,8,9,10,11,12,13,14,15)] 

trait9_mid <- na.omit(trait9)
trait9_mid$lignin <- as.numeric(trait9_mid$lignin)

trait9_pop <- aggregate(trait9_mid[,-1:-2],by=list(trait9_mid$PopulationID), FUN=mean)
rownames(trait9_pop) <- trait9_pop[,1]
trait9_pop <- trait9_pop[,-1:-2]
rm(trait9_mid)

trait9_pop[1:50,] %>% View()
trait9_pop[51:100,] %>% View()
trait9_pop[101:150,] %>% View()

# reorder the rows, CN - US - EU (along the latitude gradient)
trait9_pop <- rbind(trait9_pop[1:50,],trait9_pop[101:150,],trait9_pop[51:100,])

# Data scaling by column
apply(trait9_pop,2,function (x) {sum(is.na(x))}) 

means <- apply(trait9_pop,2,mean)
sds <- apply(trait9_pop,2,sd)
trait9_pop <- scale(trait9_pop,center=means,scale=sds) %>% as.data.frame()

rm(means)
rm(sds)

apply(trait9_pop,2,function (x) {sum(is.na(x))}) 

### exploration
plot(trait9_pop$Biomass ~ trait9_pop$SLA, data = trait9_pop)
with(trait9_pop,text(trait9_pop$Biomass ~ trait9_pop$SLA, 
                     labels=rownames(trait9_pop),pos=4))

# = = = = = = = = = = = = = = = = = = = = 
# Selecting the suitable clustering method

# Reference:
# https://uc-r.github.io/hc_clustering
# https://www.statology.org/hierarchical-clustering-in-r/#:~:text=Hierarchical%20Clustering%20in%20R%3A%20Step-by-Step%20Example%201%20Load,Linkage%20Method%20to%20Use.%20...%20Weitere%20Artikel...%20

# Since we don’t know beforehand which method will produce the best clusters, 
# we can write a short function to perform hierarchical clustering using several different methods.

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
ac <- function(x) {
  agnes(trait9_pop, method = x)$ac
}
sapply(m, ac)

rm(m)
rm(ac)

# = = = = = = = = = = = = = = = = = = = = 
# Hierarchical clustering

trait9_hier <- dist(trait9_pop, method = "euclidean") %>% 
  hclust(method = "ward.D2")

# have a look
plot(trait9_hier)
plot(trait9_hier,labels=row.names(trait9_pop), main='Hierarchical clustering (ward+Euc)')
plot(trait9_hier, hang=-1, cex = 0.45, 
     labels=row.names(trait9_pop), 
     main='Hierarchical clustering (ward+Euc)')

# calculating the gap statistic: the right number of clusters
trait9_gap <- clusGap(trait9_pop, FUN = hcut, nstart = 25, K.max = 6, B = 1000)
fviz_gap_stat(trait9_gap)

png("trait9_gap.png",width=5000,height=4000,res=800)
fviz_gap_stat(trait9_gap)
dev.off()

# Assigning the right Cluster number to the populations
cut_trait9 <- cutree(trait9_hier, k = 5)
table(cut_trait9)

# Making a table that includes for the Cluster Number for 150 populations
trait9_pop_final <- cbind(trait9_pop, cluster = as.factor(cut_trait9))
count(trait9_pop_final, cluster)

trait9_pop_final <- 
  trait9_pop_final %>% 
  mutate(range = case_when(
    str_detect(row.names(trait9_pop_final),'CN') ~ "CN",
    str_detect(row.names(trait9_pop_final),'US') ~ "US",
    str_detect(row.names(trait9_pop_final),'EU') ~ "EU"
  ))

trait9_pop_final$cluster <- as.factor(trait9_pop_final$cluster)
trait9_pop_final$range <- as.factor(trait9_pop_final$range)

colnames(trait9_pop_final) <- 
  c("Standing biomass","SLA","Leaf chlorophyll",
    "Leaf toughness","Leaf C/N","Leaf tannin","Leaf alkaloid",
    "Leaf lignin","Leaf flavonoids","cluster","range")

# = = = = = = = = = = = = = = = = = = = = 
# The dendrogram for 150 populations

# Adding the ploidy informaiton for all 150 populations (AtLeastOne8X and ProportionOf8X)
trait9_pop_final <- 
  trait9_pop_final %>% 
  tibble::rownames_to_column("PopulationID") %>%
  left_join(test_pop, by='PopulationID')

# select the columns that need to be shown as bars in the dendrogram (trait9_dend_bars, the bars dataframe)
trait9_dend_bars <- 
  trait9_pop_final %>% 
  select(PopulationID, cluster,range, AtLeastOne8X, ProportionOf8X)

# reordered the bars dataframe based on the order of the tips in the dend - labels(trait9_dend)
trait9_dend_bars <-
  trait9_dend_bars %>% 
  .[match(labels(trait9_dend),trait9_dend_bars$PopulationID), ]

# Define the color for mapping three ranges
count(trait9_dend_bars, range)
color_map <- c('CN' = "#000000", 'EU' = "#D3D3D3", 'US' = "#808080")
trait9_dend_bars <- 
  trait9_dend_bars %>% 
  mutate(range_color = color_map[range])
# # check
# trait9_dend_bars$range_color == labels_colors(trait9_dend)

# Define the colorfor mapping AtLeastOne8X
count(trait9_dend_bars, AtLeastOne8X)
color_map <- c('TRUE' = "#000000", 'FALSE' = "white")
trait9_dend_bars <- 
  trait9_dend_bars %>% 
  mutate(AtLeastOne8X_color = color_map[as.character(AtLeastOne8X)])

# Define the color for mapping ProportionOf8X
count(trait9_dend_bars, ProportionOf8X)
color_map <- c('5' = "#000000", 
               '4' = "#333333",
               '3' = "#545454",
               '2' = "#7F7F7F",
               '1' = "#A9A9A9",
               '0' = "white")
trait9_dend_bars <- 
  trait9_dend_bars %>% 
  mutate(ProportionOf8X_color = color_map[as.character(ProportionOf8X)])

# Making the dendrograms with dendextend package

trait9_dend <- trait9_hier %>% as.dendrogram()
trait9_dend <- rotate(trait9_dend, 1:150) %>% 
  color_branches(k=5,
                 col = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#AF58BA")) %>% #"#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#AF58BA"
  set("labels_cex", 0.3)

labels_colors(trait9_dend) <-
  c("#000000", "#D3D3D3", "#808080")[sort_levels_values(
    as.numeric(trait9_pop_final$range)[order.dendrogram(trait9_dend)])]

# png("dend_update.png", width=8000,height=4000,res=800)
pdf("dend_update.pdf", width=10, height=5)
par(mar = c(5,5,5,1))
plot(trait9_dend, 
     main = "Hierarchical clustering for transformed traits matrix", 
     horiz =  FALSE,  nodePar = list(cex = .005))
colored_bars(colors = trait9_dend_bars[,c("ProportionOf8X_color", "AtLeastOne8X_color", "range_color")], 
             dend = trait9_dend, 
             rowLabels = c("Number of 8X individuals","At least one 8X","Range"),
             cex.rowLabels = 0.6,
             y_shift = -2,
             sort_by_labels_order = FALSE)
legend("topright", legend = levels(trait9_pop_final$cluster), 
       fill = c("#1F77B4FF","#AF58BA","#2CA02CFF","#D62728FF","#FF7F0EFF"),
       cex = 0.6) #,"#AF58BA"
legend("topleft", legend = levels(trait9_pop_final$range), 
       fill = c("#000000", "#D3D3D3", "#808080"),
       cex = 0.6)
legend("topleft", legend = levels(factor(trait9_pop_final$AtLeastOne8X, levels = c('TRUE', 'FALSE'))), 
       fill = c("#000000", "white"),
       cex = 0.6)
legend("topleft", legend = levels(factor(trait9_pop_final$ProportionOf8X, levels = c('5', '4', '3', '2', '1', '0'))), 
       fill = c("#000000","#333333","#545454","#7F7F7F","#A9A9A9","white"),
       cex = 0.6)
dev.off()

# need to edit in the Inkscape

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 2. Linear Discriminants Analysis (LDA) for all 150 populations ####

trait9_lda <- lda(cluster ~ ., trait9_pop_final[,-c(1,12,13,14)])
trait9_lda

ggord(ord_in = trait9_lda, grp_in = trait9_pop_final$cluster)

(trait9_figure <- 
    ggord(ord_in = trait9_lda, grp_in = trait9_pop_final$cluster,
          repel = TRUE, force = 10, force_pull = 12, txt = 3,labcol = "black", # label parameter
          vec_ext = 6, veclsz =  0.4, # arrows parameter (extension and size)
          poly = FALSE, polylntyp = "solid", alpha_el = 0.3, # elli parameter
          size = 2, alpha = 0.3, # points parameter
          family = "serif", xlims = c(-10,8), ylims = c(-8,8)) + 
    geom_point(size = 2,alpha = 0.3,
               aes(colour = trait9_pop_final$cluster,
                   fill = trait9_pop_final$cluster,
                   shape = trait9_pop_final$cluster))+
    scale_color_manual(values = c("#1F77B4FF","#AF58BA","#2CA02CFF","#D62728FF","#FF7F0EFF")) +#,"#AF58BA"
    scale_fill_manual(values = c("#1F77B4FF","#AF58BA","#2CA02CFF","#D62728FF","#FF7F0EFF")) +#,"#AF58BA"
    scale_shape_manual(values = c(19,15,17,23,25)) + #,25
    #annotate("text", x = -6, y = -3, label = "trait9",
    #family = "serif", colour = "black", size = 6) +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_hline(yintercept=0, linetype="dotted") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),#Remove Grid
          panel.grid.minor = element_blank(),#Remove Grid
          text = element_text(family="serif",face="bold"), 
          panel.border = element_rect(size = 1.2),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 12, color = "black", face = "bold"),
          axis.title.x = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10))+
    xlab("LD1(51.00%)") +
    ylab("LD2(41.02%)")
)

pdf("trait9_figure.pdf", width=6,height=4)
trait9_figure
dev.off()

trait9_lda
trait9_lda$means
trait9_lda$scaling

# These statistics were included in the Supplementary Table 7.
write.csv(trait9_lda$means, file="trait9_lda_means.csv")
write.csv(trait9_lda$scaling, file="trait9_lda_scaling.csv")

# Camparing the means for traits values of between each pair of clusters
names(trait9_pop_final)
count(trait9_pop_final,cluster)

for (i in 4:12){
  print(i)
  print(paste0("Multiple comparison for ",names(trait9_pop_final)[i]))
  model <- lm(trait9_pop_final[,i] ~ cluster,data = trait9_pop_final)
  a <- Anova(model, type=3)
  b <- testInteractions(model, pairwise="cluster")
  c <- emmeans(model, "cluster",data = trait9_pop_final)
  write.csv(a,file = paste0(names(trait9_pop_final)[i],"_a.csv"))
  write.csv(b,file = paste0(names(trait9_pop_final)[i],"_b.csv"))
  write.csv(c,file = paste0(names(trait9_pop_final)[i],"_c.csv"))
  }

rm(model,a,b,c)

# These statistics were included in the Supplementary Table 7.

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 3. Radial plots for all 150 populations ####

data_radar <- aggregate(trait9_pop_final[,-10:-11],by=list(trait9_pop_final$cluster), FUN=mean)
# check count numbers
count(trait9_pop_final,cluster)

c("#1F77B4FF","#AF58BA","#2CA02CFF","#D62728FF","#FF7F0EFF")
col_names <- names(data_radar)[-1] 
col_names

# Cluster1
data_cluster1 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[1,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[1,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[1,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[1,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[1,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[1,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[1,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[1,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[1,10]))

colnames(data_cluster1) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("cluster1.png",width=4000,height=4000,res=800)
radarchart(data_cluster1, 
           axistype=0,
           seg = 5,
           pcol = 1, # Set the color of the polygon border
           pfcol = "#1F77B4FF", # Set fill color to NA (no fill)
           plwd = 1, # Set the line width of the polygon border
           cglcol="lightgrey", # Set the color of the grid lines
           cglty=1, 
           cglwd=0.8,
           vlcex=0.8, # Set the axis label size
           family = "serif",
           title="Cluster1")
dev.off()

# Cluster2
data_cluster2 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[2,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[2,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[2,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[2,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[2,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[2,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[2,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[2,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[2,10]))

colnames(data_cluster2) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("cluster2.png",width=4000,height=4000,res=800)
radarchart(data_cluster2, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#AF58BA",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster2")
dev.off()

# Cluster3
data_cluster3 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[3,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[3,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[3,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[3,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[3,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[3,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[3,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[3,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[3,10]))

colnames(data_cluster3) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("cluster3.png",width=4000,height=4000,res=800)
radarchart(data_cluster3, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#2CA02CFF",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster3")
dev.off()

# Cluster4
data_cluster4 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[4,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[4,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[4,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[4,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[4,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[4,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[4,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[4,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[4,10]))

colnames(data_cluster4) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("cluster4.png",width=4000,height=4000,res=800)
radarchart(data_cluster4, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#D62728FF",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster4")
dev.off()

# Cluster5
data_cluster5 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[5,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[5,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[5,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[5,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[5,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[5,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[5,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[5,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[5,10]))

colnames(data_cluster5) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("cluster5.png",width=4000,height=4000,res=800)
radarchart(data_cluster5, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#FF7F0EFF",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster5") # line width
dev.off()

# = = = = = = = = = = = = = = = = = = = = 
# Overlap Radial plots together

data_cluster_overlap <- 
  data.frame(
    Standing_biomass = 
      c(max(data_radar[,2]), min(data_radar[,2]), 
        data_radar[1,2], data_radar[2,2], data_radar[3,2], data_radar[4,2], data_radar[5,2]),
    SLA = 
      c(max(data_radar[,3]),min(data_radar[,3]),
        data_radar[1,3],data_radar[2,3],data_radar[3,3],data_radar[4,3],data_radar[5,3]),
    Leaf_chlorophyll = 
      c(max(data_radar[,4]),min(data_radar[,4]),
        data_radar[1,4],data_radar[2,4],data_radar[3,4],data_radar[4,4],data_radar[5,4]),
    Leaf_toughness = 
      c(max(data_radar[,5]),min(data_radar[,5]),
        data_radar[1,5],data_radar[2,5],data_radar[3,5],data_radar[4,5],data_radar[5,5]),
    Leaf_CN = 
      c(max(data_radar[,6]),min(data_radar[,6]),
        data_radar[1,6],data_radar[2,6],data_radar[3,6],data_radar[4,6],data_radar[5,6]),
    Tannins = 
      c(max(data_radar[,7]),min(data_radar[,7]),
        data_radar[1,7],data_radar[2,7],data_radar[3,7],data_radar[4,7],data_radar[5,7]),
    Alkaloids = 
      c(max(data_radar[,8]),min(data_radar[,8]),
        data_radar[1,8],data_radar[2,8],data_radar[3,8],data_radar[4,8],data_radar[5,8]),
    Lignins = 
      c(max(data_radar[,9]),min(data_radar[,9]),
        data_radar[1,9],data_radar[2,9],data_radar[3,9],data_radar[4,9],data_radar[5,9]),
    Flavonoids = 
      c(max(data_radar[,10]),min(data_radar[,10]),
        data_radar[1,10],data_radar[2,10],data_radar[3,10],data_radar[4,10],data_radar[5,10]))

# find the suitable colours
col2rgb("#1F77B4FF")
# red     31
# green  119
# blue   180
col2rgb("#AF58BA")
# red    175
# green   88
# blue   186
col2rgb("#2CA02CFF")
# red     44
# green  160
# blue    44
col2rgb("#D62728FF")
# red    214
# green   39
# blue    40
col2rgb("#FF7F0EFF")
# red    255
# green  127
# blue    14

areas <- 
  c(
    rgb(31/255, 119/255, 180/255, 0.20),
    rgb(175/255, 88/255, 186/255, 0.20),
    rgb(44/255, 160/255, 44/255, 0.20),
    rgb(214/255, 39/255, 40/255, 0.20),
    rgb(255/255, 127/255, 14/255, 0.20)
  )

areas <- c("#c9f1ff","#eddaf0","#d0ead0","#f6cfcf","#ffe3c9")

png("cluster_overlap_plot.png",width=4000,height=4000,res=800)
pdf("cluster_overlap_plot.pdf",width=5,height=5)
radarchart(data_cluster_overlap, 
           axistype=0,
           seg = 5,
           pcol = c("#1F77B4FF", "#AF58BA", "#2CA02CFF", "#D62728FF", "#FF7F0EFF"),
           pfcol = areas,
           pty = 32,
           plty = 1,
           plwd = 2,
           cglcol="grey80",
           cglty=2,
           cglwd=0.5,
           vlcex=0.8,
           family = "serif",
           title=NA)
dev.off()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 4. Mapping the clustering results for all 150 populations ####

pop_geo <- aggregate(data[,c(4,5)],by=list(data$PopulationID), FUN=mean)
colnames(pop_geo) <- c('PopulationID',"Latitude","Lontitude")

# selecting the columns just for mapping
data_map <- 
  trait9_pop_final %>% 
  left_join(pop_geo,by='PopulationID') %>% 
  relocate(Latitude,.after=PopulationID) %>% 
  relocate(Lontitude,.after=Latitude) %>% 
  dplyr::select(PopulationID,range,cluster,Latitude,Lontitude)

nc0
plot(nc0,max.plot = 1)

nc1
plot(nc1,max.plot = 1)

data_mapCN <- filter(data_map,range=="CN")
data_mapEU <- filter(data_map,range=="EU")
data_mapUS <- filter(data_map,range=="US")

c("#1F77B4FF", "#AF58BA", "#2CA02CFF", "#D62728FF", "#FF7F0EFF")

sf_use_s2(FALSE)

# set the rect coordinates:
summary(data_mapCN$Latitude)
summary(data_mapCN$Lontitude)

summary(data_mapEU$Latitude)
summary(data_mapEU$Lontitude)

summary(data_mapUS$Latitude)
summary(data_mapUS$Lontitude)

# china_rect <- 
#   c(xmin = 105, xmax = 125, 
#     ymin = 20, ymax = 40) 
# 
# europe_rect <-
#   c(xmin = 2, xmax = 22,
#     ymin = 42, ymax = 62)
# 
# america_rect <-
#   c(xmin = -86, xmax = -66,
#     ymin = 30, ymax = 50)

# so for a global map, the rect coordinates should be
# xlim = c(-86, 125), ylim = c(15, 66)

rm(data_mapCN,data_mapEU,data_mapUS)

ggplot() +
  geom_sf(data = nc0, color="grey80") +
  geom_sf(data = nc1, color="grey80") +
  coord_sf(xlim = c(-86, 125), ylim = c(15, 66), expand = FALSE) +
  geom_point(data = data_map, size = .8, mapping = aes(x = Lontitude, y = Latitude , fill = cluster, colour = cluster)) +
  scale_color_manual(values = c("#1F77B4FF", "#AF58BA", "#2CA02CFF", "#D62728FF", "#FF7F0EFF")) +
  scale_fill_manual(values = c("#1F77B4FF", "#AF58BA", "#2CA02CFF", "#D62728FF", "#FF7F0EFF")) +
  theme_bw() +
  theme(
    text = element_text(family="sans",face="plain"),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(size = 1),
    axis.text = element_text(size = 12, color = "black", face = "plain"),
    axis.ticks = element_line(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("all150_map.png", units="in", width=10, height=8, dpi=800)
ggsave("all150_map.pdf", units="in", width=10, height=8, dpi=800)
ggsave("all150_map.pdf", units="in", width=10, height=8, dpi=72)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 5. Hierarchical clustering for the populations include 8x individuals ####

# *******************************************
# Adding supplement analysis for only the 8x individuals
# So the first thing to do: keep only the 8x individuals from the data table,
# and fitlering out the populations that do not have any 8x individuals
# *******************************************

# Supplementary Figure 5, 6, and 7

pop_geo <- aggregate(data[,c(4,5)],by=list(data$PopulationID), FUN=mean)
colnames(pop_geo) <- c('PopulationID',"Latitude","Lontitude")

# = = = = = = = = = = = = = = = = = = 
# Tidying the data - Filtering out those not 8x individuals

P8X_test <- 
  data %>% 
  filter(Ploidy == "8X") %>% 
  .[,c(1,2,4,6,7:15)] # 7:15 are the 9 traits
str(P8X_test)
P8X_test$lignin <- as.numeric(P8X_test$lignin)

mid <- na.omit(P8X_test)
P8X_test_pop <- aggregate(mid[,c(-1:-4)],by=list(mid$PopulationID), FUN=mean)
rownames(P8X_test_pop) <- P8X_test_pop[,1]
P8X_test_pop <- P8X_test_pop[,-1]

rm(mid)

# see the number of populations for three ranges
cat(paste0(sum(str_detect(rownames(P8X_test_pop),"CN"))," populations from CN"),"\n",
    paste0(sum(str_detect(rownames(P8X_test_pop),"EU"))," populations from EU"),"\n",
    paste0(sum(str_detect(rownames(P8X_test_pop),"US"))," populations from US"))

# 46 populations from CN 
# 46 populations from EU 
# 34 populations from US

# Scaling by column
# check NAs
apply(P8X_test_pop,2,function (x) {sum(is.na(x))}) 

means <- apply(P8X_test_pop,2,mean)
sds <- apply(P8X_test_pop,2,sd)
P8X_test_pop <- scale(P8X_test_pop,center=means,scale=sds) %>% as.data.frame()
rm(means)
rm(sds)

# check NAs
apply(P8X_test_pop,2,function (x) {sum(is.na(x))}) 

# = = = = = = = = = = = = = = = = = = 
# Selecting the suitable clustering method

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
ac <- function(x) {
  agnes(P8X_test_pop, method = x)$ac
}
sapply(m, ac)

rm(m)
rm(ac)

# = = = = = = = = = = = = = = = = = = 
# Hierarchical clustering

P8X_test_hier <- dist(P8X_test_pop, method = "euclidean") %>% 
  hclust(method = "ward.D2")

plot(P8X_test_hier)
plot(P8X_test_hier,labels=row.names(P8X_test_pop), main='Hierarchical clustering (ward+Euc)')
plot(P8X_test_hier, hang=-1, cex = 0.45, 
     labels=row.names(P8X_test_pop), 
     main='Hierarchical clustering (ward+Euc)')

P8X_test_gap <- clusGap(P8X_test_pop, FUN = hcut, nstart = 25, K.max = 6, B = 1000)
fviz_gap_stat(P8X_test_gap)

png("only8X_gap.png",width=5000,height=4000,res=800)
fviz_gap_stat(P8X_test_gap)
dev.off()

cut_P8X <- cutree(P8X_test_hier, k = 4)
table(cut_P8X)

P8X_test_pop_final <- cbind(P8X_test_pop, cluster = as.factor(cut_P8X))
count(P8X_test_pop_final, cluster)

P8X_test_pop_final <- 
  P8X_test_pop_final %>% 
  mutate(range = case_when(
    str_detect(row.names(P8X_test_pop_final),'CN') ~ "CN",
    str_detect(row.names(P8X_test_pop_final),'US') ~ "US",
    str_detect(row.names(P8X_test_pop_final),'EU') ~ "EU"
  ))

P8X_test_pop_final$cluster <- as.factor(P8X_test_pop_final$cluster)
P8X_test_pop_final$range <- as.factor(P8X_test_pop_final$range)

colnames(P8X_test_pop_final) <- 
  c("Standing biomass","SLA","Leaf chlorophyll",
    "Leaf toughness","Leaf C/N","Leaf tannin","Leaf alkaloid",
    "Leaf lignin","Leaf flavonoids","cluster","range")

png("only8X_hier.png", width=5000,height=4000,res=800)
plot(P8X_test_hier,hang=-2, cex = 0.5, 
     labels=row.names(P8X_test_pop), 
     main='trait Clustering')
rect.hclust(P8X_test_hier, k = 4, border = 2:6)
dev.off()

# = = = = = = = = = = = = = = = = = = 
# The dendrogram for 126 populations

# Adding ProportionOf8X into the final dataframe
# still as a bar under the dendrogram

# adding the Propotion8X into the final dataframe
P8X_test_pop_final <- 
  P8X_test_pop_final %>% 
  tibble::rownames_to_column("PopulationID") %>%
  left_join(test_pop[,c(1,3)], by='PopulationID')

P8X_test_dend_bars <- 
  P8X_test_pop_final %>% 
  dplyr::select(PopulationID, cluster,range, ProportionOf8X)

# reordered the bars dataframe based on the order of the tips in the dend - labels(P8X_test_dend)
P8X_test_dend_bars <-
  P8X_test_dend_bars %>% 
  .[match(labels(P8X_test_dend),P8X_test_dend_bars$PopulationID), ]

# Define the color for mapping range
count(P8X_test_dend_bars, range)
color_map <- c('CN' = "#000000", 'EU' = "#D3D3D3", 'US' = "#808080")
P8X_test_dend_bars <- 
  P8X_test_dend_bars %>% 
  mutate(range_color = color_map[range])

# Define the color for mapping ProportionOf8X
count(P8X_test_dend_bars, ProportionOf8X)
color_map <- c('5' = "#000000", 
               '4' = "#333333",
               '3' = "#545454",
               '2' = "#7F7F7F",
               '1' = "#A9A9A9",
               '0' = "white")
P8X_test_dend_bars <- 
  P8X_test_dend_bars %>% 
  mutate(ProportionOf8X_color = color_map[as.character(ProportionOf8X)])

# plotting
P8X_test_dend <- P8X_test_hier %>% as.dendrogram()

P8X_test_dend <- rotate(P8X_test_dend, 1:126) %>% 
  color_branches(k=4,
                 col = c("#1F77B4FF",'#FF7F0EFF',"#2CA02CFF","#AF58BA")) %>% # ,"#FF7F0EFF"
  set("labels_cex", 0.3)

labels_colors(P8X_test_dend) <-
  c("#000000", "#D3D3D3", "#808080")[sort_levels_values(
    as.numeric(P8X_test_pop_final$range)[order.dendrogram(P8X_test_dend)])]

pdf("P8X_dend.pdf", width=10, height=5)
par(mar = c(5,5,5,1))
plot(P8X_test_dend, 
     main = "Hierarchical clustering for transformed traits matrix", 
     horiz =  FALSE,  nodePar = list(cex = .005))
colored_bars(colors = P8X_test_dend_bars[,c("ProportionOf8X_color", "range_color")], 
             dend = P8X_test_dend, 
             rowLabels = c("Number of 8X individuals","Range"),
             cex.rowLabels = 0.6,
             y_shift = -2,
             sort_by_labels_order = FALSE)
legend("topright", legend = levels(P8X_test_pop_final$cluster), 
       fill = c("#1F77B4FF","#2CA02CFF","#AF58BA",'#FF7F0EFF'),
       cex = 0.6)  # ,"#FF7F0EFF"
legend("topleft", legend = levels(P8X_test_pop_final$range), 
       fill = c("#000000", "#D3D3D3", "#808080"),
       cex = 0.6)
legend("topleft", legend = levels(factor(P8X_test_pop_final$ProportionOf8X, levels = c('5', '4', '3', '2', '1', '0'))), 
       fill = c("#000000","#333333","#545454","#7F7F7F","#A9A9A9","white"),
       cex = 0.6)
dev.off()

# need to edit in the Inkscape

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 6. The paired dendrograms, trait9_dend and P8X_dend ####

# Custom these kendo, and place them in a list
dl <- dendlist(
  trait9_dend %>% 
    set("branches_lty", 2),
  P8X_test_dend %>% 
    set("branches_lty", 2)
)

# png("combined_dend.png", width=9000,height=6000,res=600)
pdf("combined_dend.pdf", width=10, height=6)
tanglegram(dl, 
           intersecting = TRUE,
           # Should the leaves of the two dendrograms be pruned so that the
           # two trees will have the same labels?
           columns_width = c(3,4,3),
           # a vector with three elements, giving the relative sizes of the the three plots (left
           # dendrogram, connecting lines, right dendrogram)
           common_subtrees_color_lines = TRUE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd= FALSE, 
           margin_inner = 2,
           lwd = 0.8,
           lab.cex = 0.6,
           main_left = "under previous clustering for 150 populations",
           main_right = "only clustering 126 populations"
)
legend("topleft", legend = levels(trait9_pop_final$cluster), 
       fill = c("#1F77B4FF","#AF58BA","#2CA02CFF","#D62728FF","#FF7F0EFF"),
       cex = 0.6)
legend("topright", legend = levels(P8X_test_pop_final$cluster), 
       fill = c("#1F77B4FF","#2CA02CFF","#AF58BA",'#FF7F0EFF'),
       cex = 0.6)
dev.off()

rm(dl)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =   
# 7. Linear Discriminants Analysis (LDA) for 126 populations ####

P8X_test_lda <- lda(cluster ~ ., P8X_test_pop_final[,-c(1,12,13)])

ggord(ord_in = P8X_test_lda, grp_in = P8X_test_pop_final$cluster)

(P8X_test_figure <- 
    ggord(ord_in = P8X_test_lda, grp_in = P8X_test_pop_final$cluster,
          repel = TRUE, force = 10, force_pull = 12, txt = 3,labcol = "black", # label parameter
          vec_ext = 6, veclsz =  0.4, # arrows parameter (extension and size)
          poly = FALSE, polylntyp = "solid", alpha_el = 0.3, # elli parameter
          size = 2, alpha = 0.3, # points parameter
          family = "serif", xlims = c(-10,6), ylims = c(-5,8)) + 
    geom_point(size = 2,alpha = 0.3,
               aes(colour = P8X_test_pop_final$cluster,
                   fill = P8X_test_pop_final$cluster,
                   shape = P8X_test_pop_final$cluster))+
    scale_color_manual(values = c("#1F77B4FF","#2CA02CFF","#AF58BA",'#FF7F0EFF')) +#,"#AF58BA"
    scale_fill_manual(values = c("#1F77B4FF","#2CA02CFF","#AF58BA",'#FF7F0EFF')) +#,"#AF58BA"
    scale_shape_manual(values = c(19,15,17,23)) + #,25
    geom_vline(xintercept=0, linetype="dotted") +
    geom_hline(yintercept=0, linetype="dotted") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),#Remove Grid
          panel.grid.minor = element_blank(),#Remove Grid
          text = element_text(family="serif",face="bold"), 
          panel.border = element_rect(size = 1.2),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 12, color = "black", face = "bold"),
          axis.title.x = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10))+
    xlab("LD1(71.24%)") +
    ylab("LD2(17.37%)")
)

pdf("P8X_test_figure.pdf", width=6,height=4)
P8X_test_figure
dev.off()

P8X_test_lda
P8X_test_lda$means
P8X_test_lda$scaling

write.csv(P8X_test_lda$means, file="P8X_test_lda_means.csv")
write.csv(P8X_test_lda$scaling, file="P8X_test_lda_scaling.csv")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 8. Radial plots for 126 populations ####

data_radar <- aggregate(P8X_test_pop_final[,-c(1,11,12,13)],by=list(P8X_test_pop_final$cluster), FUN=mean)
# check count numbers
count(P8X_test_pop_final,cluster)

c("#1F77B4FF","#2CA02CFF","#AF58BA",'#FF7F0EFF')
col_names <- names(data_radar)[-1] 
col_names

# Cluster1
data_cluster1 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[1,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[1,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[1,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[1,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[1,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[1,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[1,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[1,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[1,10]))

colnames(data_cluster1) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("P8X_cluster1.png",width=4000,height=4000,res=800)
radarchart(data_cluster1, 
           axistype=0,
           seg = 5,
           pcol = 1,
           pfcol = "#1F77B4FF",
           plwd = 1,
           cglcol="lightgrey",
           cglty=1, 
           cglwd=0.8,
           vlcex=0.8,
           family = "serif",
           title="Cluster1")
dev.off()

# Cluster2
data_cluster2 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[2,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[2,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[2,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[2,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[2,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[2,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[2,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[2,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[2,10]))

colnames(data_cluster2) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("P8X_cluster2.png",width=4000,height=4000,res=800)
radarchart(data_cluster2, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#2CA02CFF",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster2")
dev.off()

# Cluster3
data_cluster3 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[3,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[3,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[3,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[3,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[3,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[3,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[3,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[3,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[3,10]))

colnames(data_cluster3) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("P8X_cluster3.png",width=4000,height=4000,res=800)
radarchart(data_cluster3, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#AF58BA",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster3")
dev.off()

# Cluster4
data_cluster4 <- data.frame(
  Standing_biomass = c(max(data_radar[,2]), min(data_radar[,2]), data_radar[4,2]),
  SLA = c(max(data_radar[,3]),min(data_radar[,3]),data_radar[4,3]),
  Leaf_chlorophyll = c(max(data_radar[,4]),min(data_radar[,4]),data_radar[4,4]),
  Leaf_toughness = c(max(data_radar[,5]),min(data_radar[,5]),data_radar[4,5]),
  Leaf_CN = c(max(data_radar[,6]),min(data_radar[,6]),data_radar[4,6]),
  Tannins = c(max(data_radar[,7]),min(data_radar[,7]),data_radar[4,7]),
  Alkaloids = c(max(data_radar[,8]),min(data_radar[,8]),data_radar[4,8]),
  Lignins = c(max(data_radar[,9]),min(data_radar[,9]),data_radar[4,9]),
  Flavonoids = c(max(data_radar[,10]),min(data_radar[,10]),data_radar[4,10]))

colnames(data_cluster4) <- 
  c("Biomass","SLA","Chlorophyll",
    "Toughness","C/N","Tannin","Alkaloid",
    "Lignin","Flavonoids")

png("P8X_cluster4.png",width=4000,height=4000,res=800)
radarchart(data_cluster4, axistype=0,
           seg = 5,
           pcol = 1,pfcol = "#FF7F0EFF",plwd = 1,
           cglcol="lightgrey", cglty=1, cglwd=0.8,
           vlcex=0.8, family = "serif",
           title="Cluster4")
dev.off()

rm(data_radar,data_cluster1,data_cluster2,data_cluster3,data_cluster4)

# = = = = = = = = = = = = = = = = = = 
# Overlap Radial plots together for 126 populations

data_cluster_overlap <- 
  data.frame(
    Standing_biomass = 
      c(max(data_radar[,2]), min(data_radar[,2]), 
        data_radar[1,2], data_radar[2,2], data_radar[3,2], data_radar[4,2]),
    SLA = 
      c(max(data_radar[,3]),min(data_radar[,3]),
        data_radar[1,3],data_radar[2,3],data_radar[3,3],data_radar[4,3]),
    Leaf_chlorophyll = 
      c(max(data_radar[,4]),min(data_radar[,4]),
        data_radar[1,4],data_radar[2,4],data_radar[3,4],data_radar[4,4]),
    Leaf_toughness = 
      c(max(data_radar[,5]),min(data_radar[,5]),
        data_radar[1,5],data_radar[2,5],data_radar[3,5],data_radar[4,5]),
    Leaf_CN = 
      c(max(data_radar[,6]),min(data_radar[,6]),
        data_radar[1,6],data_radar[2,6],data_radar[3,6],data_radar[4,6]),
    Tannins = 
      c(max(data_radar[,7]),min(data_radar[,7]),
        data_radar[1,7],data_radar[2,7],data_radar[3,7],data_radar[4,7]),
    Alkaloids = 
      c(max(data_radar[,8]),min(data_radar[,8]),
        data_radar[1,8],data_radar[2,8],data_radar[3,8],data_radar[4,8]),
    Lignins = 
      c(max(data_radar[,9]),min(data_radar[,9]),
        data_radar[1,9],data_radar[2,9],data_radar[3,9],data_radar[4,9]),
    Flavonoids = 
      c(max(data_radar[,10]),min(data_radar[,10]),
        data_radar[1,10],data_radar[2,10],data_radar[3,10],data_radar[4,10]))

col2rgb("#1F77B4FF")
# red     31
# green  119
# blue   180
col2rgb("#2CA02CFF")
# red     44
# green  160
# blue    44
col2rgb("#AF58BA")
# red    175
# green   88
# blue   186
col2rgb("#FF7F0EFF")
# red    255
# green  127
# blue    14

areas <- 
  c(
    rgb(31/255, 119/255, 180/255, 0.20),
    rgb(44/255, 160/255, 44/255, 0.20),
    rgb(175/255, 88/255, 186/255, 0.20),
    rgb(255/255, 127/255, 14/255, 0.20)
  )

png("P8X_cluster_overlap.png",width=4000,height=4000,res=800)
pdf("P8X_cluster_overlapt.pdf",width=5,height=5)
radarchart(data_cluster_overlap, 
           axistype=0,
           seg = 5,
           pcol = c("#1F77B4FF", "#2CA02CFF", "#AF58BA", "#FF7F0EFF"),
           pfcol = areas,  # see below
           pty = 32,
           plty = 1,
           plwd = 2,
           cglcol="grey80",
           cglty= 2, # 2
           cglwd= 0.5, # 0.5
           vlcex=0.8,
           family = "serif",
           title=NA)
dev.off()

rm(data_cluster_overlap,areas,cluster_overlap_plot,data_radar)

# need to edit in the Inkscape

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 9. Mapping the clustering results for 8x clustering results ####

pop_geo <- aggregate(data[,c(4,5)],by=list(data$PopulationID), FUN=mean)
colnames(pop_geo) <- c('PopulationID',"Latitude","Lontitude")

# selecting the columns just for mapping
data_map <- 
  P8X_test_pop_final %>% 
  left_join(pop_geo,by='PopulationID') %>% 
  relocate(Latitude,.after=PopulationID) %>% 
  relocate(Lontitude,.after=Latitude) %>% 
  dplyr::select(PopulationID,range,cluster,Latitude,Lontitude)

data_mapCN <- filter(data_map,range=="CN")
data_mapEU <- filter(data_map,range=="EU")
data_mapUS <- filter(data_map,range=="US")

c("#1F77B4FF", "#2CA02CFF", "#AF58BA", "#FF7F0EFF")

sf_use_s2(FALSE)

# set the rect coordinates:
summary(data_mapCN$Latitude)
summary(data_mapCN$Lontitude)

summary(data_mapEU$Latitude)
summary(data_mapEU$Lontitude)

summary(data_mapUS$Latitude)
summary(data_mapUS$Lontitude)

# china_rect <- 
#   c(xmin = 105, xmax = 125, 
#     ymin = 20, ymax = 40) 
# 
# europe_rect <-
#   c(xmin = 2, xmax = 22,
#     ymin = 42, ymax = 62)
# 
# america_rect <-
#   c(xmin = -86, xmax = -66,
#     ymin = 30, ymax = 50)

# so for a global map, the rect coordinates should be
# xlim = c(-86, 125), ylim = c(15, 66)

rm(data_mapCN,data_mapEU,data_mapUS)

windowsFonts()

ggplot() +
  # basemap first
  geom_sf(data = nc0, color="grey80") +
  geom_sf(data = nc1, color="grey80") +
  coord_sf(xlim = c(-86, 125), ylim = c(15, 66), expand = FALSE) +
  # ggsn::scalebar(data = nc0, location = "bottomleft",
  #                dist = 1000, dist_unit = "km", st.size=2, height=0.01, transform = TRUE, model = 'WGS84')+
  geom_point(data = data_map, size = .8, mapping = aes(x = Lontitude, y = Latitude , fill = cluster, colour = cluster)) +
  scale_color_manual(values = c("#1F77B4FF", "#2CA02CFF", "#AF58BA", "#FF7F0EFF")) +
  scale_fill_manual(values = c("#1F77B4FF", "#2CA02CFF", "#AF58BA", "#FF7F0EFF")) +
  theme_bw() +
  theme(
    text = element_text(family="sans",face="plain"),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(size = 1),
    axis.text = element_text(size = 12, color = "black", face = "plain"),
    axis.ticks = element_line(),
    # axis.title.x = element_text(size = 12),
    axis.title.x = element_blank(),
    # axis.title.y = element_text(size = 12),
    axis.title.y = element_blank(),
    legend.position = "none"
    # legend.position = "bottom",
    # legend.title = element_blank(),
    # legend.text = element_text(size = 12)
  )

ggsave("P8X_map.png", units="in", width=10, height=8, dpi=800)
ggsave("P8X_map.pdf", units="in", width=10, height=8, dpi=800)
ggsave("P8X_map.pdf", units="in", width=10, height=8, dpi=72)
