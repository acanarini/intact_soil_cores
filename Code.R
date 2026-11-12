###Load packages
library(ggplot2)
library(ggpubr)
library(reshape2)
library(forcats)
library(tidyverse)
library(dplyr)
library(ggthemes)
library(RColorBrewer)
## Load data

file <- file.choose(new = FALSE)#Multivariate_MSGROWTH
data.growth <- read.csv(file,header=T,dec=".",sep=",", row.names=1)

########################Data Handling#########################
#first plot data for all ecosystems all biomarkers together
data2<-data.growth[,c(5:28) ] #select only biomarkers data in data2
data2 <- mutate_all(data2, function(x) as.numeric(as.character(x)))
drought_label<-data.growth$Combined
# Transform this data in %

#data_percentage_d <- t(apply(data2, 1, function(x){x*100/sum(x,na.rm=T)}))
data_percentage_d <-data2
data_percentage_d<-as.data.frame(data_percentage_d)
colnames(data_percentage_d) <- c( "i15:0","a15:0" ,"i16:0","16:1w7","16:1w5","16:1wX", "16:0" ,"i17:0", "a17:0" ,"cy17:0","17:0","17:0 MeX", "17:1wX","10Me17:0","unknown17", "18:2w6,9","18:1w9cis", "18:1w9trans", "18:1wXa","18:1wXb","18:0", "10Me18:0","19:1w9","19:1wX")
data_percentage_d$Treatment<-data.growth$Treatment
data_percentage_d$Ecosystem<-data.growth$Ecosystem
data_percentage_d$Combined<-data.growth$Combined
data_percentage_d$Core<-data.growth$Core
##############dca test
library(vegan) 
#ord <- decorana(data_percentage_d[,c(1:24) ]) #input the number of biomarkers used; DCA test: The length of the first DCA axis above 4 you can use CCA, below 3 you can use PCA and RDA, between 3 and 4 you can use both
#ord #change NA values to 0 ?

##########PCA and plot#############
library("FactoMineR")
library("corrplot")
library("factoextra")
res.pca <- PCA(data_percentage_d, quali.sup=25:28, graph = TRUE,scale.unit = TRUE)#run PCA
fviz_eig(res.pca, addlabels = TRUE)#show explained variance by axis
var <- get_pca_var(res.pca)#get variables out of PCA
#### create groups for the PLFA peaks
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.character(res.km$cluster)
colors <- as.character(res.km$cluster)
grp[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)] <- c("gram pos",	"gram pos",	"gram pos",	"gram neg",	"amf",	"unsp",	"general",	"gram pos",	"gram pos",	"gram neg",	"general",	"unsp",	"unsp",	"actino",	"unsp",	"fungi",	"fungi",	"unsp",	"unsp",	"unsp",	"genera",	"actino",	"gram neg",	"unsp"
)
group.colors <- c("gram pos" = "blue", "gram neg" = "yellow", "general" ="grey", "fungi" = "red", "actino" = "white", "amf"="orange","unsp"="green")



var1<-fviz_pca_var(res.pca, 
             col.var = grp, palette = c("#99d8c9","#fdbb84","#e34a33","#a6bddb",
                                        "#c994c7","#31a354","#3182bd","#bdbdbd"),
             repel = TRUE,
             legend.title = "Group")+    # Indiviual fill color 
  theme_classic()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  theme(legend.position = "none")+
  labs(tag = "b)")

var1
#create graph for samples
sam1 <-fviz_pca_ind(res.pca, 
                    # Fill individuals by groups
                    geom.ind = "point",
                    pointshape = 21,
                    pointsize = 2,
                    fill.ind = drought_label,
                    col.ind = "black",
                    addEllipses = TRUE,
                    ellipse.type = "confidence",
                    repel = TRUE,
                    mean.point=FALSE,
                    legend.title = "Treatment",
                    palette = c("#1D2769","#66A0BD" , "#1E5235",  "#A4C269", "#943124", "#f4e093"))+      # Indiviual fill color 
  theme_classic()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+ 
  theme(legend.position = "bottom")+
  labs(x = "PC1 (30.6%)",y="PC2 (24%)")+
  labs(tag = "a)")


sam1
#extract legend
leg2 <- get_legend(sam1)
# Convert to a ggplot and print
leg2<-as_ggplot(leg2)
leg2

sam1<-sam1+theme(legend.position = "none")
sam1 <- sam1 + coord_cartesian()
sam1


#run Permanovas
#Treatment and Ecosystem as variables

model1 <- adonis2(
  data_percentage_d[, 1:24] ~ Treatment * Ecosystem,
  data = data_percentage_d,
  method = "euclidean",
  permutations = 999,
  by = "terms"
)

model1

#run model per site
library(vegan)
library(permute)

ecosystems <- c("Cropland", "Forest", "Meadow")

results_by_ecosystem <- list()

for (eco in ecosystems) {
  
  # Subset one ecosystem
  dat_eco <- subset(data_percentage_d, Ecosystem == eco)
  
  # Restricted permutations within each original core
  perm_eco <- with(
    dat_eco,
    how(nperm = 999, blocks = Core)
  )
  
  # Paired PERMANOVA: Intact vs Sieved
  model_eco <- adonis2(
    dat_eco[, 1:24] ~ Treatment,
    data = dat_eco,
    method = "euclidean",
    permutations = perm_eco
  )
  
  results_by_ecosystem[[eco]] <- model_eco
  
  cat("\n\n============================\n")
  cat("ECOSYSTEM:", eco, "\n")
  cat("============================\n")
  print(model_eco)
}

results_by_ecosystem$Cropland
results_by_ecosystem$Forest
results_by_ecosystem$Meadow




######Calculate Euclidean dissimilarity ################
dis <-vegdist(data2, method="euclid")   
groups <- factor(c(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5), rep(6,5)), labels = c("Cropland","Cropland - Sieved",
                                                                                               "Forest - IC","Forest - Sieved",       
                          
                                                                                                                                                              "Meadow - IC","Meadow - Sieved"))
## Using group centroids
mod3 <- betadisper(dis, groups, type = "centroid")
mod3
permutest(mod3, permutations = 99)
anova(mod3)
plot(mod3)
boxplot(mod3)
plot(TukeyHSD(mod3))

# betadisper data
# functions
# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}

# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}

# betadisper data
get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}

# get betadisper data
betadisper_dat <- get_betadisper_data(mod3)
beta_dist<-as.data.frame(betadisper_dat$distances)
beta_dist$Ecosystem<-data.growth$Ecosystem
beta_dist$Treatment<-data.growth$Treatment

#color=c("dark blue","light blue" ,"dark green",  "light green", "red", "orange")
color<-c("#1D2769","#66A0BD" , "#1E5235",  "#A4C269", "#943124", "#f4e093")

star_p2 <- beta_dist %>%
  filter(Ecosystem %in% c("Forest", "Meadow")) %>%
  group_by(Ecosystem) %>%
  summarise(
    group = first(group),
    distances = max(distances, na.rm = TRUE) * 1.08,
    label = "*"
  )

star_p2$distances[star_p2$Ecosystem == "Forest"] <-
  star_p2$distances[star_p2$Ecosystem == "Forest"] - 0.2

p2<-ggplot(beta_dist, aes(x = group, y=distances, fill=group)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  geom_text(data = star_p2,
            aes(x = group, y = distances, label = label),
            inherit.aes = FALSE,
            position = position_nudge(x = 0.5),
            size = 6) +
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(bquote("Distance to centroids"))+
  labs(tag = "b)")
p2



######Individual group mass specific growth rates ################

file <- file.choose(new = FALSE)#Taxon-specific-bmcgrowth
msgrowth <- read.csv(file,header=T,dec=".",sep=",", row.names=1)

msgrowth$bacteria<-msgrowth$GRAM.POS + msgrowth$GRAM.NEG + msgrowth$ACTINO

p3<-ggplot(msgrowth, aes(x = Combined, y=FUNGI, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Fungi growth rates", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))))+
  labs(tag = "c)")
p3

p4<-ggplot(msgrowth, aes(x = Combined, y=bacteria, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Bacteria growth rates", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))))+
  labs(tag = "d)")
p4


#######plot all together for Fig. 2############ 
library(patchwork)

layout <- c(
  area(1, 1,2,2),
  area(3, 1, 4,2),
  area(1, 3, 2, 4),
  area(3, 3, 4, 4),
  area(5, 1, 5, 4)
)
plot(layout)

all<- free(sam1, side = "lb") +  p3 +p2 +p4 + leg2+
  plot_layout(design = layout)

all

ggsave(plot=all,"Fig. 2.png", width = 9, height = 8)



######Supplementary folder ################
######Multivariate abundance ################
## Load data

file <- file.choose(new = FALSE)#multivariate abundance
data.ab <- read.csv(file,header=T,dec=".",sep=";", row.names=1)

#######Data Handling
#first plot data for all ecosystems all biomarkers together
data3<-data.ab[,c(5:28) ] #select only biomarkers data in data2
data3 <- mutate_all(data3, function(x) as.numeric(as.character(x)))
drought_label<-data.ab$Combined
# Transform this data in %

#data_percentage_d <- t(apply(data2, 1, function(x){x*100/sum(x,na.rm=T)}))
data_percentage <-data3
data_percentage<-as.data.frame(data_percentage)
colnames(data_percentage) <- c( "i15:0","a15:0" ,"i16:0","16:1w7","16:1w5","16:1wX", "16:0" ,"i17:0", "a17:0" ,"cy17:0","17:0","17:0 MeX", "17:1wX","10Me17:0","unknown17", "18:2w6,9","18:1w9cis", "18:1w9trans", "18:1wXa","18:1wXb","18:0", "10Me18:0","19:1w9","19:1wX")
data_percentage$Treatment<-data.ab$Treatment
data_percentage$Ecosystem<-data.ab$Ecosystem
data_percentage$Combined<-data.ab$Combined
data_percentage$Core<-data.ab$Core

##########PCA and plot#############
library("FactoMineR")
library("corrplot")
library("factoextra")
res.pca2 <- PCA(data_percentage, quali.sup=25:28, graph = TRUE,scale.unit = TRUE)#run PCA
fviz_eig(res.pca2, addlabels = TRUE)#show explained variance by axis
var2 <- get_pca_var(res.pca2)#get variables out of PCA
#### create groups for the PLFA peaks
res.km2 <- kmeans(var2$coord, centers = 3, nstart = 25)
grp2 <- as.character(res.km2$cluster)
colors2 <- as.character(res.km2$cluster)
grp2[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)] <- c("gram pos",	"gram pos",	"gram pos",	"gram neg",	"amf",	"unsp",	"general",	"gram pos",	"gram pos",	"gram neg",	"general",	"unsp",	"unsp",	"actino",	"unsp",	"fungi",	"fungi",	"unsp",	"unsp",	"unsp",	"genera",	"actino",	"gram neg",	"unsp"
)



var2<-fviz_pca_var(res.pca2, 
                   col.var = grp2, palette = c("#99d8c9","#fdbb84","#e34a33","#a6bddb",
                                              "#c994c7","#31a354","#3182bd","#bdbdbd"),
                   repel = TRUE,
                   legend.title = "Group")+    # Indiviual fill color 
  theme_classic()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
labs(tag = "c)")

#extract legend
leg4 <- get_legend(var2)
# Convert to a ggplot and print
leg4<-as_ggplot(leg4)
leg4

var2<-var2+theme(legend.position = "none")
var2

#create graph for samples
sam2 <-fviz_pca_ind(res.pca2, 
                    # Fill individuals by groups
                    geom.ind = "point",
                    pointshape = 21,
                    pointsize = 2,
                    fill.ind = drought_label,
                    col.ind = "black",
                    addEllipses = TRUE,
                    ellipse.type = "confidence",
                    repel = TRUE,
                    mean.point=FALSE,
                    legend.title = "Treatment",
                    palette = c("#1D2769","#66A0BD" , "#1E5235",  "#A4C269", "#943124", "#f4e093"))+      # Indiviual fill color 
  theme_classic()+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+ 
  theme(legend.position = "right")+
  labs(x = "PC1 (72.3%)",y="PC2 (17.6%)")+
  labs(tag = "a)")


sam2
#extract legend
leg3 <- get_legend(sam2)
# Convert to a ggplot and print
leg3<-as_ggplot(leg3)
leg3

sam2<-sam2+theme(legend.position = "none")
sam2

#Plot Fig
library(patchwork)

layout <- c(
  area(1, 1,3,3),
  area(4, 1, 5,2),
  area(4, 3, 5,3),
  area(1, 4, 3,5),
  area(4, 4, 5,5)
)
plot(layout)

all<- sam2 + var1 +leg4+ leg3+ var2  +
  plot_layout(design = layout)

all

ggsave(plot=all,"Fig. S2.png", width = 8, height = 8)


#Permanova
#Treatment and Ecosystem as variables


model2 <- adonis2(
  data_percentage[,c(1:24) ] ~ Treatment * Ecosystem,
  data = data_percentage,
  method = "euclidean",
  permutations = 999,
  by = "terms"
)

model2

library(vegan)
library(permute)

ecosystems <- c("Cropland", "Forest", "Meadow")

results_by_ecosystem <- list()

for (eco in ecosystems) {
  
  # Subset one ecosystem
  dat_eco <- subset(data_percentage, Ecosystem == eco)
  
  # Restricted permutations within each original core
  perm_eco <- with(
    dat_eco,
    how(nperm = 999, blocks = Core)
  )
  
  # Paired PERMANOVA: Intact vs Sieved
  model_eco <- adonis2(
    dat_eco[, 1:24] ~ Treatment,
    data = dat_eco,
    method = "euclidean",
    permutations = perm_eco
  )
  
  results_by_ecosystem[[eco]] <- model_eco
  
  cat("\n\n============================\n")
  cat("ECOSYSTEM:", eco, "\n")
  cat("============================\n")
  print(model_eco)
}

results_by_ecosystem$Cropland
results_by_ecosystem$Forest
results_by_ecosystem$Meadow


######Total growth, CUE and respiration - Supp ################
file <- file.choose(new = FALSE)#Resp_Grwoth_cue
resp_cue <- read.csv(file,header=T,dec=".",sep=",", row.names=1)


f1<-ggplot(resp_cue, aes(x = Combined, y=resp, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Mass specific respiration rates", paste(" ("*mg*" " *C*" " * g*" "* micC^-1 * h^-1 * ")"))))+
  labs(tag = "a)")
f1



f2<-ggplot(resp_cue, aes(x = Combined, y=msgrowth, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Mass specific growth rates", paste(" ("*mg*" " *C*" " * g*" "* micC^-1 * h^-1 * ")"))))+
  labs(tag = "b)")
f2

f3<-ggplot(resp_cue, aes(x = Combined, y=cue, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression("CUE"))+
  labs(tag = "c)")
f3

#extract legend
legend1 <- get_legend(f3)
# Convert to a ggplot and print
legend1<-as_ggplot(legend1)
legend1

f3<-f3+theme(legend.position = "none")
f3

#plot
library(patchwork)

layout <- c(
  area(1, 1,2,2),
  area(3, 1, 4,2),
  area(1, 3, 2, 4),
  area(3, 3, 4, 4)
)
plot(layout)

all<- f1 +  f3 +f2 +legend1 +
  plot_layout(design = layout)

all

ggsave(plot=all,"Fig. S3.png", width = 8, height = 6)


######Individual group mass specific - Supp ################
file <- file.choose(new = FALSE)#Taxon_specific_BMgrowth
individual <- read.csv(file,header=T,dec=".",sep=",", row.names=1)


i1<-ggplot(individual, aes(x = Combined, y=GRAM.POS, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Gram positive", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))))+
  labs(tag = "a)")
i1

i2<-ggplot(individual, aes(x = Combined, y=GRAM.NEG, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Gram negative", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))))+
  labs(tag = "c)")
i2

i3<-ggplot(individual, aes(x = Combined, y=ACTINO, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Actinobacteria", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))))+
  labs(tag = "b)")
i3

i4<-ggplot(individual, aes(x = Combined, y=AMF, fill=Combined)) +
  scale_fill_manual(values=color)+
  geom_boxplot(alpha=0.8) + 
  geom_point(shape = 21,
             color = "black", size = 2)+
  facet_grid(. ~ fct_relevel(Ecosystem,'Cropland','Forest','Meadow'), scales = "free", space = "free")+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.margin=unit(c(2,2,15,2), "mm"))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") +
  ylab(expression(atop("Arbuscular mycorrhizae", paste(" ("*mg*" " *C*" " * g*" "* C^-1 * h^-1 * ")"))))+
  labs(tag = "d)")
i4

#extract legend
legend2 <- get_legend(i4)
# Convert to a ggplot and print
legend2<-as_ggplot(legend2)
legend2

i4<-i4+theme(legend.position = "none")
i4

#plot
library(patchwork)

layout <- c(
  area(1, 1,2,2),
  area(3, 1, 4,2),
  area(1, 3, 2, 4),
  area(3, 3, 4, 4),
  area(5, 1, 5, 4)
)
plot(layout)

all<- i1 +  i2 +i3 +i4 + legend2+
  plot_layout(design = layout)

all

ggsave(plot=all,"Fig. S4.png", width = 9, height = 8)



######Statistic ################
library(rstatix)
#datasets:
#Euclidean distance to group centroids =beta_dist
#Fungi and bacteria = msgrowth
#Respiration etc = resp_cue
#individual groups = individual
beta_dist$Core<-individual$Core
#load the data set
d<-individual

d$new<-paste(d$Ecosystem,d$Core)
d$Treatment <- as.factor(d$Treatment)
d$Ecosystem <- as.factor(d$Ecosystem)
d$distance<-beta_dist$distances

#choose the value and make sure you have all the treatment in the data set
d$value<-(d$distance)
#remove outliers if necessary
#d<-d[-c(11),]

#linear mixed effect model
# Install and load necessary packages

library(nlme)
library(multcomp)
library(emmeans)
library(lmtest)
library(qqplotr)
# Assuming your data frame is named 'data' with columns 'response', 'sieving', 'site', and 'subject'

# Fit the linear mixed-effects model with a random intercept for subjects
model <- lme(log(value) ~ Treatment * Ecosystem, random = ~1 | Core, data = d)

# Print the summary of the model
summary(model)
anova(model)

# Create a QQ plot
residuals <- residuals(model)
# Create a QQ plot
qqnorm(residuals)
qqline(residuals, col = 2)


# Diagnostic plot for homoscedasticity
plot(model, pch = 16, col = "blue", which = c(1, 2))

# Residuals vs. Fitted plot for linearity
plot(model, pch = 16, col = "blue", which = 1)

# Check for outliers
hat_values <- influence(model)$hat
if (!all(is.na(hat_values))) {
  outliers <- hat_values > 2 * mean(hat_values, na.rm = TRUE)
} else {
  outliers <- rep(FALSE, length(hat_values))
}
hat_values
# Perform Shapiro-Wilk test for normality of residuals
shapiro.test(residuals(model))

# Bartlett's test for homogeneity of variances
bartlett.test(residuals(model) ~ interaction(d$Treatment, d$Ecosystem))

# Obtain the marginal means for each combination of 'sieving' and 'site'
marginal_means <- emmeans(model, ~ Treatment * Ecosystem)

# Perform pairwise comparisons with Tukey adjustment
posthoc <- pairs(marginal_means, adjust = "tukey")

# Print the summary of the post hoc test
summary(posthoc)

#The emmeans() function is used to obtain the marginal means for each combination of "sieving" and "site," and the pairs() function is used for pairwise comparisons with Tukey adjustment.
