library(ggplot2)
library(ggpubr)


#load data: 2Hincorporation.csv
file <- file.choose(new = FALSE)
data <- read.csv(file,header=T,dec=".",sep=",")

#Data description

#select data sets
data_forest<-subset(data, type=="Forest")

##INDIVIDUAL ANALYSIS FOR EACH TREATMENT
##make only one graph for IC and Sieved soil

####SUBSET 1 - INTACT CORE FOREST
data_forest_IC<-subset(data_forest, Treatment=="IC")
#get atom% and time
SF_IC<-data_forest_IC$time
SF1_IC<-SF_IC
vF1_IC<-data_forest_IC$values
#plot
plot (SF1_IC,vF1_IC, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Forest - intact cores", 
      pch=17, col="red")

#model from paper in GCB (Canarini et al., 2020)
MMcurveF1_IC<-formula(vF1_IC~Vmax+((Vin-Vmax)*(exp(-Km*SF1_IC))))
# We can also build a simple data frame from (S,v) data
kinDataF1_IC <- data.frame(SF1_IC,vF1_IC)

bestfitF1_IC <- nls(MMcurveF1_IC, kinDataF1_IC, start=list(Vin=12,Vmax=6,Km=0.05))
summary(bestfitF1_IC)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRangeF_IC <- seq(0,48.83,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineF1_IC <- predict(bestfitF1_IC,list(SF1_IC=SconcRangeF_IC))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsF1_IC <- coef(bestfitF1_IC)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataF1_IC, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted forest - intact core data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRangeF_IC,theorLineF1_IC,col="red")

#####model for A####
#get Vmax and Km
Fin_IC<-bestFitValsF1_IC[1]
Fmax_IC<-bestFitValsF1_IC[2]
Fk_IC<-bestFitValsF1_IC[3]

####generate observed model (eq2) and predicted (eq1)####

eq1_IC <- function(x){Fmax_IC  +(0.0156-Fmax_IC  )*(exp(- Fk_IC*x))}
eq2_IC <- function(x){Fmax_IC  +(Fin_IC-Fmax_IC  )*(exp(- Fk_IC*x))}

###SUBSET 2 - SIEVED FOREST

data_forest_S<-subset(data_forest, Treatment=="S")
#get atom% and time
SF_S<-data_forest_S$time
SF1_S<-SF_S
vF1_S<-data_forest_S$values
#plot
plot (SF1_S,vF1_S, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Forest - sieved soil", 
      pch=17, col="red")

#model from paper in GCB (Canarini et al., 2020)
MMcurveF1_S<-formula(vF1_S~Vmax+((Vin-Vmax)*(exp(-Km*SF1_S))))
# We can also build a simple data frame from (S,v) data
kinDataF1_S <- data.frame(SF1_S,vF1_S)

bestfitF1_S <- nls(MMcurveF1_S, kinDataF1_S, start=list(Vin=12,Vmax=6,Km=0.05))
summary(bestfitF1_S)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRangeF_S <- seq(0,48.99,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineF1_S <- predict(bestfitF1_S,list(SF1_S=SconcRangeF_S))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsF1_S <- coef(bestfitF1_S)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataF1_S, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted forest sieved soil data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRangeF_S,theorLineF1_S,col="red")

#####model for A####
#get Vmax and Km
Fin_S<-bestFitValsF1_S[1]
Fmax_S<-bestFitValsF1_S[2]
Fk_S<-bestFitValsF1_S[3]

####generate observed model (eq2) and predicted (eq1)####


eq1_S <- function(x){Fmax_S  +(0.0156-Fmax_S  )*(exp(- Fk_S*x))}
eq2_S <- function(x){Fmax_S  +(Fin_S-Fmax_S  )*(exp(- Fk_S*x))}

legend_colors <- c("Intact core" = "#1E5235", "Sieved soil" = "#A4C269")

forestgraph <- ggplot() +
  geom_point(data = data_forest_IC, aes(x = time, y = values, fill = "Intact core"),
             shape = 21, size = 3, stroke = 1) +
  geom_point(data = data_forest_S, aes(x = time, y = values, fill = "Sieved soil"),
             shape = 21, size = 3, stroke = 1) +
  stat_function(fun = eq1_IC, geom = "line", colour = "#1E5235", size = 1, linetype = 1) +
  stat_function(fun = eq2_IC, geom = "line", colour = "#1E5235", size = 1, linetype = 2) +
  stat_function(fun = eq1_S,  geom = "line", colour = "#A4C269", size = 1, linetype = 1) +
  stat_function(fun = eq2_S,  geom = "line", colour = "#A4C269", size = 1, linetype = 2) +
  labs(fill = "") +
  scale_fill_manual(values = legend_colors) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title  = element_text(size = 18),
    strip.text = element_text(size = 14)  # if you facet; otherwise harmless
    # plot.title not set, so it stays at its default size
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(title = "Forest") +
  ylab(bquote(' '^2*'H atom%')) + xlab('Time (hours)') +
  ylim(0, 15)
forestgraph
ggsave("H2incorp_forest_green.jpg", plot = last_plot(),
       width = 8, height = 8, units = "in", scale = 0.7)





#select data sets
data_meadow<-subset(data, type=="Meadow")

##INDIVIDUAL ANALYSIS FOR EACH TREATMENT
##make only one graph for IC and Sieved soil

####SUBSET 1 - INTACT CORE MEADOW
data_meadow_IC<-subset(data_meadow, Treatment=="IC")
#get atom% and time
SM_IC<-data_meadow_IC$time
SM1_IC<-SM_IC
vM1_IC<-data_meadow_IC$values
#plot
plot (SM1_IC,vM1_IC, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Meadow - intact cores", 
      pch=17, col="red")

#model from paper in GCB (Canarini et al., 2020)
MMcurveM1_IC<-formula(vM1_IC~Vmax+((Vin-Vmax)*(exp(-Km*SM1_IC))))
# We can also build a simple data frame from (S,v) data
kinDataM1_IC <- data.frame(SM1_IC,vM1_IC)

bestfitM1_IC <- nls(MMcurveM1_IC, kinDataM1_IC, start=list(Vin=10,Vmax=5.5,Km=0.05))
summary(bestfitM1_IC)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRangeM_IC <- seq(0,48.38,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineM1_IC <- predict(bestfitM1_IC,list(SM1_IC=SconcRangeM_IC))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsM1_IC <- coef(bestfitM1_IC)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataM1_IC, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted meadow - intact core data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRangeM_IC,theorLineM1_IC,col="red")

#####model for A####
#get Vmax and Km
Min_IC<-bestFitValsM1_IC[1]
Mmax_IC<-bestFitValsM1_IC[2]
Mk_IC<-bestFitValsM1_IC[3]

####generate observed model (eq4) and predicted (eq3)####

eq3_IC <- function(x){Mmax_IC  +(0.0156-Mmax_IC  )*(exp(- Mk_IC*x))}
eq4_IC <- function(x){Mmax_IC  +(Min_IC-Mmax_IC  )*(exp(- Mk_IC*x))}

###SUBSET 2 - SIEVED MEADOW

data_meadow_S<-subset(data_meadow, Treatment=="S")
#get atom% and time
SM_S<-data_meadow_S$time
SM1_S<-SM_S
vM1_S<-data_meadow_S$values
#plot
plot (SM1_S,vM1_S, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Meadow - sieved soil", 
      pch=17, col="red")

#model from paper in GCB (Canarini et al., 2020)
MMcurveM1_S<-formula(vM1_S~Vmax+((Vin-Vmax)*(exp(-Km*SM1_S))))
# We can also build a simple data frame from (S,v) data
kinDataM1_S <- data.frame(SM1_S,vM1_S)

bestfitM1_S <- nls(MMcurveM1_S, kinDataM1_S, start=list(Vin=10,Vmax=5.5,Km=0.05))
summary(bestfitM1_S)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRangeM_S <- seq(0,48.3,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineM1_S <- predict(bestfitM1_S,list(SM1_S=SconcRangeM_S))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsM1_S <- coef(bestfitM1_S)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataM1_S, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted meadow sieved soil data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRangeM_S,theorLineM1_S,col="red")

#####model for M####
#get Vmax and Km
Min_S<-bestFitValsM1_S[1]
Mmax_S<-bestFitValsM1_S[2]
Mk_S<-bestFitValsM1_S[3]

####generate observed model (eq4) and predicted (eq3)####


eq3_S <- function(x){Mmax_S  +(0.0156-Mmax_S  )*(exp(- Mk_S*x))}
eq4_S <- function(x){Mmax_S  +(Min_S-Mmax_S  )*(exp(- Mk_S*x))}

legend_colors <- c("Intact core" = "#943124", "Sieved soil" = "#f4e093")

meadowgraph <- ggplot() +
  geom_point(data = data_meadow_IC, aes(x = time, y = values, fill = "Intact core"),
             shape = 21, size = 2, stroke = 1) +
  geom_point(data = data_meadow_S, aes(x = time, y = values, fill = "Sieved soil"),
             shape = 21, size = 2, stroke = 1) +
  stat_function(fun = eq3_IC, geom = "line", colour = "#943124", size = 1, linetype = 1) +
  stat_function(fun = eq4_IC, geom = "line", colour = "#943124", size = 1, linetype = 2) +
  stat_function(fun = eq3_S,  geom = "line", colour = "#f4e093", size = 1, linetype = 1) +
  stat_function(fun = eq4_S,  geom = "line", colour = "#f4e093", size = 1, linetype = 2) +
  labs(fill = "") +
  scale_fill_manual(values = c("Intact core" = "#943124", "Sieved soil" = "#f4e093")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title  = element_text(size = 18),
    strip.text = element_text(size = 14)  # if you facet; otherwise harmless
    # plot.title not set, so it stays at its default size
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(title = "Meadow") +
  ylab(bquote(' '^2*'H atom%')) + xlab('Time (hours)') +
  ylim(0, 15)
meadowgraph
ggsave("H2incorp_meadow_red.jpg", plot = last_plot(), 
       width = 8, height = 8, units = "in", scale = 0.7)




#select data sets
data_field<-subset(data, type=="Field")

##INDIVIDUAL ANALYSIS FOR EACH TREATMENT
##make only one graph for IC and Sieved soil

####SUBSET 1 - INTACT CORE FIELD
data_field_IC<-subset(data_field, Treatment=="IC")
#get atom% and time
SA_IC<-data_field_IC$time
SA1_IC<-SA_IC
vA1_IC<-data_field_IC$values
#plot
plot (SA1_IC,vA1_IC, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Agricultural - intact cores", 
      pch=17, col="red")

#model from paper in GCB (Canarini et al., 2020)
MMcurveA1_IC<-formula(vA1_IC~Vmax+((Vin-Vmax)*(exp(-Km*SA1_IC))))
# We can also build a simple data frame from (S,v) data
kinDataA1_IC <- data.frame(SA1_IC,vA1_IC)

bestfitA1_IC <- nls(MMcurveA1_IC, kinDataA1_IC, start=list(Vin=10,Vmax=5.5,Km=0.05))
summary(bestfitA1_IC)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRangeA_IC <- seq(0,50.05,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineA1_IC <- predict(bestfitA1_IC,list(SA1_IC=SconcRangeA_IC))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsA1_IC <- coef(bestfitA1_IC)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataA1_IC, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted agricultural field - intact core data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRangeA_IC,theorLineA1_IC,col="red")

#####model for A####
#get Vmax and Km
Ain_IC<-bestFitValsA1_IC[1]
Amax_IC<-bestFitValsA1_IC[2]
Ak_IC<-bestFitValsA1_IC[3]

####generate observed model (eq2) and predicted (eq1)####

eq5_IC <- function(x){Amax_IC  +(0.0156-Amax_IC  )*(exp(- Ak_IC*x))}
eq6_IC <- function(x){Amax_IC  +(Ain_IC-Amax_IC  )*(exp(- Ak_IC*x))}

###SUBSET 2 - SIEVED FIELD

data_field_S<-subset(data_field, Treatment=="S")
#get atom% and time
SA_S<-data_field_S$time
SA1_S<-SA_S
vA1_S<-data_field_S$values
#plot
plot (SA1_S,vA1_S, 
      xlab="Time (hours)", 
      ylab="Atom % ", 
      main="Agricultural - sieved soil", 
      pch=17, col="red")

#model from paper in GCB (Canarini et al., 2020)
MMcurveA1_S<-formula(vA1_S~Vmax+((Vin-Vmax)*(exp(-Km*SA1_S))))
# We can also build a simple data frame from (S,v) data
kinDataA1_S <- data.frame(SA1_S,vA1_S)

bestfitA1_S <- nls(MMcurveA1_S, kinDataA1_S, start=list(Vin=10,Vmax=5.5,Km=0.05))
summary(bestfitA1_S)
# Build a theoretical line defining the best fit curve
# First, make a finely detailed set of points between 0 and 48, at 0.1 intervals
# These will be the atm% that are used to calculate 
# the predicted atm%
SconcRangeA_S <- seq(0,50.24,0.1)

# Then, calculate the predicted atm% using the predict function
theorLineA1_S <- predict(bestfitA1_S,list(SA1_S=SconcRangeA_S))

# Best fit values of Km and Vmax obtained by coef function, stored in bestFitVals
bestFitValsA1_S <- coef(bestfitA1_S)

# Now plot the data, the best fit line, and put the best fit coefficients in the plot
plot (kinDataA1_S, 
      xlab="Time (hours)", 
      ylab="Atom %", 
      title(main="Fitted agricultural field sieved soil data"))

# Draw the line
# lines() function adds to the existing plot as does points()
# here, we want a line of best fit, not points, so we use lines()
lines(SconcRangeA_S,theorLineA1_S,col="red")

#####model for A####
#get Vmax and Km
Ain_S<-bestFitValsA1_S[1]
Amax_S<-bestFitValsA1_S[2]
Ak_S<-bestFitValsA1_S[3]

####generate observed model (eq2) and predicted (eq1)####


eq5_S <- function(x){Amax_S  +(0.0156-Amax_S  )*(exp(- Ak_S*x))}
eq6_S <- function(x){Amax_S  +(Ain_S-Amax_S  )*(exp(- Ak_S*x))}

legend_colors <- c("Intact core" = "#1D2769", "Sieved soil" = "#66A0BD")

fieldgraph = ggplot()+
  
  geom_point(data = data_field_IC, aes(x = time, y = values, fill="Intact core"), shape = 21, size = 2, stroke = 1)+
  geom_point(data = data_field_S, aes(x = time, y = values, fill = "Sieved soil"),  shape = 21, size = 2, stroke = 1)+
  
  stat_function(fun=eq5_IC, geom="line",colour = "#1D2769", size = 1, linetype = 1)+
  
  stat_function(fun=eq6_IC, geom="line",colour = "#1D2769", size = 1, linetype = 2)+
  stat_function(fun=eq5_S, geom="line",colour = "#66A0BD", size = 1, linetype = 1)+
  
  stat_function(fun=eq6_S, geom="line",colour = "#66A0BD", size = 1, linetype = 2)+
  labs(fill = "") + 
  scale_fill_manual(values = legend_colors) +
  
  theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title  = element_text(size = 18),
    strip.text = element_text(size = 14)  # if you facet; otherwise harmless
    # plot.title not set, so it stays at its default size
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(title = "Cropland") +
  ylab(bquote(' '^2*'H atom%')) + xlab('Time (hours)') +
  ylim(0, 15)
fieldgraph
ggsave("H2incorp_cropland_blue.jpg", plot = last_plot(),
       width = 8, height = 8, units = "in", scale = 0.7)


#Plot all together for Fi S1
rest<-ggarrange(
  forestgraph, 
  meadowgraph,
  fieldgraph, 
  common.legend = TRUE,legend = "bottom",
  ncol=3,
  nrow=1
  
)
rest
#ggsave(plot=rest,"Fig. S2.png", width = 6, height = 4)

ggsave(file="Dexchangefig3.svg", plot=rest, width=10, height=3.8)


#Calculate integrals
###Integral
#get integral of the curve
require(pracma)
valueA_IC<-integrate(eq1_IC, lower = 0, upper = 48)$value
valueA_S<-integrate(eq1_S, lower = 0, upper = 48)$value

valueM_IC<-integrate(eq3_IC, lower = 0, upper = 48)$value
valueM_S<-integrate(eq3_S, lower = 0, upper = 48)$value

valueF_IC<-integrate(eq5_IC, lower = 0, upper = 48)$value
valueF_S<-integrate(eq5_S, lower = 0, upper = 48)$value



valueA_IC
valueA_S
valueM_IC
valueA_S
valueF_IC
valueF_S


###NEED REAL INUCUB TIME, get one value per sample
#get the atm%
totA_IC<-valueA_IC/48
totA_S<-valueA_S/48
totM_IC<-valueM_IC/48
totM_S<-valueM_S/48
totF_IC<-valueF_IC/48
totF_S<-valueF_S/48
totA_IC
totA_S
totM_IC
totM_S
totF_IC
totF_S


#ggsave(file="deuteriumfig3.8_HS.svg", plot=rest, width=10, height=3.8)

