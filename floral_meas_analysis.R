# Byttnerioideae Floral Measurements analysis 
# Last updated on 29 March 2025 by K Wolcott

# Load packages using pacman
#install.packages("pacman")
pacman::p_load(devtools, rgl, purrr, dplyr, glue, abind, stringr,
               SlicerMorphR, ggplot2, ggforce, concaveman, Rcpp) 

install.packages("tidyr")
install.packages("forcats")
install.packages("viridis")
install.packages("hrbrthemes")
install.packages("ggpubr")
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(ggpubr)

# Set working directory  (define your data_wd below)
data_wd = "/Users/katherinewolcott/Documents/r/byttnerioideae/lms"
setwd(data_wd)

# Generate descriptive statistics for Table 1 - flower measurements
fpath = "measurements_raw.tsv"
meas = read.table(fpath, sep="\t", header=TRUE)

# Inspect raw data
head(meas)
str(meas)

# Drop columns with qualifying variables, keep only measurements
meas_cols = na.omit(meas[,c(2,4:19)])

# Get summary stats for all columns
meas_stats = meas_cols %>%
  group_by(species) %>%
  summarise(across(everything(), 
                   list(min=min, mean=mean, median=median, max=max, sd=sd)))

# Transpose and set colnames
t = data.frame(t(meas_stats))
colnames(t) = t[1,]
t = t[-1, ] 

# Export to file - Table 1
write.table(t, "Table_1-measurements.tsv", sep="\t")


# Plot violins of Ovary Width
p = ggplot(meas, aes(x=species, y=OvaryW, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
                   geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Ovary Width") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Ovary Height
p = ggplot(meas, aes(x=species, y=OvaryH, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Ovary Height") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Style Height
p = ggplot(meas, aes(x=species, y=StyleH, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Style Height") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Style Width
p = ggplot(meas, aes(x=species, y=SepalLen, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Sepal Length") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Sepal Width
p = ggplot(meas, aes(x=species, y=SepalW, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Sepal Width (max)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Ligule Length
p = ggplot(meas, aes(x=species, y=LiguleLen, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Ligule Length") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Ligule Width
p = ggplot(meas, aes(x=species, y=LiguleW, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Ligule Width (max)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Mean Petal Curve Length
meas$PetalMean = (meas$PetalCurveL + meas$PetalCurveR) / 2
meas$PetalCurvMean = (meas$PetalCurveLCurv + meas$PetalCurveRCurv) / 2

p = ggplot(meas, aes(x=species, y=PetalMean, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Length (mm)") +
  ggtitle("Petal Curve Length (mean)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p

# Plot violins of Mean Petal Curvature
p = ggplot(meas, aes(x=species, y=PetalCurvMean, col = species)) +
  geom_violin(width=3, size=0.6, trim=FALSE, aes(fill = species), alpha=0.5) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="darkslategrey") +
  xlab("") +
  ylab("Curvature (mm^-1)") +
  ggtitle("Petal Curvature (mean)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p