# Byttnerioideae Landmark Analysis 
# Last updated on 29 March 2025 by K Wolcott

# Built with R Studio v 2023.09.1+494 (2023.09.1+494) and R v 4.4.2 (2024-10-31)
# If on Mac, first install XQuartz 2.8.5
# https://www.xquartz.org/

# If need to remove landmark files
# Open terminal, cd to lm directory (one folder per specimen with multiple lm files), and remove lms only used for measurements, not gmm
# LM names to remove: ligule_length, ligule_width, staminode_length, staminode_width
#find . -name '*ligule*' #-delete
#find . -name '*staminode*' #-delete

#1) Load packages using pacman
#install.packages("pacman")
pacman::p_load(devtools, geomorph, rgl, purrr, dplyr, glue, abind, stringr,
               SlicerMorphR, ggplot2, ggforce, concaveman, Rcpp, ape) 

#2) Set working directory  (define your data_wd below)
data_wd = "/Users/katherinewolcott/Documents/r/byttnerioideae/lms"
lms_wd = paste(data_wd, "/whole_flower", sep="")
setwd(lms_wd)

#3) Build df of all landmark files for each specimen
# Get paths for all specimens (1 folder each containing several landmark files)
dirs = list.dirs(".", full.names=FALSE, recursive=FALSE)
    
# Populate df with all landmark files per specimen
all_lms = lapply(dirs, function(y) {
    glue("Fetching landmark files for specimen: {y}")  
    path = paste(y, '/lms/', sep='')
    files = dir(path, pattern = "*.mrk.json")
    glue("Compiling landmark files: {files}") 
    lms = files %>% map_df(~data.frame(read.markups.json(
      file.path(path, .))) %>% mutate_if(is.numeric, as.character))
    rownames(lms) = NULL
    lms["species"] = unlist(strsplit(y,'_'))[1]
    lms["spec_id"] = y
    rbind(lms)
})

# Optional: Check df to make sure no specimens are missing any lms
#for(i in 1:length(all_lms)) {        
  #print(i)
  #print(nrow(all_lms[[i]][1]))             
#}

# Export to file
df = bind_rows(all_lms) # Normalize dataframe with landmarks as rows, coords, spec_id, and lm_id as cols
write.table(df, "all_lms_wout-names.tsv", sep = "\t", row.names = FALSE)

#4) Add specimen and landmark names to all_lms dataframe
master_lm_df = read.table("all_lms_wout-names.tsv", sep = "\t", header = TRUE)

# Make list of landmark names (lm curve name, number of lm points)
# hacky, but avoids file name irregularities from duplicate naming in 3D slicer
lm_names = c(rep("fil_bot", 15), rep("fil_top", 10), rep("OvaryH", 2), 
             rep("OvaryW", 2), rep("petal_curveL", 40), rep("petal_curveR", 40), 
             rep("sepal_length", 20), rep("sepal_width", 10), rep("StyleH", 2))
master_lm_df["lm_name"] = rep(lm_names, length(unique(df$spec_id)))

# Add numeric labels for 141 landmarks for each specimen 
master_lm_df["lm_id"] = rep(c(1:141), length(all_lms))
master_lm_df = master_lm_df[, c(4, 5, 6, 1, 2, 3)] # re-order columns for easier reading

# Optional: Check if there are any NA's
#master_lm_df[!complete.cases(master_lm_df),]

# Write to file File S2 - all_lms
write.table(master_lm_df, "all_lms.tsv", sep = "\t", row.names = FALSE)

#5) Format landmark coordinates as an array for geomorph
df_lms = df[c('X', 'Y', 'Z')] # Take only coodinates from all_lms df
df_lms_numeric = mutate_all(df_lms, function(x) as.numeric(as.character(x))) # Force as numeric
m = data.matrix(df_lms_numeric) # Make into an array
a = arrayspecs(A=m, p=141, k=3, sep=NULL) # Define array shapes

# Inspect data
a[,,1] # Inspect one specimen
a[1,,] # Inspect one landmark for all specimens
a[,1,] # Inspect landmark coord 1 of 3 for all specimens

#6) Define sliding landmarks for curves of flower
setwd(data_wd)
sliders = read.table("lm_curves.tsv", header=TRUE)
sliders = sliders[3:5]

#7) Make specimen info df to sort lms by groups (location, month, self compat, etc)
# Get list of specimen names from directories
fpath_idx = lapply(dirs, function(y) {
  path = paste(y, '/lms/', sep='')
  fpath_idx = list(path)
})
fpath_idx = unlist(fpath_idx)
write.table(fpath_idx, "specimen_idx.tsv", sep = "\t")

# Load in specimen idx file to map specimen names to numbers during analysis
spec_nums = read.table("specimen_idx.tsv", sep = "\t", row.names=1)
spec_nums["spec_num"] = row.names(spec_nums)
spec_nums["spec_id"] = spec_nums$x
spec_nums["spec_id"] = str_replace(spec_nums[,3], '/lms/', '')
spec_nums$spec_id<-tolower(spec_nums$spec_id) 

# Load in specimen id file that maps names to collection info
spec_ids = read.table("map_specimen_ids.tsv", header=T, sep="\t")
spec_ids$spec_id<-tolower(spec_ids$spec_id) 

# Combine files into specimen info dataframe
spec_info <- merge(spec_ids, spec_nums)

#8) Set variables of interest from specimen info dataframe as factors
genus = spec_info$genus = as.factor(spec_info$genus) 
tribe = spec_info$tribe = as.factor(spec_info$tribe) 
growth_form = spec_info$growth_form = as.factor(spec_info$growth_form) 
sepals = spec_info$sepals = as.factor(spec_info$sepals) 
staminodes = spec_info$staminodes = as.factor(spec_info$staminodes) 
staminode_shape = spec_info$staminode_shape = as.factor(spec_info$staminode_shape) 
ligules = spec_info$ligules = as.factor(spec_info$ligules) 
ligule_shape = spec_info$ligule_shape = as.factor(spec_info$ligule_shape) 
filament_condition = spec_info$filament_condition = as.factor(spec_info$filament_condition) 
col_id = spec_info$col_id = as.factor(spec_info$col_id) 
indiv = spec_info$individual =  as.factor(spec_info$individual) 
col_loc = spec_info$col_loc = as.factor(spec_info$col_loc) 
no_stamens = spec_info$no_stamens = as.factor(spec_info$no_stamens) 
no_stamens_after_fusion = spec_info$no_stamens_after_fusion = as.factor(spec_info$no_stamens_after_fusion) 
no_thecae_per_stamen = spec_info$no_thecae_per_stamen = as.factor(spec_info$no_thecae_per_stamen) 
no_thecae_per_whorl = spec_info$no_thecae_per_whorl = as.factor(spec_info$no_thecae_per_whorl) 
no_staminodes = spec_info$no_staminodes = as.factor(spec_info$no_staminodes) 
no_sepals = spec_info$no_sepals = as.factor(spec_info$no_sepals) 
living_collection = spec_info$living_collection = as.factor(spec_info$living_collection) 

#9) Run GPA to align and scale specimens with Procrustes distance with sliding LMs
Y.gpa <- gpagen(a, curves = sliders, ProcD = TRUE)
summary(Y.gpa)
plot(Y.gpa)

#10) Inspect GPA results
# Get mean shape coordinates (reference shape) after GPA alignment
ref = mshape(Y.gpa$coords, na.action = 3)
plot(ref)

# Compare a specimen to mean consensus shape with different methods
# Thin	plate	spline
plotRefToTarget(ref, target, method="TPS")
# Lollipops
plotRefToTarget(ref,	target, method="vector")
# Points
plotRefToTarget(ref,	target, method="points")

# Plot outliers
outliers = plotOutliers(Y.gpa$coords)

# Find the specimen with the closest values to the mean shape
meanspec = findMeanSpec(Y.gpa$coords)
a[,,meanspec] # Inspect mean spec
spec_nums[meanspec,] # Get the mean specimen id

#11) PCA
pca = gm.prcomp(Y.gpa$coords)
sum = summary(pca)
prop_var = sum$PC.summary[2,1:2]

# Plot Vanilla PCA
dtp = data.frame('genus'= genus, # Take first 2 components of PCA
                 'PC1'=pca$x[,1], 'PC2'=pca$x[,2]) 
# Set parameters
p = ggplot(data = dtp, aes(x = PC1, y = PC2, col = genus)) + 
          # Add polygons around each genus
          geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=genus))+
          geom_point() +  
          theme_minimal()
# Add labels
p + labs(x=glue("PC1: {round(prop_var[1], 4)*100}%"), 
         y=glue("PC2: {round(prop_var[2], 4)*100}%"))


# Plot PCA with collection location points as diff shapes and labelled polygons
dtp = data.frame('col_loc'= col_loc, 'genus' = genus, # Take first 2 components of PCA
                 'PC1'=pca$x[,1], 'PC2'=pca$x[,2]) 
# Set parameters
p = ggplot(data = dtp, aes(x = PC1, y = PC2, col = genus)) + 
          # Add polygons around each genus
          geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=genus, label=genus),
                        label.buffer = unit(1, 'mm'), con.cap = 0,
                        show.legend = FALSE) +
          geom_point(aes(shape=as.factor(col_loc))) + 
                    scale_shape_manual(values=seq(0,8))+
          guides(color=FALSE, shape=guide_legend("site")) +
          theme_minimal()
# Add labels
p + labs(x=glue("PC1: {round(prop_var[1], 4)*100}%"), 
         y=glue("PC2: {round(prop_var[2], 4)*100}%"))

#12) Make deformation meshes (manually add to PCA plots outside of R)
# Compare min and max shapes from PCA to mean consensus shape
# Mean consensus shape
ref<-mshape(Y.gpa$coords)
# PC1
plotRefToTarget(ref, pca$shapes$shapes.comp1$min) # min shape to mean
plotRefToTarget(ref, pca$shapes$shapes.comp1$max) # max shape to mean
# PC2
plotRefToTarget(ref, pca$shapes$shapes.comp2$min) # min shape to mean
plotRefToTarget(ref, pca$shapes$shapes.comp2$max) # max shape to mean

#13) Procrustes Anova
# Make a geomorph dataframe
gdf = geomorph.data.frame(Y.gpa, tribe=tribe, genus=genus, no_stamens=no_stamens, 
                          no_stamens_after_fusion=no_stamens_after_fusion,
                          no_thecae_per_stamen=no_thecae_per_stamen, 
                          no_thecae_per_whorl=no_thecae_per_whorl, sepals=sepals, 
                          staminodes=staminodes, staminode_shape=staminode_shape, 
                          ligules=ligules, ligule_shape=ligule_shape,
                          filament_condition=filament_condition,
                          no_sepals=no_sepals, no_staminodes=no_staminodes, 
                          col_id=col_id, col_loc=col_loc, indiv=indiv,
                          growth_form=growth_form, filament_condition=filament_condition)
attributes(gdf)

# Make various linear models, run anova, and plot regressions to assess fit

#A) Null model (no group effects)
fit.null = procD.lm(coords ~ log(Csize), data = gdf)
aov.null = anova(fit.null)
summary(aov.null)
capture.output(summary(aov.null), file="aov_null.txt") # Export results

# Plot allometry - These plots are independent of linear model used
# CAC - Common Allometric Component & Residual Shape Components (RSCs)
plotAllometry(fit.null, size=gdf$Csize, logsz=TRUE, method="CAC", col=genus)
#par(mar=c(5.1, 4.1, 4.1, 8.1))
#legend("topright", inset=c(-.7,0), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

# Size-shape
plotAllometry(fit.null, size=gdf$Csize, logsz=TRUE, method="size.shape", col=genus, bty='L')
legend("topright", inset=c(0,0), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

# Plot allometry - These plots depend on linear model
# Regression scores
plot(fit.null, type="regression", reg.type = "RegScore", predictor=log(gdf$Csize), col=genus)
par(mar=c(5.1, 4.1, 4.1, 8.1))
title("fit.null RegScore - coords ~ log(Csize)")
legend("topright", inset=c(0,0.1), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

# Prediction line
plot(fit.null, type="regression", reg.type = "PredLine", predictor=log(gdf$Csize), col=genus)
par(mar=c(5.1, 4.1, 4.1, 8.1))
title("fit.null PredLine - coords ~ log(Csize)")
legend("topright", inset=c(0,0.3), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

# PLS analysis on shape and size
pls = two.b.pls(log(gdf$Csize), gdf$coords, print.progress=FALSE)
print(pls)
plot(pls, col=genus, bty='L')
par(mar=c(5.1, 4.1, 4.1, 8.1))
legend("topright", inset=c(0,0.6), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#B) Common fit model (see if group effects are common/shared)
fit.common = procD.lm(coords ~ log(Csize) + genus, data = gdf, SS.type = "III", RRPP=TRUE)
aov.common = anova(fit.common)
summary(aov.common)
capture.output(summary(aov.common), file="aov_common.txt") # Export results

# Plot model fit
plotAllometry(fit.common, size = gdf$Csize, logsz=TRUE, method="PredLine", col=genus, bty="L")
legend("topright", inset=c(0,0), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)
plotAllometry(fit.common, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus, bty="L")
legend("topright", inset=c(0,0.5), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#B) Unique fit model (see if group effects are unique)
fit.unique = procD.lm(coords ~ log(Csize) * genus, data = gdf, SS.type = "III", RRPP=TRUE)
aov.unique = anova(fit.unique)
summary(aov.unique)
capture.output(summary(aov.unique), file="aov_unique.txt") # Export results

# Plot model fit
plotAllometry(fit.unique, size = gdf$Csize, logsz=TRUE, method="PredLine", col=genus)
plotAllometry(fit.unique, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus)
legend("topright", inset=c(0,0.5), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#C) Species fit model 
fit.species = procD.lm(coords ~ genus, data = gdf, SS.type = "III", RRPP=TRUE)
aov.species = anova(fit.species)
summary(aov.species)
capture.output(summary(aov.species), file="aov_species.txt") # Export results

# Plot model fit
plotAllometry(fit.species, size = gdf$Csize, logsz=TRUE, method="PredLine", col=genus)
plotAllometry(fit.species, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus)
legend("topright", inset=c(0,0.7), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#D) Tribe fit model 
fit.tribe = procD.lm(coords ~ tribe, data = gdf, SS.type = "III", RRPP=TRUE)
aov.tribe = anova(fit.tribe)
summary(aov.tribe)
capture.output(summary(aov.tribe), file="aov_tribe.txt") # Export results

# Plot model fit
plotAllometry(fit.tribe, size = gdf$Csize, logsz=TRUE, method="PredLine", col=genus)
plotAllometry(fit.tribe, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus)
legend("topright", inset=c(0,0.5), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#E) Complex fit model 1 (group effects of tribe, genus, individual plant, collection location)
fit.complex1 = procD.lm(coords ~ log(Csize) * tribe/genus + indiv + col_loc, data = gdf, SS.type = "III", RRPP=TRUE)
aov.complex1 = anova(fit.complex1)
summary(aov.complex1)
capture.output(summary(aov.complex1), file="aov_complex1.txt") # Export results

plotAllometry(fit.complex1, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus)
legend("topright", inset=c(0,0.5), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#F) Complex fit model  2 (group effects of tribe, genus, anther states, structure shapes)
fit.complex2 = procD.lm(coords ~ log(Csize) * tribe/genus + staminodes + staminode_shape + ligules + ligule_shape + filament_condition, data = gdf, SS.type = "III", RRPP=TRUE)
aov.complex2 = anova(fit.complex2)
summary(aov.complex2)
capture.output(summary(aov.complex2), file="aov_complex2.txt") # Export results

plotAllometry(fit.complex2, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus)
legend("topright", inset=c(0,0.5), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#G) Complex fit model  3 - no size (group effects of tribe, genus, individual plant, col_loc, anther states, structure shapes)
fit.complex3 = procD.lm(coords ~ tribe/genus + staminodes + staminode_shape + ligules + ligule_shape + filament_condition, data = gdf, SS.type = "III", RRPP=TRUE)
aov.complex3 = anova(fit.complex3)
summary(aov.complex3)
capture.output(summary(aov.complex3), file="aov_complex3.txt") # Export results

plotAllometry(fit.complex3, size = gdf$Csize, logsz=TRUE, method="RegScore", col=genus)
legend("topright", inset=c(0,0.5), legend=levels(genus), xpd=TRUE, fill = 1:length(genus), cex=0.5)

#14) Compare model fit using ANOVA
aov.compare.lms = anova(fit.null, fit.common, fit.unique, fit.complex1, fit.complex2, fit.complex3, SS.type = "III")
summary(aov.compare.lms)
capture.output(summary(aov.compare.lms),file="aov_compare_lms.txt")

#15) Pairwise tests for group differences
PW = pairwise(fit = fit.unique, fit.null = fit.common, groups = genus)
summary(PW, confidence = 0.95)
capture.output(summary(PW, confidence = 0.95),file="pw_fit_unique_v_common.txt")

PW = pairwise(fit = fit.unique, fit.null = fit.null, groups = genus)
summary(PW, confidence = 0.95)
capture.output(summary(PW, confidence = 0.95),file="pw_fit_unique_v_null.txt")

#16) Include phylogeny in model fit
# Use this one for plotting a pretty figure
bytt_tree_fig = ape::read.tree(text = "((Guazuma_ulmifolia_Lam.:2.6, 
                           (Theobroma_cacao_L.:0.8,Herrania_umbratica_R.E._Schult.:.76):0.8):1.1, 
                           ((Byttneria_microphylla_Jacq.:1.6,Ayenia_euphrasiifolia_Griseb.:2.1):4.6,Commersonia_bartramia_L._Merr.:5.1):0.2);")
# Use this one for analysis (names must match your data)
bytt_tree = ape::read.tree(text = "((Guazuma:2.6, 
                           (Theobroma:0.8,Herrania:.76):0.8):1.1, 
                           ((Byttneria:1.6,Ayenia:2.1):4.6,Commersonia:5.1):0.2);")
plot.phylo(bytt_tree_fig, cex=0.9, no.margin = TRUE, x.lim=13, y.lim=7.8, label.offset=.1)
add.scale.bar(x=0, y=0.8, cex=0.8)
edgelabels(bytt_tree_fig$edge.length, bg="transparent", col="black", frame="none", adj=c(0.4,-0.5), font=1, cex=0.6)

# Get mean shape for each species
coords.sub = coords.subset(Y.gpa$coords,genus)
means = lapply(coords.sub, mshape)
means = simplify2array(means)

# Run GPA on mean shape per species
mean.gpa = gpagen(means, curves = sliders, ProcD = TRUE)
summary(mean.gpa)
plot(mean.gpa)

# Check model fit of mean species shapes against phylo
gdf_phy = geomorph.data.frame(mean.gpa, phy=bytt_tree)
fit.phy = procD.pgls(coords ~ log(Csize), bytt_tree, SS.type="III", data=gdf_phy)
aov.phy = anova(fit.phy)
summary(aov.phy)
capture.output(summary(aov.phy), file="aov_phy.txt") # Export results

# Plot model fit
predict(fit.phy)
dev.off()
plotAllometry(fit.phy, size = gdf_phy$Csize, logsz=TRUE, method="PredLine", col=1:length(gdf_phy$phy$tip.label))
plotAllometry(fit.phy, size = gdf_phy$Csize, logsz=TRUE, method="RegScore", col=1:length(gdf_phy$phy$tip.label))
legend("topright", inset=c(0,0.5), legend=levels(as.factor(gdf_phy$phy$tip.label)), xpd=TRUE, fill = 1:length(gdf_phy$phy$tip.label), cex=0.5)

attributes(fit.phy) # Note the PGLS object
fit.phy$LM$Cov # the projection matrix derived from the phylogenetic covariance matrix
fit.phy$pgls.fitted # the PGLS fitted values 
fit.phy$GM$pgls.fitted # The same fitted values, in a 3D array

#17) Test for morphological disparity
# Morphological disparity for entire data set
morphol.disparity(coords ~ 1, groups = NULL, data=gdf, iter=999, print.progress = FALSE)
capture.output(morphol.disparity(coords ~ 1, groups = NULL, data=gdf, iter=999, print.progress = FALSE)
               ,file="morphol_disparity_overall.txt")

morphol.disparity(coords ~ tribe, groups=~genus, data=gdf, iter=999, print.progress = FALSE)
capture.output(morphol.disparity(coords ~ tribe, groups=~genus, data=gdf, iter=999, print.progress = FALSE)
                ,file="morphol_disparity_coords_x_tribe.txt")

morphol.disparity(coords ~ log(Csize) * tribe, groups=~genus, data=gdf, iter=999, print.progress = FALSE)
capture.output(morphol.disparity(coords ~ tribe, groups=~genus, data=gdf, iter=999, print.progress = FALSE)
               ,file="morphol_disparity_coords_x_size_x_tribe.txt")

# Incoporating phylogeny
#morphol.disparity(fit.phy, groups=~as.factor(phy$tip.label), data=gdf_phy, iter=999, print.progress = FALSE)
#capture.output(morphol.disparity(coords ~ tribe, groups=~genus, data=gdf, iter=999, print.progress = FALSE)
               #,file="morphol_disparity_coords_x_size_x_tribe.txt")
