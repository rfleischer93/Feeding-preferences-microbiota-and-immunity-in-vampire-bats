

# R code for the manuscript 


#   "Human land-use changes induce dietary and gut microbial shifts linked to immunity in vampire bats with implications for host health"

#    by Ramona Fleischer1†*, Christie Jones2,3†, Paula Ledezma-Campos4, Gábor Á. Czirják5, Simone Sommer1, Thomas R. Gillespie2,3,6, Amanda Vicente-Santos6*
#    † Shared first author
#    *  Corresponding authors: ramona.fleischer@uni-ulm.de; amanda.vicente@ou.edu 



# Load packages

x <- c("plyr", "tidyverse", "vegan", "phyloseq", "microbiome", "RColorBrewer", 
    "rstatix", "ggbiplot", "car", "ggplot2", "ape", "vegan", "ggrepel", 
    "DT", "ggpubr", "decontam", "reshape2", "dplyr", "remotes", "btools", "effects",
    "glmmTMB", "ggeffects","sjPlot", "sjlabelled", "sjmisc", "exactRankTests", "nlme", 
    "compositions", "ggrepel", "eulerr", "lme4", "MuMIn")


lapply(x, function(y) {
  # check if installed, if not install
  if (!y %in% installed.packages()[, "Package"])
    install.packages(y)
  
  # load package
  try(require(y, character.only = T), silent = T)
})



# Load data

# load individual dero meta data (n=48)
micro_meta <- read.csv("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\meta_dero.csv", sep=",", header = TRUE)
micro_meta$sex <- recode_factor(micro_meta$sex, F = "F", H = "F", M = "M") # Changing the factor labeling for sex

# load cave-level land-use data (n=11)
cave_land <- read.csv("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\cave_data.csv")



# Distribution of delta carbon values


#################
### Figure S2 ###
#################

# order feeding preference correctly for figure legend

micro_meta$dC13_VPDB <- as.numeric(as.character(micro_meta$dC13_VPDB))

micro_meta$feeding_2c <- factor(micro_meta$feeding_2c, levels = c("wildlife","cattle"))


feeding_preferences_figure <- ggplot((micro_meta %>% drop_na(dC13_VPDB)), aes(y= reorder(sample_name,-dC13_VPDB), x = dC13_VPDB, fill=feeding_2c)) + 
  geom_bar(stat ="identity") +
  xlab(expression(delta^"13"*"C (\u2030)"))+
  scale_fill_manual(name = "feeding preference", labels = c("wildlife blood","cattle blood"), values=c("chartreuse4","#ff5e74"))+  
  ylab("Sample ID")+
  theme_bw(base_size=13)+
  geom_vline(xintercept=-11,linetype=2)



#################
### Figure S3 ###
#################

carbon_nitrogen<-ggplot((micro_meta %>% drop_na(dC13_VPDB)), aes(x = dC13_VPDB, y = dN15_air, group=disturbance_2c)) +
  geom_point(aes(x = dC13_VPDB, y = dN15_air, shape= disturbance_2c, color=feeding_2c), size=2.5) +
  stat_ellipse(aes(x=dC13_VPDB, y=dN15_air, linetype = disturbance_2c),type = "norm") +
  theme_bw() +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  scale_color_manual(name = "feeding preference", labels = c("wildlife blood", "cattle blood"), values=c("chartreuse4","#ff5e74"))+  
  scale_shape_manual(name = "roosting cave 
location", labels = c("disturbed forest", "pristine forest"), values = c(15, 17)) + 
  theme_bw(base_size = 14) +
  scale_linetype_manual(name = "roosting cave 
location", labels = c("disturbed forest", "pristine forest"), values = c(1, 3)) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) 




# Disturbance and feeding preference

# Land-use PCA with cave-level feeding preference

cave_pca <- prcomp(cave_land[,4:7], scale = T)
summary(cave_pca) 


#################
### Figure S6 ###
#################

pca_biplot <- ggbiplot(cave_pca, obs.scale = 1, var.scale = 1, group = cave_land$disturbance_2C, varname.size = 4, ellipse = TRUE,var.axes = T) +
  geom_point(aes(color = cave_land$disturbance_2C, shape = cave_land$disturbance_2C), size = 2.5) +
  scale_color_manual(name = "roosting cave 
location", labels = c("disturbed forest", "pristine forest"), values = c("#fca311","#0077b6"))+  
  scale_shape_manual(name = "roosting cave 
location", labels = c("disturbed forest", "pristine forest"), values = c(15, 17)) + 
  xlab("PC1 (81.4%)") +
  ylab("PC2 (13.2%)") + 
  theme_bw(base_size = 13) +
  geom_text_repel(label=c("DA", "LP", "TM", "TP", "ED", "LA", "BI", "MA", "TI", "GA", "TA")) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())




## Test for an effect between delta carbon and disturbance roost cave location (with PC1 as land-use indicator variable)

# This contains cave spatial data, the PC results for each cave/sample, and sample 13C isotope ratio
cave_micro_si_spatial<-read.csv("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\cavePCA.csv")
 

# GLMM regression analysis of 13C with PC1 
lm_C<- lmer(dC13_VPDB ~ PC1 + (1|cave), data = cave_micro_si_spatial)
summary(lm_C)

# Create a null model with only random effects to run against the mixed model (lm vs lmer function)
null<-lm(dC13_VPDB ~ PC1, data = cave_micro_si_spatial)
summary(null)

r.squaredGLMM(lm_C, null)


# GLMM regression analysis on 15N
lm_N<- lmer(dN15_air ~ PC1  + (1|cave), data = cave_micro_si_spatial)
summary(lm_N) #not significant



# GLMM with 13C, 15N, and PC1 (land-use)
lm_C_N_PC<- lmer(dC13_VPDB ~ PC1 + dN15_air + (1|cave), data = cave_micro_si_spatial)
summary(lm_C_N_PC) #not significant





########################################################################################################################
#################################################### Microbiota ########################################################
########################################################################################################################


## Load phyloseq object

dero_faecal <- readRDS("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\dero_faecal.RDS")
dero_faecal 


## Rarefraction curves

#To check whether sequencing depth has a strong influence on how many taxa we find we can plot rarefaction curves. Rarefaction curves show how many ASVs are detected (y axis) when you randomly subsample more and more reads from the sample (x axis). The appropriate threshold for rarefying is therefore after the curve has flattened, as this indicates that increasing the read depth past this point will not significantly add many more ASVs to the sample.

# generate rarefraction curve function
require(scales)
require(reshape2)
# function to generate rarefaction curves (copied from https://github.com/mahendra-mariadassou/phyloseq-extended/blob/master/R/graphical_methods.R)

ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default 'NULL'. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          'sample_variables(physeq)' ) or taxonomic rank (among the set
  ##          returned by 'rank_names(physeq)').
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          'link{ggplot}', but it can be modified afterward with an
  ##          additional layer using 'scale_color_manual'.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed.
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}



#################
### Figure S5 ###
#################

rarefaction_figure <- ggrare(dero_faecal, step = 1000, se= FALSE)+ xlim(c(0,30000)) +
  geom_vline(xintercept=16000)+
  theme_bw() + 
  theme(legend.position ="none")


#dero_rare <- rarefy_even_depth(dero_faecal, sample.size = 16000, rngseed = 123)
#saveRDS(dero_rare, "dero_rare.RDS")

dero_rare <- readRDS("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\dero_rare.RDS")




################################################## ALPHA DIVERSITY #######################################################

# Here we work with the unrarefied data

dero_faecal 

# now lets add each alpha diversity metric as a column in our sample data within our phyloseq object

sample_data(dero_faecal)$Alpha_observed<-estimate_richness(dero_faecal, measures = "Observed")$Observed
sample_data(dero_faecal)$Alpha_shannon<-estimate_richness(dero_faecal, measures = "Shannon")$Shannon
sample_data(dero_faecal)$Alpha_faiths<-estimate_pd(dero_faecal)$PD

head(sample_data(dero_faecal)) # three new columns added on

sample_data(dero_faecal)$feeding_2c <- factor(sample_data(dero_faecal)$feeding_2c, levels = c("wildlife","cattle"))
sample_data(dero_faecal)$disturbance_2C <- factor(sample_data(dero_faecal)$disturbance_2C, levels = c("p","d"))


#### Feeding preference 

my_theme1 <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

P1 <- ggplot(sample_data(dero_faecal), aes(x = feeding_2c, y = Alpha_observed, fill = feeding_2c)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 14) +
  scale_fill_manual(name = "feeding preference", labels = c("wildlife-blood","cattle-blood"), values=c("chartreuse4","#ff5e74"))+  
  geom_jitter(aes(fill = feeding_2c), width = 0.1, alpha = 0.6) +
  xlab("") +
  ylab("Observed ASVs") +
  my_theme1 + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  geom_signif(comparisons = list(c("cattle", "wildlife")), map_signif_level = TRUE, annotation = c("*"))+
  ylim(0, 295)

P2 <- ggplot(sample_data(dero_faecal), aes(x = feeding_2c, y = Alpha_shannon, fill = feeding_2c)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_bw(base_size = 14) +
  scale_fill_manual(name = "feeding preference", labels = c("wildlife-blood","cattle-blood"), values=c("chartreuse4","#ff5e74"))+  
  geom_jitter(aes(fill = feeding_2c), width = 0.1, alpha = 0.6) +  
  xlab("") +
  ylab("Shannon Index") +
  my_theme1 + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  geom_signif(comparisons = list(c("cattle", "wildlife")), map_signif_level = TRUE, annotation = c("ns"))+
  ylim(0, 2.5)

P3 <- ggplot(sample_data(dero_faecal), aes(x = feeding_2c, y = Alpha_faiths, fill = feeding_2c)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 14) +
  scale_fill_manual(name = "feeding preference", labels = c("wildlife-blood","cattle-blood"), values=c("chartreuse4","#ff5e74"))+  
  geom_jitter(aes(fill = feeding_2c), width = 0.1, alpha = 0.6)  +
  xlab("") +
  ylab("Faith's PD") +
  my_theme1+ 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  geom_signif(comparisons = list(c("cattle", "wildlife")), map_signif_level = TRUE, annotation = c("(*)"))+
  ylim(0, 29)




#################
### Figure 1a ###
#################

figure1a <- ggarrange(P1,P2,P3, ncol =3, common.legend = TRUE, legend = "none", labels = c("a")) 



######### Alpha diversity statistics across all bats ########


metadata_df <- data.frame(sample_data(dero_faecal))
metadata_df$feeding_2c <- factor(metadata_df$feeding_2c, levels = c("cattle","wildlife"))
metadata_df$disturbance_2C <- factor(metadata_df$disturbance_2C, levels = c("d","p"))


# edit columns for figures

colnames(metadata_df)[colnames(metadata_df) == "sex"] = "Sex"
colnames(metadata_df)[colnames(metadata_df) == "feeding_2c"] = "Feeding_pref"
colnames(metadata_df)[colnames(metadata_df) == "disturbance_2C"] = "Cave_location"
colnames(metadata_df)[colnames(metadata_df) == "fecal_storage"] = "Fecal_storage"


##### Alpha observed

m <- lm(log10(Alpha_observed) ~ Feeding_pref + Sex + Seq_depth + Fecal_storage, data = metadata_df)
summary(m)

m.obs <- plot_model(m, show.values = TRUE, value.offset = .3) + ggtitle("Number of ASVs") + theme_bw(base_size = 12) 


m <- lm(log10(Alpha_observed) ~ Cave_location + Sex + Seq_depth + Fecal_storage, data = metadata_df)
summary(m)

m.obs1 <- plot_model(m, show.values = TRUE, value.offset = .3) + ggtitle("Number of ASVs") + theme_bw(base_size = 12) 


##### Alpha Shannon

m <- lm(Alpha_shannon ~ Feeding_pref + Sex + Seq_depth + Fecal_storage, data = metadata_df)
summary(m)

m.sh <- plot_model(m, show.values = TRUE, value.offset = .3) + ggtitle("Shannon Index") + theme_bw(base_size = 12) 


m <- lm(Alpha_shannon ~ Cave_location + Sex + Seq_depth + Fecal_storage, data = metadata_df)
summary(m)

m.sh1 <- plot_model(m, show.values = TRUE, value.offset = .3) + ggtitle("Shannon Index") + theme_bw(base_size = 12) 


##### Alpha Faiths

m <- lm(log10(Alpha_faiths) ~ Feeding_pref + Sex + Seq_depth + Fecal_storage, data = metadata_df)
summary(m)

m.pd <- plot_model(m, show.values = TRUE, value.offset = .3) + ggtitle("Faith's PD") + theme_bw(base_size = 12) 


m <- lm(log10(Alpha_faiths) ~ Cave_location + Sex + Seq_depth + Fecal_storage, data = metadata_df)
summary(m)

m.pd1 <- plot_model(m, show.values = TRUE, value.offset = .3) + ggtitle("Faith's PD") + theme_bw(base_size = 12) 



#################
### Figure S4 ###
#################

alpha_feeding <- ggarrange(m.obs, m.sh, m.pd, nrow =3)

alpha_disturbance <- ggarrange(m.obs1, m.sh1, m.pd1, nrow =3)

alphadiv_figure <- ggarrange(alpha_feeding, alpha_disturbance, ncol=2, labels=c("a","b"))





################################################## BETA DIVERSITY ##################################################


dero_rare <- readRDS("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\dero_rare.RDS")

metadata_df <- data.frame(sample_data(dero_rare))

sample_data(dero_rare)$disturbance_2C <- factor(sample_data(dero_rare)$disturbance_2C, levels = c("p","d"))
sample_data(dero_rare)$feeding_2c <- factor(sample_data(dero_rare)$feeding_2c, levels = c("wildlife","cattle"))

metadata_df$disturbance_2C <- factor(metadata_df$disturbance_2C, levels = c("p","d"))
metadata_df$feeding_2c <- factor(metadata_df$feeding_2c, levels = c("wildlife","cattle"))

################ Jaccard ####################

## phyloseq's ordinate() function will generate a distance matrix and then ordinate it into n - 1 axes

set.seed(1)
dero_mds_jaccard <- phyloseq::ordinate(physeq = dero_rare, method = "MDS", distance = "jaccard")


## how much variation does each axis represent? 

dero_mds_jaccard$values$Relative_eig[1:5] # proportion of variation explained by each axis

## plot ordination (using axes 1 and 2)

## extract axis info

dero_mds_jaccard_vectors <- dero_mds_jaccard$vectors
dero_mds_jaccard_vectors <- as.data.frame(dero_mds_jaccard_vectors)

## only keep axis 1 and 2 for the plot

dero_df_axis <- dero_mds_jaccard_vectors[,1:2]

## add feeding and habitat variables
## check if both df's are ordered the same

metadata_df$feature_id
row.names(dero_df_axis)

metadata_df$feeding_2c
metadata_df$disturbance_2C

## now add

dero_df_axis$feeding_2c <- metadata_df$feeding_2c
dero_df_axis$disturbance_2C <- metadata_df$disturbance_2C

dero_df_axis$disturbance_2C <- factor(dero_df_axis$disturbance_2C, levels = c("p","d"))

## ggplot with ellipses around the feeding preferences

p.jacc <- ggplot(dero_df_axis)+
  geom_point(aes(x=Axis.1, y=Axis.2, color= feeding_2c, shape= disturbance_2C), size = 3) + 
  stat_ellipse(aes(x=Axis.1, y=Axis.2, color= feeding_2c),type = "norm", size = 1) +
  theme_bw(base_size = 15) + 
  scale_color_manual(name = "feeding preference", labels = c("wildlife blood", "cattle blood"), values = c("chartreuse4","#ff5e74"))+
  scale_shape_manual(name = "roosting cave location", labels = c("pristine forest","disturbed forest"), values = c(17, 15)) + 
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Jaccard MDS") + 
  labs(x = "Axis 1 [22.5 %]", y ="Axis 2 [12.0 %]")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=16))




############## Bray Curtis ###########

set.seed(1)
dero_mds_bray <- phyloseq::ordinate(physeq = dero_rare, method = "MDS", distance = "bray")


## how much variation does each axis represent? 

dero_mds_bray$values$Relative_eig[1:5] # proportion of variation explained by each axis

## plot ordination (using axes 1 and 2)

## extract axis info

dero_mds_bray_vectors <- dero_mds_bray$vectors
dero_mds_bray_vectors <- as.data.frame(dero_mds_bray_vectors)

## only keep axis 1 and 2 for the plot

dero_df_axis <- dero_mds_bray_vectors[,1:2]

## add feeding and habitat variables
## check if both df's are ordered the same

metadata_df$feature_id
row.names(dero_df_axis)

metadata_df$feeding_2c
metadata_df$disturbance_2C

## now add

dero_df_axis$feeding_2c <- metadata_df$feeding_2c
dero_df_axis$disturbance_2C <- metadata_df$disturbance_2C


## ggplot with ellipses around the feeding preferences

p.bray <- ggplot(dero_df_axis)+
  geom_point(aes(x=Axis.1, y=Axis.2, color= feeding_2c, shape= disturbance_2C), size = 3) + 
  stat_ellipse(aes(x=Axis.1, y=Axis.2, color= feeding_2c),type = "norm", size = 1) +
  theme_bw(base_size = 15) + 
  scale_color_manual(name = "feeding preference", labels = c("wildlife blood", "cattle blood"), values = c("chartreuse4","#ff5e74"))+
  scale_shape_manual(name = "roosting cave location", labels = c("pristine forest","disturbed forest"), values = c(17, 15)) + 
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Bray MDS") + 
  labs(x = "Axis 1 [28.8 %]", y ="Axis 2 [16.1 %]")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=16))





################## Unifrac ################## 

set.seed(1)
dero_mds_unifrac<- phyloseq::ordinate(physeq = dero_rare, method = "MDS", distance = "unifrac")

## how much variation does each axis represent? 

dero_mds_unifrac$values$Relative_eig[1:5] # proportion of variation explained by each axis


## plot ordination (using axes 1 and 2)
## extract axis info

dero_mds_unifrac_vectors <- dero_mds_unifrac$vectors
dero_mds_unifrac_vectors <- as.data.frame(dero_mds_unifrac_vectors)

## only keep axis 1 and 2 for the plot

dero_df_axis <- dero_mds_unifrac_vectors[,1:2]

## add feeding and habitat variables
## check if both df's are ordered the same

metadata_df$feature_id
row.names(dero_df_axis)

metadata_df$feeding_2c
metadata_df$disturbance_2C

## now add

dero_df_axis$feeding_2c <- metadata_df$feeding_2c
dero_df_axis$disturbance_2C <- metadata_df$disturbance_2C


str(dero_df_axis$feeding_2c)
levels(dero_df_axis$feeding_2c) <- c("wildlife","cattle")

## ggplot with ellipses around the feeding preferences

p.uni <- ggplot(dero_df_axis)+
  geom_point(aes(x=Axis.1, y=Axis.2, color= feeding_2c, shape= disturbance_2C), size = 3) + 
  stat_ellipse(aes(x=Axis.1, y=Axis.2, color= feeding_2c),type = "norm", size = 1) +
  theme_bw(base_size = 15) + 
  scale_color_manual(name = "feeding preference", labels = c("wildlife blood", "cattle blood"), values = c("chartreuse4","#ff5e74"))+
  scale_shape_manual(name = "roosting cave location", labels = c("pristine forest","disturbed forest"), values = c(17, 15)) + 
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("UniFrac MDS") + 
  labs(x = "Axis 1 [20.2 %]", y ="Axis 2 [8.5 %]")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=16)) #+ 
geom_text(x=-0.53, y=0.375, label="R2 = 0.03
  p = 0.1290", size = 4)





################## weighted unifrac ############

set.seed(1)
dero_mds_wunifrac <- phyloseq::ordinate(physeq = dero_rare, method = "MDS", distance = "wunifrac")


## how much variation does each axis represent? 

dero_mds_wunifrac$values$Relative_eig[1:5] # proportion of variation explained by each axis


## plot ordination (using axes 1 and 2)
## extract axis info

dero_mds_wunifrac_vectors <- dero_mds_wunifrac$vectors
dero_mds_wunifrac_vectors <- as.data.frame(dero_mds_wunifrac_vectors)

## only keep axis 1 and 2 for the plot

dero_df_axis <- dero_mds_wunifrac_vectors[,1:2]

## add feeding and habitat variables
## check if both df's are ordered the same

metadata_df$feature_id
row.names(dero_df_axis)

metadata_df$feeding_2c
metadata_df$disturbance_2C

## now add

dero_df_axis$feeding_2c <- metadata_df$feeding_2c
dero_df_axis$disturbance_2C <- metadata_df$disturbance_2C

str(dero_df_axis$feeding_2c)
levels(dero_df_axis$feeding_2c) <- c("wildlife","cattle")

## ggplot with ellipses around the feeding preferences

p.wuni <- ggplot(dero_df_axis)+
  geom_point(aes(x=Axis.1, y=Axis.2, color = feeding_2c, shape= disturbance_2C), size = 3) + 
  stat_ellipse(aes(x=Axis.1, y=Axis.2, color = feeding_2c),type = "norm", size = 1) +
  theme_bw(base_size = 15) + 
  scale_color_manual(name = "feeding preference", labels = c("wildlife blood", "cattle blood"), values = c("chartreuse4","#ff5e74"))+  
  scale_shape_manual(name = "roosting cave location", labels = c("pristine forest","disturbed forest"), values = c(17, 15)) + 
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Weighted UniFrac MDS") + 
  labs(x = "Axis 1 [55.4 %]", y ="Axis 2 [18.4 %]")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=16)) #+ 
geom_text(x=-0.3, y=0.245, label=" R2 = 0.06
  p = 0.0440", size = 4)



## PERMANOVAs

faecal_jaccard_dist <- phyloseq::distance(dero_rare, method = "jaccard")
faecal_bray_dist <- phyloseq::distance(dero_rare, method = "bray")
faecal_unifrac_dist <- phyloseq::distance(dero_rare, method = "unifrac")
faecal_wunifrac_dist <- phyloseq::distance(dero_rare, method = "wunifrac")

metadata_df<-data.frame(sample_data(dero_rare)) # extract metadata as dataframe


## we include feeding preference and habitat in separate models, since the are not independent from each other.

## feeding preference

set.seed(1)
ad_jaccard <- vegan::adonis(faecal_jaccard_dist ~ feeding_2c + sex + fecal_storage, data = metadata_df) 
ad_jaccard
perm_jaccard <- permustats(ad_jaccard)
summary(perm_jaccard)

set.seed(1)
ad_bray <- vegan::adonis(faecal_bray_dist ~ feeding_2c + sex + fecal_storage, data = metadata_df) 
ad_bray
perm_bray<-permustats(ad_bray)
summary(perm_bray)

set.seed(1)
ad_uni <- vegan::adonis(faecal_unifrac_dist ~ feeding_2c + sex + fecal_storage, data = metadata_df) 
ad_uni
perm_uni<-permustats(ad_uni)
summary(perm_uni)

set.seed(1)
ad_wuni <- vegan::adonis(faecal_wunifrac_dist ~ feeding_2c + sex + fecal_storage, data = metadata_df) 
ad_wuni
perm_wuni<-permustats(ad_wuni)
summary(perm_wuni)



## disturbance

set.seed(1)
ad_jaccard <- vegan::adonis(faecal_jaccard_dist ~ disturbance_2C + sex + fecal_storage, data = metadata_df) 
ad_jaccard
perm_jaccard <- permustats(ad_jaccard)
summary(perm_jaccard)

set.seed(1)
ad_bray <- vegan::adonis(faecal_bray_dist ~ disturbance_2C + sex + fecal_storage, data = metadata_df) 
ad_bray
perm_bray<-permustats(ad_bray)
summary(perm_bray)

set.seed(1)
ad_uni <- vegan::adonis(faecal_unifrac_dist ~ disturbance_2C + sex + fecal_storage, data = metadata_df) 
ad_uni
perm_uni<-permustats(ad_uni)
summary(perm_uni)

set.seed(1)
ad_wuni <- vegan::adonis2(faecal_wunifrac_dist ~ disturbance_2C + sex + fecal_storage, data = metadata_df)
ad_wuni
perm_wuni<-permustats(ad_wuni)
summary(perm_wuni)



#################
### Figure 1b ###
#################

figure1b <- ggarrange(p.uni, p.wuni, common.legend = TRUE, legend = "right", labels = c("b"))



#################
### Figure 1 ####
#################

Figure1 <- ggarrange(figure1a, figure1b, common.legend = TRUE, legend = "right", nrow = 2)




############## Heterogeneity among samples ################

## feeding preferences

mod <- betadisper(faecal_jaccard_dist, metadata_df$feeding_2c)
mod
anova(mod)
plot(mod)

mod <- betadisper(faecal_bray_dist, metadata_df$feeding_2c)
mod
anova(mod)
plot(mod)

mod <- betadisper(faecal_unifrac_dist, metadata_df$feeding_2c)
mod
anova(mod)
plot(mod)

mod <- betadisper(faecal_wunifrac_dist, metadata_df$feeding_2c)
mod
anova(mod)
plot(mod)

# no differences


## roosting cave location

mod <- betadisper(faecal_jaccard_dist, metadata_df$disturbance_2C)
mod
anova(mod)
plot(mod)

mod <- betadisper(faecal_bray_dist, metadata_df$disturbance_2C)
mod
anova(mod)
plot(mod)

mod <- betadisper(faecal_unifrac_dist, metadata_df$disturbance_2C)
mod
anova(mod)
plot(mod)

mod <- betadisper(faecal_wunifrac_dist, metadata_df$disturbance_2C)
mod
anova(mod)
plot(mod)

# no differences




################################################## SHARED ASVs #######################################################


## Roosting cave locations


# Code adapted from:
#https://microbiome.github.io/tutorials/core_venn.html

# convert to relative abundances
pseq.rel <- microbiome::transform(dero_rare, "compositional")

# make a list of roosting cave locations
disturbance_roost <- unique(as.character(meta(pseq.rel)$disturbance_2C))


list_all <- c() # an empty object to store information

for (n in disturbance_roost){ # for each variable n in disturbance_roost
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, disturbance_2C == n) # Choose sample from feeding by n
  
  core_m <- core_members(ps.sub, 
                         detection = 0/100, 
                         prevalence = 0/100, include.lowest = FALSE)
  print(paste0("No. of taxa in ", n, " : ", length(core_m))) # print core taxa identified in each feeding group.
  list_all[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}


# reorder
list_all <- list_all[c("p","d")]

# rename columns
names(list_all)[1] <- "roosting cave 
in pristine forest"  
names(list_all)[2] <- "roosting cave 
in disturbed forest"  

# set colours
mycols <- c("#0077b6","#F5BB00")

# plot
ven_b <- plot(venn(list_all), fills = mycols)




## Feeding preferences

# make a list of roosting cave locations
feeding_pref <- unique(as.character(meta(pseq.rel)$feeding_2c))

list_all <- c() # an empty object to store information

for (n in feeding_pref){ # for each variable n in feeding_pref
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, feeding_2c == n) # Choose sample from feeding by n
  
  core_m <- core_members(ps.sub, 
                         detection = 0/100, 
                         prevalence = 0/100, include.lowest = FALSE)
  print(paste0("No. of taxa in ", n, " : ", length(core_m))) # print core taxa identified in each feeding group.
  list_all[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}


# rename columns
names(list_all)[1] <- "wildlife-blood"  
names(list_all)[2] <- "cattle-blood"  

# set colours
mycols <- c(wildlife="#86B181", cattle="#FF8587") 

# plot
ven_a <- plot(venn(list_all), fills = mycols)



###############
### Figure8 ###
###############

venn_figure <- ggarrange(ven_a, ven_b, labels = c("a","b"), nrow=2)



####################################################################################################################################
############################################################## Ancom ###############################################################
####################################################################################################################################

# load Ancom source code from
#https://github.com/FrederickHuangLin/ANCOM
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/

source("ancom_v2.1.R")


# Load unrarefied phyloseq object

dero_faecal 

## data preparation

otu_data <- data.frame(otu_table(dero_faecal))

otu_data <- tibble::rownames_to_column(otu_data, "feature_id")
otu_id = otu_data$feature_id
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

meta_data <- sample_data(dero_faecal)

colnames(meta_data)[1] <- "Sample.ID"


# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = "feeding_2c"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "feeding_2c"; p_adj_method = "BH"; alpha = 0.05
adj_formula =  "disturbance_2C+sex+fecal_storage"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

res$out
list <- res$out
list0.9 <- subset(list, detected_0.9 == "TRUE")
list0.9
list0.8 <- subset(list, detected_0.8 == "TRUE")
list0.8
list0.7 <- subset(list, detected_0.7 == "TRUE")
list0.7
list0.6 <- subset(list, detected_0.6 == "TRUE")
list0.6

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))

# Cutoff values for declaring differentially abundant taxa
cut_off_y = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off_y) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off_y["detected_0.6"], label = "W[0.6]")

fig = res$fig +  
  geom_hline(yintercept =cut_off_y["detected_0.6"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig 


ggsave("ancom_feeding_0.6.pdf")

# save "fig" and "res" results


write.csv(res$fig$data, "ancom_fig.csv")
write.csv(res$out, "ancom_res.csv")



## Volcano plot

ancom_res <- read.csv("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\ancom_res.csv")

### add taxonomix info of ASVs

tax_data <- data.frame(tax_table(dero_faecal))
head(tax_data)
tail(tax_data)

tax_data <- tibble::rownames_to_column(tax_data, "feature_id")
tax_data$feature_id

names(ancom_res)[1] <- "feature_id"  #rename first column in acom table so we can merge both dfs
ancom_res$feature_id

ancom_res_tax <- left_join(ancom_res, tax_data, by="feature_id")

ancom_res_tax <- read.csv("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\ancom_res_tax.csv")
names(ancom_res_tax)


#################
### Figure 2a ###
#################



ancom_fig_a <- ggplot(ancom_res_tax, aes(x = clr, y=W))+
  geom_point(data=ancom_res_tax %>% filter(zero_ind=="No"))+
  geom_point(data=ancom_res_tax %>% filter(detected_0.6=="TRUE") %>% filter(zero_ind=="No"), aes(x = clr, y=W))+
  geom_hline(yintercept = cut_off_y["detected_0.6"], linetype = "dashed") +
  geom_text_repel(data=ancom_res_tax %>% filter(detected_0.6=="TRUE") %>% filter(zero_ind=="No"),
                  # Filter data first
                  aes(label=Genus), size=rel(4.7), fontface = "italic")+
  theme_bw()+
  labs(x="CLR mean difference", y="W statistic")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill="none")+
  theme(text = element_text(size = 16), legend.position = "none")



## Check in which group the ASV is more / less abundant in absolute counts

# first extract otu-data and metadata
otu <- otu_table(dero_faecal)
meta <- sample_data(dero_faecal)

# add +1 to all values in the dataframe, then we are able to compare abundances later properly 
otu_1 <- otu + 1

# transpose, because we need samples as rows, as in the metadata file
otu_trans <- t(as.data.frame(otu_1))
otu_trans <- as.data.frame(otu_trans)

# check both files are in the same order
rownames(otu_trans)
rownames(meta)

# now add feeding column to otu-table
otu_trans$feeding_2c <- as.data.frame(meta)$feeding_2c

# check each ASV by feeding-status
otu_trans_cattle <- subset(otu_trans, feeding_2c == "cattle")
otu_trans_wildilfe <- subset(otu_trans, feeding_2c == "wildlife")


# OTU 131daaddbb4eba63347cf2c984e775f1
# Helicobacter

sum(otu_trans_cattle$"131daaddbb4eba63347cf2c984e775f1")   # 37 132
sum(otu_trans_wildilfe$"131daaddbb4eba63347cf2c984e775f1") # 866

# OTU 29f5d7605cf4510c64c7b4fb1585949d
# Edwardsiella

sum(otu_trans_cattle$"29f5d7605cf4510c64c7b4fb1585949d")   # 39 798
sum(otu_trans_wildilfe$"29f5d7605cf4510c64c7b4fb1585949d") # 196 408



## Boxplot of ASVs identified in Ancom at 60% 

## Add ASVs found at 60% in ancom

dero_tax <- tax_table(dero_faecal)

names <- rownames(dero_tax)
test<- as.matrix(cbind(dero_tax, names))
tax_table(dero_faecal) <- test # map back

my_ancom_asvs <- subset_taxa(dero_faecal, names=="131daaddbb4eba63347cf2c984e775f1"|
                               names=="29f5d7605cf4510c64c7b4fb1585949d")



my_ancom_asvs <- subset_taxa(dero_faecal, names=="131daaddbb4eba63347cf2c984e775f1"|
                               names=="29f5d7605cf4510c64c7b4fb1585949d")


my_ancom_asvs.df <- data.frame(otu_table(my_ancom_asvs))

# add +1 to all values in the dataframe, then we are able to compare abundances later properly with log transformation
my_ancom_asvs.df_1 <- my_ancom_asvs.df + 1  # first save into new object, double check is correct
my_ancom_asvs.df <- my_ancom_asvs.df_1      # then rename

my_ancom_asvs.df  = cbind(as(my_ancom_asvs.df, "data.frame"), as(tax_table(my_ancom_asvs)[rownames(my_ancom_asvs.df), ], "matrix"))

my_ancom_asvs.df = subset(my_ancom_asvs.df, select = -c(Kingdom,Phylum, Class, Order, Genus, Species, names))
my_ancom_asvs.df = subset(my_ancom_asvs.df, select = -c(Family) )

my_ancom_asvs.df<-data.frame(t(my_ancom_asvs.df))
head(my_ancom_asvs.df)
tail(my_ancom_asvs.df)

# check if rows are ordered similar in metadata
metadata<-data.frame(sample_data(my_ancom_asvs))
head(metadata)
tail(metadata)

# both metadata and asv-df are in the same order, thus we can bind together
metadata1<-cbind(metadata, my_ancom_asvs.df)

# Change colnames of ASV columns to correct and short names
colnames(metadata1)[colnames(metadata1) %in% c("X131daaddbb4eba63347cf2c984e775f1", "X29f5d7605cf4510c64c7b4fb1585949d")] <- c("Helicobacter", "Edwardsiella")

#convert to long format
data_long <- gather(metadata1, ASV_ID, Abundance,  Helicobacter:Edwardsiella, factor_key=TRUE)

table(data_long$ASV_ID) # check it worked

data_long$Abundance <- as.numeric(data_long$Abundance)


data_long <- read.csv("D:\\Ramona_F\\3PhD Project\\VampireBats\\publication_files\\ancom_data_long.csv")


data_long$feeding_2c
data_long$feeding_2c <- factor(data_long$feeding_2c, levels = c("wildlife","cattle"))


#################
### Figure 2b ###
#################

ancom_fig_b <- ggplot(data_long , aes(x = feeding_2c, y = Abundance))+
  geom_boxplot(outlier.shape = NA, aes(fill=feeding_2c))+
  scale_fill_manual(name = "feeding preference", labels = c("wildlife-blood","cattle-blood"), values=c("chartreuse4","#ff5e74"))+  
  geom_jitter(width=0.05, size = 1.5, alpha = 0.6)+
  facet_wrap(~ASV_ID, scales = "free")+
  scale_y_log10()+
  theme_bw(base_size=16)+
  labs(x="",y="log(Total Abundance)") +
  my_theme1 +
  theme(legend.text = element_text(size = 15), axis.title.y = element_text(size = 14)) + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text = element_text(face = "italic"))+
  geom_signif( comparisons = list(c("cattle", "wildlife")), map_signif_level = TRUE)+ 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "none")





########################################################################################################################
################################################# Immune markers  ######################################################
########################################################################################################################


immune_meta <- data.frame(sample_data(dero_faecal))

# Distributions

hist(immune_meta$ig_G) 
hist(immune_meta$wbc_40x) 
hist((immune_meta$neutrophils/immune_meta$lymphocytes)) 
hist(immune_meta$lysozyme) 

# Make log transformations as needed
# Add column of log WBC, NL and Lysozyme

immune_meta <- immune_meta %>% mutate(log.WBC = log(wbc_40x), log.NL = log(neutrophils/lymphocytes), log.Lysozyme = log(lysozyme))

# Check distributions again

hist(immune_meta$log.WBC) 
hist((immune_meta$log.NL)) 
hist(immune_meta$log.Lysozyme) 

## Plots by feeding preference

P.IgG <- ggplot(immune_meta, aes(x = feeding_2c, y = ig_G, color = feeding_2c))+
  geom_boxplot() + theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = feeding_2c), width = 0.2) 

P.WBC <- ggplot(immune_meta, aes(x = feeding_2c, y = log.WBC, color = feeding_2c))+
  geom_boxplot()+theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = feeding_2c), width = 0.2) 

P.NL <- ggplot(immune_meta, aes(x = feeding_2c, y = log.NL, color = feeding_2c))+
  geom_boxplot()+theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = feeding_2c), width = 0.2) 

P.LYS <- ggplot(immune_meta, aes(x = feeding_2c, y = log.Lysozyme, color = feeding_2c))+
  geom_boxplot()+theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = feeding_2c), width = 0.2) 

ggarrange(P.IgG,P.WBC,P.NL, P.LYS, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom") + 
  theme_classic()+
  ggtitle("Immune markers between diet preferences")



## Plots by disturbance at roosting cave location

P.IgG1 <- ggplot(immune_meta, aes(x = disturbance_2C, y = ig_G, color = disturbance_2C))+
  geom_boxplot() + theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = disturbance_2C), width = 0.2) 

P.WBC1 <- ggplot(immune_meta, aes(x = disturbance_2C, y = log.WBC, color = disturbance_2C))+
  geom_boxplot()+theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = disturbance_2C), width = 0.2) 

P.NL1 <- ggplot(immune_meta, aes(x = disturbance_2C, y = log.NL, color = disturbance_2C))+
  geom_boxplot()+theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = disturbance_2C), width = 0.2) 

P.LYS1 <- ggplot(immune_meta, aes(x = disturbance_2C, y = log.Lysozyme, color = disturbance_2C))+
  geom_boxplot()+theme_classic()+
  stat_compare_means(method = "t.test")  +
  geom_jitter(aes(col = disturbance_2C), width = 0.2) 

ggarrange(P.IgG1,P.WBC1,P.NL1, P.LYS1, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom") + 
  theme_classic()+
  ggtitle("Immune markers between habitats")



# Linear regressions of IM and feeding preference

lm.igg <- lm(ig_G ~ dC13_VPDB + sex, data = immune_meta)
plot(lm.igg)
summary(lm.igg)

lm.wbc <- lm(log.WBC ~ dC13_VPDB + sex, data = immune_meta)
summary(lm.wbc) 

lm.nl <- lm(log.NL ~ dC13_VPDB + sex, data = immune_meta)
summary(lm.nl)

lm.lys <- lm(log.Lysozyme ~ dC13_VPDB + sex, data = immune_meta)
summary(lm.lys)





## Now investigate ASVs identified by Ancom against immune markers

X131daaddbb4eba63347cf2c984e775f1 # Helicobacter sp.
X29f5d7605cf4510c64c7b4fb1585949d # Edwardsiella sp.


## Add ASVs found at 60% in ancom
dero_tax <- tax_table(dero_faecal)

names <- rownames(dero_tax)
test<- as.matrix(cbind(dero_tax, names))
tax_table(dero_faecal) <- test # map back

my_ancom_asvs <- subset_taxa(dero_faecal, names=="131daaddbb4eba63347cf2c984e775f1"|
                               names=="29f5d7605cf4510c64c7b4fb1585949d")


my_ancom_asvs.df <- data.frame(otu_table(my_ancom_asvs))

# add +1 to all values in the dataframe, then we are able to compare abundances later properly with log transformation
my_ancom_asvs.df_1 <- my_ancom_asvs.df + 1  # first save into new object, double check is correct
my_ancom_asvs.df <- my_ancom_asvs.df_1      # then rename

my_ancom_asvs.df  = cbind(as(my_ancom_asvs.df, "data.frame"), as(tax_table(my_ancom_asvs)[rownames(my_ancom_asvs.df), ], "matrix"))

my_ancom_asvs.df<-data.frame(t(my_ancom_asvs.df))
head(my_ancom_asvs.df)
tail(my_ancom_asvs.df)

# check if rows are ordered similar in metadata
metadata<-data.frame(sample_data(my_ancom_asvs))
head(metadata)
tail(metadata)

# both metadata and asv-df are in the same order, thus we can bind together
metadata1<-cbind(metadata, my_ancom_asvs.df)

# Change colnames of ASV columns to correct and short names
colnames(metadata1)[colnames(metadata1) %in% c("X131daaddbb4eba63347cf2c984e775f1", "X29f5d7605cf4510c64c7b4fb1585949d")] <- c("Helicobacter", "Edwardsiella")

#convert to long format
data_long <- gather(metadata1, ASV_ID, Abundance,  "Helicobacter":"Edwardsiella", factor_key=TRUE)

table(data_long$ASV_ID) # check it worked

data_long$Abundance <- as.numeric(data_long$Abundance)
data_long$ASV_ID


# add immune info to microbiome sample_data
immune_meta

# check both files are in the same order
rownames(data_long)
rownames(immune_meta) <- immune_meta$feature_id
rownames(immune_meta)

# now add relevant columns to df
data_long$feeding_2c <- as.data.frame(immune_meta)$feeding_2c
data_long$ig_G <- as.data.frame(immune_meta)$ig_G
data_long$log.WBC <- as.data.frame(immune_meta)$log.WBC
data_long$log.NL <- as.data.frame(immune_meta)$log.NL
data_long$log.Lysozyme <- as.data.frame(immune_meta)$log.Lysozyme



# change order of asvs in a new df

ggplot(data_long, aes(x = disturbance_2C, y = Abundance))+
  geom_boxplot(outlier.shape = NA, aes(fill=disturbance_2C))+
  geom_jitter(width=0.05, size = 1, alpha = 0.6)+
  facet_wrap(~ASV_ID, scales = "free", ncol = 2)+
  scale_y_log10()+
  theme_bw(base_size=16)+
  stat_compare_means(method = "wilcox.test", size=4.8, aes(label = sprintf("Wilcoxon, p = %5.4f", as.numeric(..p.format..))))+
  labs(x= " ", y =expression(Log[10] * " ASV abundance")) +
  theme(legend.position=c(0.6,0.15))+
  theme(legend.text = element_text(size = 17))+ggtitle("B")+theme(strip.text.x = element_text(face = "bold.italic"))


ggplot(data_long, aes(x = feeding_2c, y = Abundance))+
  geom_boxplot(outlier.shape = NA, aes(fill=feeding_2c))+
  geom_jitter(width=0.05, size = 1, alpha = 0.6)+
  facet_wrap(~ASV_ID, scales = "free", ncol = 2)+
  scale_y_log10()+
  theme_bw(base_size=16)+
  stat_compare_means(method = "wilcox.test", size=4.8, aes(label = sprintf("Wilcoxon, p = %5.4f", as.numeric(..p.format..))))+
  labs(x= " ", y =expression(Log[10] * " ASV abundance")) +
  theme(legend.position=c(0.6,0.15))+
  theme(legend.text = element_text(size = 17))+ggtitle("B")+theme(strip.text.x = element_text(face = "bold.italic"))




helicobacter_df <- subset(data_long, ASV_ID =="Helicobacter")
hist(helicobacter_df$Abundance)
hist(log10(helicobacter_df$Abundance))



lm.igg.c <- lm(ig_G ~ log(helicobacter_df$"Abundance"), data = helicobacter_df)
summary(lm.igg.c)

lm.wbc.c <- lm(log.WBC ~ log(helicobacter_df$"Abundance"), data = helicobacter_df)
summary(lm.wbc.c)

lm.nl.c <- lm(log.NL ~ log(helicobacter_df$"Abundance"), data = helicobacter_df)
summary(lm.nl.c)

lm.lys.c <- lm(log.Lysozyme ~ log(helicobacter_df$"Abundance"), data = helicobacter_df)
summary(lm.lys.c)





edwardsiella_df <- subset(data_long, ASV_ID =="Edwardsiella")
hist(edwardsiella_df$Abundance)
hist(log10(edwardsiella_df$Abundance))



lm.igg.c <- lm(ig_G ~ log(edwardsiella_df$"Abundance"), data = edwardsiella_df)
summary(lm.igg.c)
plot(lm.igg.c)
confint(lm.igg.c) # 95% CI for the coefficients
predict(lm.igg.c, type="response") # predicted values
residuals(lm.igg.c, type="deviance") # residuals 


lm.wbc.c <- lm(log.WBC ~ log(edwardsiella_df$"Abundance"), data = edwardsiella_df)
summary(lm.wbc.c)

lm.nl.c <- lm(log.NL ~ log(edwardsiella_df$"Abundance"), data = edwardsiella_df)
summary(lm.nl.c)

lm.lys.c <- lm(log.Lysozyme ~ log(edwardsiella_df$"Abundance"), data = edwardsiella_df)
summary(lm.lys.c)






# create x-label since it should be partially italics

my_x_title <- expression(paste("log(Total Abundance) ", italic(" Edwardsiella sp.")))

# function to plot relevant values at top of the plot

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    theme_bw() +
    stat_smooth(method = "lm", col = "red") +
    xlab(my_x_title) +
    ylab("IgG") +
    theme(axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 14)) +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       " p =",signif(summary(fit)$coef[2,4], 2)))
}


#p <- ggplotRegression(lm(ig_G ~ log(edwardsiella_df$"Abundance"), data = edwardsiella_df))




### Linear regression ###


# first extract otu-data
otu <- otu_table(dero_faecal)

# add +1 to all values in the dataframe, then we are able to compare abundances later properly 
otu_1 <- otu + 1

# transpose, because we need samples as rows, as in the metadata file
otu_trans <- t(as.data.frame(otu_1))
otu_trans <- as.data.frame(otu_trans)

# add immune info to microbiome sample_data
immune_meta

# check both files are in the same order
rownames(otu_trans)
rownames(immune_meta) <- immune_meta$feature_id
rownames(immune_meta)

# now add relevant columns to df
otu_trans$feeding_2c <- as.data.frame(immune_meta)$feeding_2c
otu_trans$ig_G <- as.data.frame(immune_meta)$ig_G
otu_trans$log.WBC <- as.data.frame(immune_meta)$log.WBC
otu_trans$log.NL <- as.data.frame(immune_meta)$log.NL
otu_trans$log.Lysozyme <- as.data.frame(immune_meta)$log.Lysozyme

# create x-label since it should be partially italics
my_x_title <- expression(paste("log(Total Abundance) ", italic(" Edwardsiella sp.")))

my_x_title2 <- expression(paste("log(Total Abundance) ", italic(" Helicobacter sp.")))


### 1) Helicobacter sp. vs IM

hist(otu_trans$"131daaddbb4eba63347cf2c984e775f1")
hist(log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"))

# Linear regressions

lm.igg.c <- lm(ig_G ~ log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), data = otu_trans)
summary(lm.igg.c)

lm.wbc.c <- lm(log.WBC ~ log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), data = otu_trans)
summary(lm.wbc.c)

lm.nl.c <- lm(log.NL ~ log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), data = otu_trans)
summary(lm.nl.c)

lm.lys.c <- lm(log.Lysozyme ~ log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), data = otu_trans)
summary(lm.lys.c)



### Edwardsiella sp. vs IM

otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"
log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d")

# Linear regressions

lm.igg.c <- lm(ig_G ~ log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), data = otu_trans)
summary(lm.igg.c) # tendency to higher IgG levels and high abundance of Edwardsiella sp.

lm.wbc.c <- lm(log.WBC ~ log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), data = otu_trans)
summary(lm.wbc.c)

lm.nl.c <- lm(log.NL ~ log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), data = otu_trans)
summary(lm.nl.c)

lm.lys.c <- lm(log.Lysozyme ~ log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), data = otu_trans)
summary(lm.lys.c)



## regression figures

immune_1 <- ggplot(otu_trans, aes(x = log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), y = ig_G)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=0.6, size= 4, label = "adj. R2=0.065, p=0.05")+
  labs(x=my_x_title, y="IgG")+
  theme_bw(base_size = 14)+
  my_theme1

immune_2 <- ggplot(otu_trans, aes(x = log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), y = log.WBC)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=-2.3, size= 3.5, label = "adj. R2=0.024, p=0.16")+
  labs(x="", y="log10 WBC")+
  theme_bw(base_size = 11)+
  my_theme1

immune_3 <- ggplot(otu_trans, aes(x = log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), y = log.NL)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=-1.2, size= 3.5, label = "adj. R2=-0.001, p=0.45")+
  labs(x=my_x_title, y="log10 NL")+
  theme_bw(base_size = 11)+
  my_theme1

immune_4 <- ggplot(otu_trans, aes(x = log(otu_trans$"29f5d7605cf4510c64c7b4fb1585949d"), y = log.Lysozyme)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=-0.3, size= 3.5, label = "adj. R2=0.001, p=0.31")+
  labs(x="", y="log10 lysozyme")+
  theme_bw(base_size = 11)+
  my_theme1

immune_5 <- ggplot(otu_trans, aes(x = log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), y = ig_G)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=5.7, y=0.6, size= 3.5, label = "adj. R2=-0.016, p=0.57")+
  labs(x="", y="IgG")+
  theme_bw(base_size = 11)+
  my_theme1

immune_6 <- ggplot(otu_trans, aes(x = log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), y = log.WBC)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=-2.3, size= 3.5, label = "adj. R2=0.015, p=0.2")+
  labs(x="", y="log10 WBC")+
  theme_bw(base_size = 11)+
  my_theme1

immune_7 <- ggplot(otu_trans, aes(x = log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), y = log.NL)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=-1.2, size= 3.5, label = "adj. R2=-0.006, p=0.26")+
  labs(x=my_x_title2, y="log10 NL")+
  theme_bw(base_size = 11)+
  my_theme1

immune_8 <- ggplot(otu_trans, aes(x = log(otu_trans$"131daaddbb4eba63347cf2c984e775f1"), y = log.Lysozyme)) + 
  geom_point() +
  #stat_smooth(method = "lm", col = "black", linetype = "dashed") +
  geom_text(x=6.2, y=-0.3, size= 3.5, label = "adj. R2=0.013, p=0.21")+
  labs(x="", y="log10 lysozyme")+
  theme_bw(base_size = 11)+
  my_theme1


immune_a <- ggarrange(immune_6, immune_7, immune_8,immune_5,  labels = c("a)"), ncol=4)
immune_b <- ggarrange(immune_2, immune_3, immune_4, labels = c("b)"), ncol=4)



#################
### Figure S7 ###
#################

ggarrange(immune_a, immune_b, ncol=1)


################
### Figure 2 ###
################

ggarrange(ancom_fig_a, ancom_fig_b, immune_1, labels= c("a", "b", "c"), ncol=2, nrow=2, widths = c(1,1.3))





## Most abundant bacterial taxa

## Comparison feeding preferences

dero_merged_feeding <- merge_samples(dero_faecal,"feeding_2c", fun = mean) # merge samples by storage medium
dero_merged_feeding <- tax_glom(dero_merged_feeding, taxrank = "Family")   # agglomerate taxa to Genus level (all taxa that belong to same genus become one ASV)
dero_merged_feeding # now only one sample representing feeding preferences


# change to compositional format (i.e., relative abundance)
dero_merged_feeding <- microbiome::transform(dero_merged_feeding, transform = "compositional") 

composition_df <- psmelt(dero_merged_feeding)    # psmelt() 'melts' all information in a phyloseq object into a long dataframe                   

# same, but on Family level
composition_summary <- composition_df %>%              # take the dataframe 'composition_df"
  group_by(Sample, Family) %>%                        # group it by Sample and Phylum (NOTE: change 'Phylum' to another taxonomic rank if you want to summarize)
  summarise(Relative_abundance = sum(Abundance)) %>%  # summarize the sum abundance [of each family per sample]
  arrange(Sample, desc(Relative_abundance)) %>%       # arrange by (descending) relative abundance
  mutate(Rank=rank(desc(Relative_abundance)))         # add a column with rank

################
### Table S3 ###
################

composition_summary<-subset(composition_summary, Rank <=10) # only keep the top 10

# view in table format
datatable(composition_summary)



## Comparison roosting cave location

dero_merged_habitat <- merge_samples(dero_faecal,"disturbance_2C", fun = mean) # merge samples by storage medium
dero_merged_habitat <- tax_glom(dero_merged_habitat, taxrank = "Family")   # agglomerate taxa to Genus level (all taxa that belong to same genus become one ASV)
dero_merged_habitat # now only one sample representing feeding preferences


# change to compositional format (i.e., relative abundance)
dero_merged_habitat <- microbiome::transform(dero_merged_habitat, transform = "compositional") 

composition_df <- psmelt(dero_merged_habitat)    # psmelt() 'melts' all information in a phyloseq object into a long dataframe                   

# same, but on Family level
composition_summary <- composition_df %>%              # take the dataframe 'composition_df"
  group_by(Sample, Family) %>%                        # group it by Sample and Phylum (NOTE: change 'Phylum' to another taxonomic rank if you want to summarize)
  summarise(Relative_abundance = sum(Abundance)) %>%  # summarize the sum abundance [of each family per sample]
  arrange(Sample, desc(Relative_abundance)) %>%       # arrange by (descending) relative abundance
  mutate(Rank=rank(desc(Relative_abundance)))         # add a column with rank

################
### Table S3 ###
################

composition_summary<-subset(composition_summary, Rank <=10) # only keep the top 10

# view table
datatable(composition_summary)




## Most abundant bacterial taxa: ASV prevalence vs. abundance

# Separate by feeding preference & roosting cave location

dero_pristine_unrare <-  subset_samples(dero_faecal, disturbance_2C == "p")
dero_disturbed_unrare <-  subset_samples(dero_faecal, disturbance_2C == "d")

dero_wildlife_unrare <- subset_samples(dero_faecal, feeding_2c == "wildlife")
dero_cattle_unrare <- subset_samples(dero_faecal, feeding_2c == "cattle")


###### wildlife blood ######

prevalence_df_wildife <- data.frame(microbiome::prevalence(dero_wildlife_unrare)) # a function called 'prevalence' from the package 'microbiome' that calculates prevalence for each ASV, and we turn this into a dataframe to work with

prevalence_df_wildife$Abundance<-taxa_sums(dero_wildlife_unrare)  # add a column with abundance
prevalence_df_wildife$Relative_abundance<-prevalence_df_wildife$Abundance/sum(prevalence_df_wildife$Abundance)  # calculate relative abundance
prevalence_df_wildife$Phylum<-data.frame(tax_table(dero_wildlife_unrare))$Phylum  # add phylum as a column
prevalence_df_wildife$Class<-data.frame(tax_table(dero_wildlife_unrare))$Class  # add class as a column
prevalence_df_wildife$Family<-data.frame(tax_table(dero_wildlife_unrare))$Family  # add family as a column
prevalence_df_wildife$Genus<-data.frame(tax_table(dero_wildlife_unrare))$Genus  # add genus as a column
prevalence_df_wildife$Species<-data.frame(tax_table(dero_wildlife_unrare))$Species  # add species as a column

names(prevalence_df_wildife)[1]<-"Prevalence"  # change name of first column to 'Prevalence'



# colourvector
colorlist <- c("#968a00",
               "#6e58e0",
               "#74ce3b",
               "#593bb5",
               "#019115",
               "#ed60e0",
               "#7ea100",
               "#0262e0",
               "#f9ab0b",
               "#0183fe",
               "#f7bc4d",
               "#85219f",
               "#00aa5a",
               "#bf1ea7",
               "#7cda8c",
               "#ef1485",
               "#009052",
               "#c4007c",
               "#01cdb8",
               "#ff4a4d",
               "#00a49c",
               "#ff4879",
               "#006132",
               "#ff61ab",
               "#b1d26e",
               "#9d8bff",
               "#c6cd64",
               "#624098",
               "#ffa54a",
               "#015499",
               "#be6c00",
               "#9da5ff",
               "#b84700",
               "#80baff",
               "#a71d00",
               "#029bb6",
               "#ff6e48",
               "#008772",
               "#c60045",
               "#007b68",
               "#9b0d74",
               "#dbc669",
               "#eba2ff",
               "#5e6900",
               "#ffa8f4",
               "#7d7200",
               "#dcb9f7",
               "#9f7200",
               "#564c7c",
               "#ff9463",
               "#8a7aad",
               "#a70920",
               "#82925d",
               "#ff9bd0",
               "#65500e",
               "#ff99a9",
               "#5d5a27",
               "#942c43",
               "#c5b67f",
               "#813b58",
               "#ffb493",
               "#953010",
               "#b57963",
               "#784f00",
               "#7a452b")


# scatterplot of prevalence against abundance
p.wildlife <- ggplot(prevalence_df_wildife, aes(x = Prevalence, y = Relative_abundance, fill = Phylum, label = Family)) + # fill colour by Phylum, add text labesl 
  geom_point(aes(size = Abundance), pch = 21, alpha = 0.7) +  # size by abundance, make slightly transparent (alpha = 0.7)
  scale_size(range=c(3,8)) +  # set limits to the point size (between 3 and 8) so that points aren't too small
  scale_fill_manual(values = colorlist) +                 
  guides(fill = guide_legend(override.aes = list(size = 4))) +  # this makes the legend symbols bigger
  geom_text_repel( data= subset(prevalence_df_wildife, Relative_abundance > 0.03)) + # add Family labels to ASVs with relative abundance over 3%
  theme_bw(base_size = 12) +
  ylab("relative abundance of ASVs")+
  ggtitle("mainly wildlife-blood feeding")



###### cattle blood ######

prevalence_df_cattle <- data.frame(microbiome::prevalence(dero_cattle_unrare)) # a function called 'prevalence' from the package 'microbiome' that calculates prevalence for each ASV, and we turn this into a dataframe to work with

prevalence_df_cattle$Abundance<-taxa_sums(dero_cattle_unrare)  # add a column with abundance
prevalence_df_cattle$Relative_abundance<-prevalence_df_cattle$Abundance/sum(prevalence_df_cattle$Abundance)  # calculate relative abundance
prevalence_df_cattle$Phylum<-data.frame(tax_table(dero_cattle_unrare))$Phylum  # add phylum as a column
prevalence_df_cattle$Class<-data.frame(tax_table(dero_cattle_unrare))$Class  # add class as a column
prevalence_df_cattle$Family<-data.frame(tax_table(dero_cattle_unrare))$Family  # add family as a column
prevalence_df_cattle$Genus<-data.frame(tax_table(dero_cattle_unrare))$Genus  # add genus as a column
prevalence_df_cattle$Species<-data.frame(tax_table(dero_cattle_unrare))$Species  # add species as a column

names(prevalence_df_cattle)[1]<-"Prevalence"  # change name of first column to 'Prevalence'


# scatterplot of prevalence against abundance
p.cattle <- ggplot(prevalence_df_cattle, aes(x = Prevalence, y = Relative_abundance, fill = Phylum, label = Family)) + # fill colour by Phylum, add text labesl 
  geom_point(aes(size = Abundance), pch = 21, alpha = 0.7) +  # size by abundance, make slightly transparent (alpha = 0.7)
  scale_size(range=c(3,8)) +  # set limits to the point size (between 3 and 8) so that points aren't too small
  scale_fill_manual(values = colorlist) +                 
  guides(fill = guide_legend(override.aes = list(size = 4))) +  # this makes the legend symbols bigger
  geom_text_repel( data= subset(prevalence_df_cattle, Relative_abundance > 0.03)) + # add Family labels to ASVs with relative abundance over 3%
  theme_bw(base_size = 12) +
  ylab("relative abundance of ASVs")+
  ggtitle("mainly cattle-blood feeding")




###### roosting cave in pristine forest ######

prevalence_df_pristine <- data.frame(microbiome::prevalence(dero_pristine_unrare)) # a function called 'prevalence' from the package 'microbiome' that calculates prevalence for each ASV, and we turn this into a dataframe to work with

prevalence_df_pristine$Abundance<-taxa_sums(dero_pristine_unrare)  # add a column with abundance
prevalence_df_pristine$Relative_abundance<-prevalence_df_pristine$Abundance/sum(prevalence_df_pristine$Abundance)  # calculate relative abundance
prevalence_df_pristine$Phylum<-data.frame(tax_table(dero_pristine_unrare))$Phylum  # add phylum as a column
prevalence_df_pristine$Class<-data.frame(tax_table(dero_pristine_unrare))$Class  # add class as a column
prevalence_df_pristine$Family<-data.frame(tax_table(dero_pristine_unrare))$Family  # add family as a column
prevalence_df_pristine$Genus<-data.frame(tax_table(dero_pristine_unrare))$Genus  # add genus as a column
prevalence_df_pristine$Species<-data.frame(tax_table(dero_pristine_unrare))$Species  # add species as a column

names(prevalence_df_pristine)[1]<-"Prevalence"  # change name of first column to 'Prevalence'


# scatterplot of prevalence against abundance
p.pristine <- ggplot(prevalence_df_pristine, aes(x = Prevalence, y = Relative_abundance, fill = Phylum, label = Family)) + # fill colour by Phylum, add text labesl 
  geom_point(aes(size = Abundance), pch = 21, alpha = 0.7) +  # size by abundance, make slightly transparent (alpha = 0.7)
  scale_size(range=c(3,8)) +  # set limits to the point size (between 3 and 8) so that points aren't too small
  scale_fill_manual(values = colorlist) +                 
  guides(fill = guide_legend(override.aes = list(size = 4))) +  # this makes the legend symbols bigger
  geom_text_repel( data= subset(prevalence_df_pristine, Relative_abundance > 0.03)) + # add Family labels to ASVs with relative abundance over 3%
  theme_bw(base_size = 12) +
  ylab("relative abundance of ASVs")+
  ggtitle("roosting cave in pristine forest")



###### roosting cave in disturbed forest ######

prevalence_df_disturbed <- data.frame(microbiome::prevalence(dero_disturbed_unrare)) # a function called 'prevalence' from the package 'microbiome' that calculates prevalence for each ASV, and we turn this into a dataframe to work with

prevalence_df_disturbed$Abundance<-taxa_sums(dero_disturbed_unrare)  # add a column with abundance
prevalence_df_disturbed$Relative_abundance<-prevalence_df_disturbed$Abundance/sum(prevalence_df_disturbed$Abundance)  # calculate relative abundance
prevalence_df_disturbed$Phylum<-data.frame(tax_table(dero_disturbed_unrare))$Phylum  # add phylum as a column
prevalence_df_disturbed$Class<-data.frame(tax_table(dero_disturbed_unrare))$Class  # add class as a column
prevalence_df_disturbed$Family<-data.frame(tax_table(dero_disturbed_unrare))$Family  # add family as a column
prevalence_df_disturbed$Genus<-data.frame(tax_table(dero_disturbed_unrare))$Genus  # add genus as a column
prevalence_df_disturbed$Species<-data.frame(tax_table(dero_disturbed_unrare))$Species  # add species as a column

names(prevalence_df_disturbed)[1]<-"Prevalence"  # change name of first column to 'Prevalence'


# scatterplot of prevalence against abundance
p.disturbed <- ggplot(prevalence_df_disturbed, aes(x = Prevalence, y = Relative_abundance, fill = Phylum, label = Family)) + # fill colour by Phylum, add text labesl 
  geom_point(aes(size = Abundance), pch = 21, alpha = 0.7) +  # size by abundance, make slightly transparent (alpha = 0.7)
  scale_size(range=c(3,8)) +  # set limits to the point size (between 3 and 8) so that points aren't too small
  scale_fill_manual(values = colorlist) +                 
  guides(fill = guide_legend(override.aes = list(size = 4))) +  # this makes the legend symbols bigger
  geom_text_repel( data= subset(prevalence_df_disturbed, Relative_abundance > 0.03)) + # add Family labels to ASVs with relative abundance over 3%
  theme_bw(base_size = 12) +
  ylab("relative abundance of ASVs")+
  ggtitle("roosting cave in disturbed forest")



###################
#### Figure S7 ####
###################

ggarrange(p.wildlife, p.cattle, p.pristine, p.disturbed, nrow=2, ncol=2, labels = c("a","b","c","d"), common.legend = TRUE, legend = "right")


