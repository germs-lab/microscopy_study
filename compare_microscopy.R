
setwd("~/Box Sync/Github/germs-lab/microscopy_study/")
##########   Load library   ############
library(ape)
library(ggplot2)
library(grid)
library(plyr)
library(reshape)
library(reshape2)
library(gridExtra)
library(corrplot)
library(vegan)
library(phyloseq)

#//////////////////////////////////////////////////////////////////////////////////////////#
##########  Read data, USE saved object ############
#//////////////////////////////////////////////////////////////////////////////////////////#
# physeq <- readRDS("water-project-object_physeq_before_normalization.rds")
# #Change OTU name
# otu_id <- data.frame(paste("OTU-", 1:ntaxa(physeq), sep=""),taxa_names(physeq))
# colnames(otu_id) <- c("OTU name", "Qiime assigned OTU name")
# write.csv(otu_id, "internal_info.csv")
# taxa_names(physeq) <- paste("OTU-", 1:ntaxa(physeq), sep="")

# #two group
# sample_data(physeq)$two_group <- sample_data(physeq)$new_group
# sample_data(physeq)$two_group <- gsub("PlimHC", "Plim", sample_data(physeq)$two_group)
# sample_data(physeq)$two_group <- gsub("PlimLC", "Plim", sample_data(physeq)$two_group)
# sample_data(physeq)$two_group <- gsub("NlimHC", "Nlim", sample_data(physeq)$two_group)
# sample_data(physeq)$two_group <- gsub("NlimLC", "Nlim", sample_data(physeq)$two_group)

# # Rarefy
# set.seed(1)
# rare <- rarefy_even_depth(physeq)
# saveRDS(rare, "water-project-object-physeq-rarefied.rds")
# # relative abundance
# GPr = transform_sample_counts(rare, function(x) x/sum(x))
# GPf = filter_taxa(GPr, function(x) var(x) > 1e-05, TRUE)

# # microscopy data
# phyto <- read.csv("phytoplankton ds FINAL.csv")


# # this area to calculate percent of eukaryotic autotrophs in microscoy method
# phyto64 <- subset(phyto, sampleid %in% sample_names(rare))
# sum(colSums(phyto64[,4:ncol(phyto64)]))
# onlycyno <- phyto64[,grepl("Cyanophyta", names(phyto64))]
# 1- sum(colSums(onlycyno)) / sum(colSums(phyto64[,4:ncol(phyto64)]))


# #only cyanobacteria
# onlycyno <- phyto[,grepl("Cyanophyta", names(phyto))]

# optical <- cbind(sampleid = phyto$sampleid, onlycyno)
# #64 samples
# optical64 <- subset(optical, sampleid %in% sample_names(rare))

# write.csv(optical64, file="microscopy_data.csv", row.names = F)

# colnames(optical64)[1] <- "sample ID"
# write.csv(optical64, file="Supplementary_table2.csv", row.names = F)


rare <- readRDS("water-project-object-physeq-rarefied.rds")
optical64 <- read.csv("microscopy_data.csv")

# relative abundance
GPr = transform_sample_counts(rare, function(x) x/sum(x))
GPf = filter_taxa(GPr, function(x) var(x) > 1e-05, TRUE)


# supplementary table 1
b <- colnames(optical64)[2:28]
physeqdf <- psmelt(GPr)
#head(physeqdf)
a <- as.character(unique(physeqdf$Rank6))
c <- a[gsub("g__","", a) %in% gsub("Cyanophyta_","", b)]
rlist <- list(a,b,c)
attributes(rlist) =  list(names= names(rlist))
mem <- ldply(rlist, rbind)
tmem <- t(mem)
no_g_tmem <- gsub("g__","",tmem)
no_na_tmem <- subset(no_g_tmem, no_g_tmem[,1] != "")
colnames(no_na_tmem) <- c("Genera identified through 16S rRNA sequencing","Genera identified through microscopy","Genera identified through both methods")
write.csv(no_na_tmem, file="supplementary_table1.csv")


# # one_phylum <- subset_taxa(GPr, Rank2 == "p__Cyanobacteria")
# physeqdf <- psmelt(one_phylum)
# unique(physeqdf$Rank6)
# sample_sums(one_phylum)
# range(sample_sums(one_phylum))

# opti_meta <- read.csv("~/Box Sync/2016/6June/16s/meta_w_new_group.csv")
# se <- data.frame(sample_sums(one_phylum))
# op <- data.frame(opti_meta$sampleid, opti_meta$cyanophyta.ra)
# merged <- merge(se, op, by.x = "row.names", by.y = "opti_meta.sampleid")
# ggplot(merged, aes(x=sample_sums.one_phylum., y= opti_meta.cyanophyta.ra))+geom_point()
# subset(merged, sample_sums.one_phylum. > 0.4 & opti_meta.cyanophyta.ra < 0.5)
# subset(merged, sample_sums.one_phylum. < 0.1 & opti_meta.cyanophyta.ra > 0.9)

# out <- c("2014223012", "2014188002","2014202005","2014202012")
# subset(opti_meta, sampleid %in% out)

# taxa_names(one_phylum)

# #major phylum
# major_phylum <- subset_taxa(GPr, Rank2 == "p__Acidobacteria" |Rank2 == "p__Actinobacteria" |Rank2 == "p__Armatimonadetes"|Rank2 == "p__Bacteroidetes"|Rank2 == "p__Chlorobi"|Rank2 == "p__Cyanobacteria"|Rank2 == "p__Planctomycetes"|Rank2 == "p__Proteobacteria"|Rank2 == "p__Verrucomicrobia" )
# range(sample_sums(major_phylum))


# major_phylum <- subset_taxa(rare, Rank2 == "p__Cyanobacteria")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Actinobacteria")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Proteobacteria")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Bacteroidetes")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Acidobacteria")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Armatimonadetes")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Chlorobi")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Planctomycetes")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
# major_phylum <- subset_taxa(rare, Rank2 == "p__Verrucomicrobia")
# sum(sample_sums(major_phylum)) / sum(sample_sums(rare))



#//////////////////////////////////////////////////////////////////////////////////////////#
##########  Figure 1 diversity and richness   ##########  
#//////////////////////////////////////////////////////////////////////////////////////////#
#usre rarefyed data

#diversity for sequencing
otu <- t(as.data.frame(otu_table(rare)))

di_16 <- diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))


#diversity for optical
cyano_optical <- optical64
rownames(cyano_optical) <- cyano_optical$sampleid
cyano_optical$sampleid <- NULL
di_mi <- diversity(cyano_optical[,-c(1)] , index = "shannon", MARGIN = 1, base = exp(1))

#merge and plot
merged <- merge(di_16, di_mi,by.x="row.names", by.y="row.names")
mean(merged$x)
range(merged$x)
mean(merged$y)
range(merged$y)
colnames(merged)=c("sample_id", "Sequencing Method","Microscopy Method")
melted <- melt(merged, id = c("sample_id"))

pdf("~/Box Sync/2017/7July/water/figures_tables/Figure1-diversity.pdf", width=6, height=6)
ggplot(melted, aes(x=variable, y=value))+geom_boxplot()+labs(x="Method",y="Shannon's Diversity")+theme_bw()
dev.off()





#//////////////////////////////////////////////////////////////////////////////////////////#
##########     FIgure 2  Compare All Cyanobacteria       ##########  
#//////////////////////////////////////////////////////////////////////////////////////////#

# we relative abundance

# all Cyanobacteria
sco <- data.frame(optical64$sampleid, rowSums(optical64[,2:ncol(optical64)]))
colnames(sco) <- c("sampleid", "microcystis_scopy")
range(sco$microcystis_scopy)
mean(sco$microcystis_scopy)
#sequncing this is toxic
#sub <- subset_taxa(rare, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena"| Rank6 == "g__Dolichospermum")

#relative abundance
sub <- subset_taxa(GPr, Rank2 == "p__Cyanobacteria")
sample_data(GPr)[c("2014188003"),]
otu_table(sub)[,c("2014188003")]

one_sample <- prune_samples(sample_names(sub)=="2014188003",sub)
glom <- tax_glom(one_sample, "Rank6")
48.6/62.7
tax_table(GPr)
seq16s <- data.frame(sample_names(sub),sample_sums(sub))
colnames(seq16s) <- c("sampleid", "microcystis_16s")

merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
cor(merged$microcystis_16s, merged$microcystis_scopy)


## Get fit, and make a variable with labels for points outside conf. interval
fit <- lm(merged$microcystis_scopy ~ merged$microcystis_16s)
dat <- predict(fit, interval="confidence")
merged$color <- ifelse(merged$microcystis_scopy < dat[,"upr"] & merged$microcystis_scopy > dat[,"lwr"], "black", "red")




pdf("~/Box Sync/2017/7July/water/figures_tables/Figure2_cyano_compare.pdf", width=6, height=6)
ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point(colour = merged$color)+labs(x="Cyanobacteria 16S rRNA gene relative abundance", y="Cyanobacteria biomass per liter of water (mg/L)")+theme_bw()
dev.off()
cor.test(merged$microcystis_16s, merged$microcystis_scopy)
outlier <- subset(merged, microcystis_16s > 0.4 & microcystis_scopy <50)
nrow(subset(merged,color=="red"))

sample_data(GPr)[outlier$sampleid,]
sample_data(GPr)[c("2014188002","2014202005","2014223012"),]
     # sampleid microcystis_16s microcystis_scopy
# 2  2014188002       0.4309904          1.160961
# 34 2014202005       0.5060703         11.864625
# 61 2014223012       0.4846645          1.682105

outlier <- subset(merged, microcystis_16s < 0.3 & microcystis_scopy > 350)
sample_data(GPr)[c("2014210001"),]
#     sampleid microcystis_16s microcystis_scopy
#48 2014210001       0.2690096          365.6038


s_data<- as.data.frame(sample_data(GPr))
re_merged <- merge(merged, s_data, by.x = "sampleid", by.y = "row.names")
cor(re_merged$microcystis_16s, re_merged$cyanophyta.ra)
#ggplot(re_merged, aes(microcystis_16s, cyanophyta.ra))+ geom_smooth(method="lm")+geom_point()+labs(x="Cyanobacteria 16S gene relative abundance", y="Cyanobacteria biomass (mg/L)")
note <- c("2014188002","2014202005","2014223012")
subset(re_merged, sampleid %in% note)


subset(re_merged,microcystis_16s < 0.1 & cyanophyta.ra > 0.9 )

range(re_merged$microcystis_scopy)






#//////////////////////////////////////////////////////////////////////////////////////////#
##########       Overall toxic genus       ##########  
#//////////////////////////////////////////////////////////////////////////////////////////#
# how many toxic genus in the sample?
#only cyano
#now, toxic genus
onlycyno <- optical64
toxic_optical <- onlycyno[,grepl("Microcystis|Planktothrix|Anabaena|Cylindrospermopsis|Chroococcus|Pseudanabaena|Pseudanabaenaceae|Synechococcaceae", names(onlycyno))]
#this is overall ratio
sum(rowSums(toxic_optical)) / sum(rowSums(onlycyno[,2:ncol(onlycyno)]))
range(rowSums(toxic_optical)/rowSums(onlycyno[,2:ncol(onlycyno)]))


#sequncing
sub <- subset_taxa(rare, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena"| Rank6 == "g__Dolichospermum")

sum(rowSums(t(otu_table(sub)))) / sum(rowSums(t(otu_table(rare))))
range(  rowSums(t(otu_table(sub))) / rowSums(t(otu_table(rare)))  )



# # ##### below now used

# diver <- readRDS("~/Box Sync/2017/1January/water/water64/water-project-object_physeq_before_normalization.rds")
# #two group
# sample_data(diver)$two_group <- sample_data(diver)$new_group
# sample_data(diver)$two_group <- gsub("PlimHC", "Plim", sample_data(diver)$two_group)
# sample_data(diver)$two_group <- gsub("PlimLC", "Plim", sample_data(diver)$two_group)
# sample_data(diver)$two_group <- gsub("NlimHC", "Nlim", sample_data(diver)$two_group)
# sample_data(diver)$two_group <- gsub("NlimLC", "Nlim", sample_data(diver)$two_group)

# otu <- t(as.data.frame(otu_table(diver)))
# sa <-sample_data(diver)
# #dim(otu)
# di_16 <- diversity(otu, index = "shannon", MARGIN = 1, base = exp(1))
# #length(di_16)
# cyano_optical <- optical
# rownames(cyano_optical) <- cyano_optical$sampleid
# cyano_optical$sampleid <- NULL
# di_mi <- diversity(cyano_optical[,-c(1)] , index = "shannon", MARGIN = 1, base = exp(1))
# #length(di_mi)
# merged <- merge(di_16, di_mi,by.x="row.names", by.y="row.names")
# #colnames(merged)
# pdf("Figure2a.pdf", width=6, height=6)
# ggplot(merged, aes(x=x , y=y))+ geom_smooth(method="lm")+geom_point()+labs(x="Sequencing Method (16S)", y="Optical Method")
# dev.off()
# merged_sa <- merge(merged, sa, by.x="Row.names", by.y="row.names")
# merged_pl <- merged_sa[,c(2,3,21)]
# melted <- melt(merged_pl, id=c("two_group" ))
# mean(subset(melted, label == "Sequencing Method")$value)
# mean(subset(melted, label == "Optical Method")$value)
# #ggplot(merged_pl, aes(x=two_group, y=x))+geom_boxplot()+labs(x="Sequencing Method (16S)", y="Richness")
# fit <- aov(x ~ two_group, data=merged_pl)
# #plot(fit)
# summary(fit)

# fit <- aov(y ~ two_group, data=merged_pl)
# #plot(fit)
# summary(fit)

# melted <- melt(merged_pl, id=c("two_group" ))
# melted$label <- factor(melted$variable, labels=c("Sequencing Method","Optical Method"))
# pdf("Figure2b.pdf", width=6, height=6)
# ggplot(melted, aes(y=value, x=two_group))+geom_boxplot()+facet_grid( ~label )+labs(x="", y="Shanon's diversity")
# dev.off()

# #diversity contains richness and evenness
# #richness presence/absence
# physeq <- readRDS("~/Box Sync/2017/1January/water/water64/water-project-object_physeq_before_normalization.rds")
# #two group
# sample_data(physeq)$two_group <- sample_data(physeq)$new_group
# sample_data(physeq)$two_group <- gsub("PlimHC", "Plim", sample_data(physeq)$two_group)
# sample_data(physeq)$two_group <- gsub("PlimLC", "Plim", sample_data(physeq)$two_group)
# sample_data(physeq)$two_group <- gsub("NlimHC", "Nlim", sample_data(physeq)$two_group)
# sample_data(physeq)$two_group <- gsub("NlimLC", "Nlim", sample_data(physeq)$two_group)
# otu <- as.data.frame(otu_table(physeq))
# library(vegan)
# rare <- rarefy(t(otu), 10000)
# sam <- sample_data(physeq)
# merged <- merge(rare, sam,by.x="row.names", by.y="row.names")
# pdf("Figure2c.pdf", width=6, height=6)
# ggplot(merged, aes(x=two_group, y = x ))+geom_boxplot()+labs(x="Sequencing Method (16S)", y="Richness")
# dev.off()

# #colnames(merged)
# pdf("Figure2a.pdf", width=6, height=6)
# ggplot(merged, aes(x=x , y=y))+ geom_smooth(method="lm")+geom_point()+labs(x="Sequencing Method (16S)", y="Optical Method")
# dev.off()
# merged_sa <- merge(merged, sa, by.x="Row.names", by.y="row.names")
# merged_pl <- merged_sa[,c(2,3,21)]
# melted <- melt(merged_pl, id=c("two_group" ))
# mean(subset(melted, label == "Sequencing Method")$value)
# mean(subset(melted, label == "Optical Method")$value)
# #ggplot(merged_pl, aes(x=two_group, y=x))+geom_boxplot()+labs(x="Sequencing Method (16S)", y="Richness")
# fit <- aov(x ~ two_group, data=merged_pl)
# #plot(fit)
# summary(fit)

# fit <- aov(y ~ two_group, data=merged_pl)
# #plot(fit)
# summary(fit)

# melted <- melt(merged_pl, id=c("two_group" ))
# melted$label <- factor(melted$variable, labels=c("Sequencing Method","Optical Method"))
# pdf("Figure2b.pdf", width=6, height=6)
# ggplot(melted, aes(y=value, x=two_group))+geom_boxplot()+facet_grid( ~label )+labs(x="", y="Shanon's diversity")
# dev.off()


# #evenness how equal abundance 

#//////////////////////////////////////////////////////////////////////////////////////////#
##########        Figure  each toxic genus     ##########  
#//////////////////////////////////////////////////////////////////////////////////////////#
##########  overall abundance ##########  
# 1. Microcystis
# how many microcystis among cyanobacteria?
#sequencing

# all cyanobacteria
only_cyano  <- subset_taxa(rare, Rank2 == "p__Cyanobacteria")


#only microcystis
only_microcystis <- subset_taxa(rare,  Rank6 == "g__Microcystis" )
sum(rowSums(otu_table(only_microcystis))) / sum(rowSums(otu_table(only_cyano)))
#range
relative_one_genus = subset_taxa(GPr,  Rank6 == "g__Microcystis" )
range(sample_sums(relative_one_genus))

optical64[1:10,1:10]
sum(optical64[,grepl("Microcystis", names(optical64))]) / sum(rowSums(optical64[2:ncol(optical64)]))
range(optical64[,grepl("Microcystis", names(optical64))])

optical64[1:10,1:10]
sum(optical64[,grepl("Microcystis", names(optical64))]) / sum(rowSums(optical64[2:ncol(optical64)]))
range(optical64[,grepl("Microcystis", names(optical64))])

# 2. Planktothrix
only_one_genus <- subset_taxa(rare,  Rank6 == "g__Planktothrix" )
sum(rowSums(otu_table(only_one_genus))) / sum(rowSums(otu_table(only_cyano)))
#range
relative_one_genus = subset_taxa(GPr,  Rank6 == "g__Planktothrix" )
range(sample_sums(relative_one_genus))

# 3. Dolichospermum (Anabaena): 
only_one_genus <- subset_taxa(rare,  Rank6 == "g__Dolichospermum" | Rank6 == "g__Anabaena" )
sum(rowSums(otu_table(only_one_genus))) / sum(rowSums(otu_table(only_cyano)))


optical64
sum(optical64[,grepl("Anabaena", names(optical64))]) / sum(rowSums(optical64[2:ncol(optical)]))
range(optical64[,grepl("Anabaena", names(optical64))])

#4 Cylindrospermopsis
only_one_genus <- subset_taxa(rare,  Rank6 == "g__Cylindrospermopsis")
sum(rowSums(otu_table(only_one_genus))) / sum(rowSums(otu_table(only_cyano)))

sum(optical64[,grepl("Cylindrospermopsis", names(optical64))]) / sum(rowSums(optical64[2:ncol(optical)]))
range(optical64[,grepl("Cylindrospermopsis", names(optical64))])

#5 Chroococcus
only_one_genus <- subset_taxa(rare,  Rank6 == "g__Chroococcus")
sum(rowSums(otu_table(only_one_genus))) / sum(rowSums(otu_table(only_cyano)))

sum(optical64[,grepl("Chroococcus", names(optical64))]) / sum(rowSums(optical64[2:ncol(optical)]))
range(optical64[,grepl("Chroococcus", names(optical64))])

############################################################
##########  Supplementary Figure 1  compare two method   ##########  
############################################################


# # all toxic genus
# rel_opti <- optical64
# #for (i in 1:nrow(optical64)){
# #	rel_opti[i,2:ncol(rel_opti)] = optical64[i,2:ncol(rel_opti)] / sum( optical64[i,2:ncol(rel_opti)] )
# #}
# #rowSums(rel_opti[, 2:ncol(rel_opti)])
# toxic_optical <- rel_opti[,grepl("Microcystis|Planktothrix|Anabaena|Cylindrospermopsis|Chroococcus|Pseudanabaena|Pseudanabaenaceae|Synechococcaceae", names(rel_opti ))]
# sco <- data.frame(rel_opti$sampleid, rowSums(toxic_optical))
# colnames(sco) <- c("sampleid", "microcystis_scopy")

# #sequncing
# #sub <- subset_taxa(rare, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena"| Rank6 == "g__Dolichospermum")

# # sequencing data, relative abundance
# sub <- subset_taxa(GPr, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena"| Rank6 == "g__Dolichospermum")

# seq16s <- data.frame(sample_names(sub),sample_sums(sub))
# colnames(seq16s) <- c("sampleid", "microcystis_16s")

# merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
# cor(merged$microcystis_16s, merged$microcystis_scopy)
# pdf("~/Box Sync/2017/7July/water/figures_tables/all_toxic.pdf", width=6, height=6)
# ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Cyanobacteria 16S gene relative abundance", y="Toxi genus biomass  per liter of water(mg/L)")
# dev.off()
# cor.test(merged$microcystis_16s, merged$microcystis_scopy)
# subset(merged, microcystis_16s > 0.35 & microcystis_scopy <50)



# 1. Microcystis
micro <-  subset_taxa(GPr , Rank6 == "g__Microcystis")
seq16s <- data.frame(sample_names(micro),sample_sums(micro))
colnames(seq16s) <- c("sampleid", "microcystis_16s")
taxa_sums(micro)

sco <- data.frame(optical64$sampleid, optical64$Cyanophyta_Microcystis)
colnames(sco) <- c("sampleid", "microcystis_scopy")

merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
cor(merged$microcystis_16s, merged$microcystis_scopy)
#pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure1.pdf", width=6, height=6)
p1 <- ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Microcystis 16S rRNA gene ", y="Microcystis biomass per liter of water (mg/L)")+theme_bw()+ggtitle("A")
#dev.off()
cor.test(merged$microcystis_16s, merged$microcystis_scopy)

## under estimated risk assesment 
#subset(merged, microcystis_16s > 250 & microcystis_scopy < 50)
#     sampleid microcystis_16s microcystis_scopy
#32 2014202002             280          11.42902
#54 2014216002             412          10.70978
#63 2014224003             655          19.79966


# 2. Planktothrix
micro <-  subset_taxa(GPr , Rank6 == "g__Planktothrix")
seq16s <- data.frame(sample_names(micro),sample_sums(micro))
colnames(seq16s) <- c("sampleid", "microcystis_16s")

sco <- data.frame(optical64$sampleid, optical64$Cyanophyta_Planktothrix)
colnames(sco) <- c("sampleid", "microcystis_scopy")

merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
cor.test(merged$microcystis_16s, merged$microcystis_scopy)

range(otu_table(micro)[c("OTU-8058"),])


merged <- subset(merged,microcystis_16s < 0.5)
#pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure_Planktothrix.pdf", width=6, height=6)
p2 <- ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Planktothrix 16S rRNA gene", y="Planktothrix biomass per liter of water (mg/L)")+theme_bw()+ggtitle("B")
#dev.off()


# # how many in N limiting condition?
# micro <-  subset_taxa(rare, Rank6 == "g__Planktothrix")
# nlim <- prune_samples(sample_data(micro)$two_group == "Nlim", micro)
# sum(sample_sums(nlim)) /sum(sample_sums(micro))
# #dim(merged)


# 3. Dolichospermum (Anabaena): 
micro <-  subset_taxa(GPr ,  Rank6 == "g__Dolichospermum" | Rank6 == "g__Anabaena" )
seq16s <- data.frame(sample_names(micro),sample_sums(micro))
colnames(seq16s) <- c("sampleid", "microcystis_16s")

sco <- data.frame(optical64$sampleid, optical64$Cyanophyta_Anabaena)
colnames(sco) <- c("sampleid", "microcystis_scopy")

merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
cor(merged$microcystis_16s, merged$microcystis_scopy)
cor.test(merged$microcystis_16s, merged$microcystis_scopy)
#pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure2.pdf", width=6, height=6)
p3 <- ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Dolichospermum 16S gene", y="Dolichospermum biomass per liter of water (mg/L)")+theme_bw()+ggtitle("C")
#dev.off()
subset(merged, microcystis_16s > 0.4)
sample_data(GPr)["2014188003",]

micro <-  subset_taxa(GPr,  Rank6 == "g__Dolichospermum" | Rank6 == "g__Anabaena" )
sample_sums(micro)
sum(subset(optical64, sampleid == "2014188003")[2:ncol(optical64)])
#40.9/90

#4 Cylindrospermopsis
micro <-  subset_taxa(GPr,  Rank6 == "g__Cylindrospermopsis" )
seq16s <- data.frame(sample_names(micro),sample_sums(micro))
colnames(seq16s) <- c("sampleid", "microcystis_16s")

sco <- data.frame(optical64$sampleid, optical64$Cyanophyta_Cylindrospermopsis)
colnames(sco) <- c("sampleid", "microcystis_scopy")
merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
cor.test(merged$microcystis_16s, merged$microcystis_scopy)
p4 <- ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Cylindrospermopsis 16S gene", y="Cylindrospermopsis biomass per liter of water  (mg/L)")+theme_bw()+ggtitle("D")

#5 Chroococcus
micro <-  subset_taxa(GPr,  Rank6 == "g__Chroococcus" )
seq16s <- data.frame(sample_names(micro),sample_sums(micro))
colnames(seq16s) <- c("sampleid", "microcystis_16s")

sco <- data.frame(optical64$sampleid, optical64$Cyanophyta_Chroococcus)
colnames(sco) <- c("sampleid", "microcystis_scopy")
merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
cor.test(merged$microcystis_16s, merged$microcystis_scopy)
p5 <- ggplot(merged, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Chroococcus 16S gene", y="Chroococcus biomass per liter of water (mg/L)")+theme_bw()+ggtitle("E")


pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure1.pdf", width=13, height=8)
grid.arrange(p1,p2,p3,p4,p5, ncol = 3)
dev.off()




#zoom in
# sub <- subset(merged, microcystis_16s < 250)
# #dim(sub)
# pdf("~/Box Sync/2017/5May/water/Figure2b.pdf", width=6, height=6)
# ggplot(sub, aes(microcystis_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()+labs(x="Microcystin 16S gene abundance", y="Microcystin biomass (mg/L)")
# dev.off()
# cor.test(sub$microcystis_16s, sub$microcystis_scopy)
#4000, 0.5
#20000, 0.722



#temp <- t( as.data.frame( otu_table(micro)) )
#merged <- merge(temp, sco, by.x = "row.names", by.y = "sampleid")


# # cor(merged[,-c(1)])
# head(merged[,-c(1)])



# ggplot(merged, aes(New.CleanUp.ReferenceOTU99005, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()
# ggplot(merged, aes(New.CleanUp.ReferenceOTU94718, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()


# test<-unique(genus_16s_nog[genus_16s_nog %in% micro_name2])
# fin_table = data.frame()
# for (x in test){
	# sp = paste0("g__", x)
	# phy_sp = paste0("Cyanophyta_",x)
	# if( !(sp %in% c("g__Pseudanabaenaceae","g__Synechococcaceae")  )   ){
		# temp <- subset_taxa(physeq, Rank6 == sp)
		# seq16s <- data.frame(sample_names(micro),sample_sums(micro))
		# colnames(seq16s) <- c("sampleid", "seq_16s")
		
		# phyto[,phy_sp]
		# sco <- data.frame(phyto$sampleid, phyto[,phy_sp])
		# colnames(sco) <- c("sampleid", "microcystis_scopy")
		
		# merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
		# p1 = ggplot(merged, aes(seq_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()
		# temp = data.frame(cor(merged$microcystis_scopy, merged$seq_16s), x)
		# colnames(temp) = c("cor","genus")
		# fin_table = rbind(fin_table, temp)
		# filename = paste0(x,".pdf")
		# pdf(filename, width=6, height=6)
		# print(p1)
		# dev.off()
	# }
		
# }



# x = "Planktothrix"
# sp = paste0("g__", x)
# temp <- subset_taxa(physeq, Rank6 == sp)
# seq16s <- data.frame(sample_names(micro),sample_sums(micro))
# colnames(seq16s) <- c("sampleid", "seq_16s")
# sco <- data.frame(phyto$sampleid, phyto$Cyanophyta_Planktothrix)
# colnames(sco) <- c("sampleid", "microcystis_scopy")
# merged <- merge(seq16s, sco, by.x = "sampleid", by.y = "sampleid")
# ggplot(merged, aes(seq_16s, microcystis_scopy))+ geom_smooth(method="lm")+geom_point()
# cor(merged[,-c(1)])
# cor.test(merged$seq_16s, merged$microcystis_scopy)

# x = "Anabaena"
# ge = paste0("Cyanophyta_", x)
# sco <- data.frame(phyto$sampleid, phyto$Cyanophyta_Anabaena)


#//////////////////////////////////////////////////////////////////////////////////////////////////////////#
##########    Figure 3 Microcystis 33 OTU heatmap   ##########  
#//////////////////////////////////////////////////////////////////////////////////////////////////////////#
##########  Figure 3  : heatmap, add additional line ##########  
library(scales)

micro <-  subset_taxa(rare , Rank6 == "g__Microcystis")
log_micro <- micro
otu_table(log_micro) <- otu_table(micro)+1
sample_data(micro)$name_group <- paste(rownames(sample_data(micro)), sample_data(micro)$two_group, sep="_")

png("~/Box Sync/2017/7July/water/figures_tables/Figure3_heatmap.png", width=12, height= 6,unit="in", res=300)
plot_heatmap(micro, "RDA", "none", sample.label= "lake.name", low="#66CCFF", high="#000033", na.value="white")
#plot_heatmap(micro, "PCoA", "unifrac", "name_group") 
dev.off()


micro_otu <- as.data.frame(tax_table(micro))
rownames(micro_otu)

otu_info <- read.csv("~/Box Sync/2017/7July/water/internal_info.csv")
head(otu_info)

micro_name <- as.character(subset(otu_info, OTU.name %in% rownames(micro_otu))$Qiime.assigned.OTU.name)
write(micro_name,"~/Box Sync/2017/7July/water/micro_name.txt")


# plot_heatmap(micro, "PCoA", "bray")
# plot_heatmap(micro, "PCoA", "unifrac", "new_group") #good
# plot_heatmap(micro, "MDS", "unifrac", "new_group") #good
# plot_heatmap(micro, "NMDS", "bray") #change everytime
# plot_heatmap(micro, "NMDS", "jaccard") #change everytime
# 


library(gplots)
heatmap.2(otu_table(micro), density.info="none", trace="none", dendrogram=c("row"), 
            symm=F,symkey=T,symbreaks=T, scale="none")

#ratio microcystis
range(optical[,c("Cyanophyta_Microcystis")] /  rowSums(optical[,2:ncol(optical)]))


# 2. Planktothrix
micro <-  subset_taxa(GPr , Rank6 == "g__Planktothrix")
gpt <- prune_taxa(names(sort(taxa_sums(micro),TRUE)[1:20]), micro)

sample_data(micro)$name_group <- paste(rownames(sample_data(micro)), sample_data(micro)$two_group, sep="_")
plot_heatmap(micro, "PCoA", "unifrac", "name_group") 

names(sort(taxa_sums(micro), TRUE)[1:10])
names(sort(sample_sums(micro), TRUE)[1:10])
# optical$sampleid
# optical$Cyanophyta_Microcystis
# opt <- data.frame(optical$sampleid,optical$Cyanophyta_Microcystis)
# opt <- data.frame(optical$Cyanophyta_Microcystis)
# row.names(opt) <- optical$sampleid
# opt2 <- cbind(opt,opt)

# row.names(opt) <- opt$optical.sampleid
# optical.Cyanophyta_Microcystis
# opt2 <- data.frame(opt$optical.Cyanophyta_Microcystis, opt$optical.Cyanophyta_Microcystis)
# opt2 <- as.matrix(opt2)
# row.names(opt2) <- opt$optical.sampleid

# heatmap(as.matrix(opt2))
# heatmap(as.matrix(optical[,-c(1)] ) )
# as.matrix(optical, row.names=1)
# otu <- t(otu_table(micro))
# heatmap(otu)

# merged <- merge(otu, opt, by.x="row.names", by.y="optical.sampleid")
# tmer <- t(merged[,-c(1)])
# heatmap(tmer)
# row.names(merged) <- merged$Row.names
# merged_tr <- t(as.matrix(merged[,-c(1)]))
# heatmap(merged_tr)






#,trans=log_trans(10)
#otu_table(gpt)[which(otu_table(gpt)==0)]=1
# plot_bar(micro, fill = "Rank7")
# plot_bar(micro, fill = row.names(tax_table(micro)))
# row.names(tax_table(micro))

# pdf("figure5_headmap.pdf", width=25, height= 12)
# #p1 <- plot_bar(micro, fill = "OTU")
# p1 <- plot_tree(micro, color="OTU", size="Abundance")
# print(p1)
# dev.off()
# plot_heatmap(micro, "PCoA", "bray")
# plot_heatmap(micro, "PCoA", "unifrac", "new_group") #good
# plot_heatmap(micro, "MDS", "unifrac", "new_group") #good
# plot_heatmap(micro, "NMDS", "bray") #change everytime
# plot_heatmap(micro, "NMDS", "jaccard") #change everytime
# plot_heatmap(micro, "RDA", "none")

# heatmap(otu_table(micro))
# pdf("figure5_heatmap.pdf", width=6, height= 6)
# #png("figure5_heatmap.png", width=6, height= 6, unit="in", res=300)
# plot_heatmap(micro, "PCoA", "unifrac", "new_group") #good
# dev.off()


# p1$data
# library("data.table")
# newtab = data.table(p1$data)
# setorder(newtab, Abundance)
# p1$data <- newtab
# print(p1)


# test <- readRDS("~/Box Sync/2017/1January/water/water64/water-project-object_physeq_before_normalization.rds")
# test.df <- as.data.frame(otu_table(test))
# rarecurve(test.df)




#//////////////////////////////////////////////////////////////////////////////////////////////////////////#
##########    Figure 4 toxic and non-toxic    ##########  
#//////////////////////////////////////////////////////////////////////////////////////////////////////////#\
library(ape)

########## three OTU ##########
#tree = read.tree("~/Box Sync/2017/4April/water1/test.dnd")
#tree = read.tree("~/Box Sync/2017/4April/water1/test.clustal.tree")
tree = read.tree("~/Box Sync/2017/4April/water1/test.fast.tree")
tree$tip.label <- c("OTU-7860","OTU-7831","OTU-7829","Microcystis aeruginosa PCC7806","OTU-7862","Microcystis wesenbergii NIES112","Microcystis wesenbergii NIES-107","Microcystis aeruginosa PCC7005","Microcystis aeruginosa NIES-98","OTU-7863","Microcystis elabens NIES42","Microcystis holsatica NIES43","OTU-7859","Microcystis sp. AWT139","Microcystis aeruginosa PCC7820","Microcystis aeruginosa strain NIES-843","Microcystis aeruginosa NIES89")
#collab <- c("#EB6841","#EB6841","#EB6841","#EB6841","#EB6841","#000000","green","black","blue","black","green","blue","#EB6841","#EB6841","green","blue","blue")

collab <- c("green","black","green","#EB6841","black","blue","#EB6841","#EB6841","blue","black","blue","blue","green","#EB6841","#EB6841","#EB6841","#EB6841")

#tree_rooted <- root(tree, outgroup="U40336.1",resolve.root=T)
i#s.rooted(tree_rooted)
pdf("~/Box Sync/2017/7July/water/figures_tables/Figure4.pdf",width=8,height=8)
plot(tree, type='phylogram', tip.color=collab)
legend("bottomright", inset=.05, c("Toxic Microcystis","Non-toxic Microcystis","Three most abundant OTUs"), fill=c("#EB6841","blue", "green"), horiz=F)
dev.off()

# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
# library(ggtree)
# plot_tree(tree)

# ggtree(tree)
# library(ggplot2)
# library(Biostrings)
# nwk <- system.file("extdata","~/Box Sync/2017/4April/water1/test.dnd", package="ggtree")
# tr <- read.tree(nwk)


########## all OTU ##########
tree = read.tree("~/Box Sync/2017/4April/water1/test_all.dnd")
collab <- ifelse(tree$tip.label =="U40339.1" | tree$tip.label == "U40331.2" | tree$tip.label == "U40333.2"|tree$tip.label == "U03403.1"|tree$tip.label == "U03402.1"|tree$tip.label == "U40338.1"|tree$tip.label == "NR_074314.1", "red" , ifelse(tree$tip.label =="U40337.2" |tree$tip.label =="U40334.1"|tree$tip.label =="U40335.1"|tree$tip.label =="U40336.1","blue",
 ifelse(tree$tip.label =="3237" |tree$tip.label =="4351168"|tree$tip.label =="3704350","green","black" )
 ))

pdf("test_all.pdf",width=8,height=23)
plot(tree, type='c', tip.color=collab)
dev.off()

##########  Figure toxic genus ##########  
toxic_genus = c("Microcystis","Planktothrix","Anabaena")
tax <- as.data.frame(tax_table(physeq))
subset(tax, Rank6 == "g__Nostoc")
colnames(tax_table(physeq))


sub <- subset_taxa(physeq, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena")
tax_table(sub)


## 6 toxic genus in seq
toxic <- subset_taxa(physeq, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis"| Rank6 == "g__Nostoc")

## 64 sample of optical
opti64 <- subset(optical, sampleid %in% sample_names(physeq))
dim(opti64)

	dim(optical)
head(optical)

physeq
sample_names(physeq)

########################################
##########  Figure 5 all phylum ##########  
########################################

# pdf("~/Box Sync/2017/7July/water/figures_tables/group_all.pdf", width=12, height= 6)
# plot_bar(GPf, fill = "Rank2")+facet_grid(new_group~.)
# dev.off()
#plot_bar(one_phylum, fill = "Rank3")+facet_grid(new_group~.)

for_plot <- GPf
tax_table(for_plot)[is.na(tax_table(for_plot))] <-"Unassigned"
tax_table(for_plot) <- gsub("f__\\[Chromatiaceae]","Chromatiaceae",tax_table(for_plot))
tax_table(for_plot) <- gsub("p__","",tax_table(for_plot))
tax_table(for_plot) <- gsub("\\bf__\\b","Unassigned",tax_table(for_plot))
tax_table(for_plot) <- gsub("f__","",tax_table(for_plot))

physeqdf <- psmelt(for_plot)
head(physeqdf)
unique(physeqdf$Rank6)

#ggplot(physeqdf, aes(x = reorder(physeqdf$Sample, -ifelse(is.na(physeqdf$Rank2), 0, ifelse(physeqdf$Rank2=="Cyanobacteria" , physeqdf$Abundance, 0))), y= Abundance, fill=Rank2))+ geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Sample", y="Relative Abundance")+scale_fill_discrete(name="Phylum")

pdf("~/Box Sync/2017/7July/water/figures_tables/Figure5-bar_chart_all.pdf", width=12, height= 10)
ggplot(physeqdf, aes(x = reorder(physeqdf$lake.name, -ifelse(is.na(physeqdf$Rank2), 0, ifelse(physeqdf$Rank2=="Cyanobacteria" , physeqdf$Abundance, 0))), y= Abundance, fill=Rank2))+ geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Sample", y="Relative Abundance")+scale_fill_discrete(name="Phylum")
dev.off()


####################################################
##########  Supplementary Figure 3 compare phylum ##########  
####################################################
one_phylum <- subset_taxa(for_plot, Rank2 == "Actinobacteria")
physeqdf_actino <- psmelt(one_phylum)
#pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure4-bar_chart_actinobacteria.pdf", width=12, height= 6)
p1 <- ggplot(physeqdf_actino, aes(x = reorder(physeqdf_actino$lake.name, -ifelse(is.na(physeqdf_actino$Rank2), 0, ifelse(physeqdf_actino$Rank2=="Actinobacteria" , physeqdf_actino$Abundance, 0))), y= Abundance, fill=Rank5))+ geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Sample", y="Relative Abundance")+ggtitle("A. Actinobacteria")
#dev.off()

one_phylum <- subset_taxa(for_plot, Rank2 == "Proteobacteria")
physeqdf_pro <- psmelt(one_phylum)
#pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure5-bar_chart_proteobacteri.pdf", width=12, height= 6)
p2 <- ggplot(physeqdf_pro, aes(x = reorder(physeqdf_pro$lake.name, -ifelse(is.na(physeqdf_pro$Rank2), 0, ifelse(physeqdf_pro$Rank2=="Proteobacteria" , physeqdf_pro$Abundance, 0))), y= Abundance, fill=Rank5))+ geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="Sample", y="Relative Abundance")+ggtitle("B. Proteobacteria")
#dev.off()

pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure3.pdf", width=13, height=14)
grid.arrange(p1,p2, ncol = 1)
dev.off()

#find most popular family among Proteobacteria
one_phylum <- subset_taxa(GPf, Rank2 == "p__Proteobacteria")
glom <- tax_glom(one_phylum, "Rank5")
sort(taxa_sums(glom))


major_phylum <- subset_taxa(rare, Rank5 == "f__Pelagibacteraceae")
sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
major_phylum <- subset_taxa(rare, Rank5 == "f__Comamonadaceae")
sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
major_phylum <- subset_taxa(rare, Rank5 == "f__Rhodocyclaceae")
sum(sample_sums(major_phylum)) / sum(sample_sums(rare))
major_phylum <- subset_taxa(rare, Rank5 == "f__Methylophilaceae")
sum(sample_sums(major_phylum)) / sum(sample_sums(rare))

one_genus <- subset_taxa(rare, Rank5 == "f__Pelagibacteraceae")
one_genus <- subset_taxa(rare, Rank5 == "g__SAR11")

one <- prune_samples(sample_names(GPf) == "2014202011",GPf)
one <- prune_samples(sample_names(GPf) == "2014188003",GPf)
one <- prune_samples(sample_names(GPf) == "2014216006",GPf)
one2 <- subset_taxa(one, Rank2 == "p__Cyanobacteria")
plot_bar(one2, fill="Rank6")

one <- prune_samples(sample_names(GPf) == "2014202012",GPf)
one <- prune_samples(sample_names(GPf) == "2014203008",GPf)
one <- prune_samples(sample_names(GPf) == "2014202012",GPf)
one2 <- subset_taxa(one, Rank2 == "p__Actinobacteria")
one2 <- subset_taxa(one, Rank2 == "p__Proteobacteria")
plot_bar(one2, fill="Rank6")


one3 <- subset_taxa(physeq, Rank6 == "g__Anabaena")
one3 <- subset_taxa(physeq, Rank6 == "g__Dolichospermum")
one3 <- subset_taxa(physeq, Rank6 == "g__Planktothrix")

tax_table(one3)
sample_sums(one3)

sample_data(one3)

# correlation planktothrix

#use relative abundance
GPr
# microscopy data
phyto <- read.csv("~/Box Sync/2016/9September/water/phytoplankton ds FINAL.csv")
#only cyanobacteria
onlycyno <- phyto[,grepl("Cyanophyta", names(phyto))]
optical <- cbind(sampleid = phyto$sampleid, onlycyno)


fin_table <- data.frame()
genus <- c("Planktothrix", "Dolichospermum","Microcystis","Cylindrospermopsis","Chroococcus")
genera <- "Planktothrix"
"Microcystis"

for (genera in genus){
seq_genus <- paste0("g__", genera)
one_genus <- subset_taxa(GPr, Rank6 == seq_genus)
sum <- data.frame(colnames(otu_table(one_genus)), colSums(otu_table(one_genus)))
colnames(sum) <- c("sampleid", "seq_value")
average<- ave(sum$seq_value)
mi <- min(sum$seq_value)
ma <- max(sum$seq_value)
 if (genera == "Dolichospermum"){
 	genera = "Anabaena"
 }
opti_genus <- paste0("Cyanophyta_", genera)
optiB <- optical[,c("sampleid", opti_genus)]
colnames(optiB) <- c("sampleid", "opti_value")

merged <- merge(sum, optiB, by.x = "sampleid", by.y = "sampleid")

p=ggplot(merged, aes(x=seq_value, y=opti_value))+geom_point()+geom_smooth()
temp_file <- paste0("cor_", genera)
filename <- paste0(temp_file, ".pdf")
pdf(filename)
plot(p)
dev.off()
co_value <- cor(merged$seq_value, merged$opti_value)
#0.9
temp <- data.frame(genera,co_value, average[1], mi,ma)
fin_table <- rbind(fin_table,temp)
}
heatmap(otu_table(one_genus))
plot_heatmap(one_genus, "PCoA", "unifrac", "name_group") 
merged[order(merged$seq_value),]

#plank
four_plank <- prune_samples(sample_names(GPr) == "2014202011"|sample_names(GPf) == "2014195009"|sample_names(GPr) == "2014210005"|sample_names(GPr) == "2014210002",GPr)
four_plank2 <- subset_taxa(four_plank, Rank6 == "g__Planktothrix")
plot_heatmap(four_plank2, "PCoA", "unifrac") 

heatmap(otu_table(four_plank2))

png("heatmap_plank.png")
plot_heatmap(four_plank2, "NMDS", "bray") 
dev.off()

#microcysits

one_genus <- subset_taxa(GPr, Rank6 == "g__Microcystis")
sum <- data.frame(colnames(otu_table(one_genus)), colSums(otu_table(one_genus)))
colnames(sum) <- c("sampleid", "seq_value")
most_high <- sum[rev(order(sum$seq_value)),]$sampleid[1:4]

four_plank <- prune_samples(sample_names(one_genus) %in% most_high,one_genus)
plot_heatmap(four_plank, "PCoA", "unifrac") 
png("heatmap_mc.png")
plot_heatmap(four_plank, "NMDS", "bray") 
dev.off()


################################################################################
##########  Supplementary Figure 2 Table  Stat shows Cyanobacteria vs Actinobacteria and proteobacteria ##########  
################################################################################

one_phylum <- subset_taxa(GPr, Rank2 == "p__Cyanobacteria")
cyano <- as.data.frame(sample_sums(one_phylum))
one_phylum <- subset_taxa(GPr, Rank2 == "p__Actinobacteria")
actino <- as.data.frame(sample_sums(one_phylum))
one_phylum <- subset_taxa(GPr, Rank2 == "p__Proteobacteria")
proteo <- as.data.frame(sample_sums(one_phylum))

fin_table <- data.frame(cyano,actino,proteo)
colnames(fin_table) <- c("Cyanobacteria","Actinobacteria","Proteobacteria")

p1 <- ggplot(fin_table, aes(Cyanobacteria, Actinobacteria))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Cyanobacteria 16S gene", y="Relative abundance of Actinobacteria 16S gene")+theme_bw()+ggtitle("A. Actinobacteria")
p2 <- ggplot(fin_table, aes(Cyanobacteria, Proteobacteria))+ geom_smooth(method="lm")+geom_point()+labs(x="Relative abundance of Cyanobacteria 16S gene", y="Relative abundance of Proteobacteria 16S gene")+theme_bw()+ggtitle("B. Proteobacteria")

pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure2.pdf", width=13, height=8)
grid.arrange(p1,p2, ncol = 2)
dev.off()

cor.test(fin_table$Cyanobacteria, fin_table$Actinobacteria)
cor.test(fin_table$Cyanobacteria, fin_table$Proteobacteria)
cor(fin_table)
 
# # one_phylum <- subset_taxa(rare, Rank2 == "p__Cyanobacteria")
# cyano <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Actinobacteria")
# actino <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Proteobacteria")
# proteo <- as.data.frame(sample_sums(one_phylum))

# one_phylum <- subset_taxa(rare, Rank2 == "p__Cyanobacteria")
# cyano <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Armatimonadetes")
# arm <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Bacteroidetes")
# bac <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Chlorobi")
# ch <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Planctomycetes")
# pl <- as.data.frame(sample_sums(one_phylum))
# one_phylum <- subset_taxa(rare, Rank2 == "p__Verrucomicrobia")
# ve <- as.data.frame(sample_sums(one_phylum))

# fin_table <- data.frame(cyano,actino,proteo, arm,bac,ch,pl,ve)
# cor.test(fin_table[,1], fin_table[,8])

# fin_table <- data.frame(cyano,actino,proteo)
# colnames(fin_table) <- c("Cyanobacteria","Actinobacteria","Proteobacteria")


# # library(corrplot)
# # M <- cor(fin_table)
# # corrplot(M, method="circle")
# # corrplot.mixed(M)
# install.packages("PerformanceAnalytics")
# library(PerformanceAnalytics)
# pdf("~/Box Sync/2017/7July/water/figures_tables/Supplementary_Figure3.pdf", width=6, height=6)
# chart.Correlation(fin_table)
# dev.off()
# #add other factors
# fin_table <- cbind(fin_table,sample_data(rare))

# # 
# cor(fin_table[,c(1,2,3,7,8,9,10,11,12,13,14,15,16,17,18)], method="pearson")
# cor(fin_table, method="spearman")


# fin_table <- data.frame()
# toxic <- subset_taxa(GPf, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena"| Rank6 == "g__Dolichospermum")
# sample_sums(toxic)
# one_phylum <- subset_taxa(GPr, Rank2 == "p__Actinobacteria")
# glom <- tax_glom(one_phylum, "Rank5")
# tax_table(glom)
# taxa_sums(glom)
# mean(otu_table(glom)[33,])
# range(otu_table(glom)[33,])
# mean(otu_table(glom)[26,])
# range(otu_table(glom)[26,])
# mean(otu_table(glom)[7,])
# range(otu_table(glom)[7,])
# totu <- t(otu_table(glom))
# fin_table <- cbind(toxic=sample_sums(toxic), totu)
# cor(fin_table)
# cor.test(fin_table[,1], fin_table[,34])
# cor.test(fin_table[,1], fin_table[,8])

# one_phylum <- subset_taxa(GPr, Rank2 == "p__Proteobacteria")
# glom <- tax_glom(one_phylum, "Rank5")
# tax_table(glom)
# totu <- t(otu_table(glom))
# fin_table <- cbind(toxic=sample_sums(toxic), totu)
# co <- cor(fin_table)
# range(co[-1,1])
########################################################################
########################################################################
########################################################################

##########  Figure 0: ven diagram, compare microscopy vs. 16s, how many genus catch? ############
optical_tax <- colnames(optical)[-c(1)]
write(optical_tax, "optical_tax.txt")
micro_name <- names(optical[,-c(1)])
micro_name1<-gsub(".*_","",micro_name)
micro_name2 <- gsub("\\..*","",micro_name1)

tax <- as.data.frame( tax_table(physeq))
class(tax)
dim(tax)
head(tax)

uniq_tax <- tax[!duplicated(tax),]
dim(uniq_tax)
head(uniq_tax)

genus_tax <- uniq_tax[,-c(7)]
head(genus_tax )
uniq_genus_tax <- genus_tax[!duplicated(genus_tax),]
dim(uniq_genus_tax)
head(uniq_genus_tax)
class(uniq_genus_tax)
dat <- as.matrix(uniq_genus_tax)
dat[is.na(dat)] <- "not"
uniq_genus_tax <- as.data.frame(dat)
uniq_genus_tax <- dat

uniq_tax_16s <- as.matrix(uniq_genus_tax$Rank6[!duplicated(uniq_genus_tax$Rank6)])
tax_16s <- as.matrix(as.data.frame(tax_table(physeq))$Rank6)
write(tax_16s, "16s_tax.txt")
length(uniq_tax_16s)
for (i in 1:nrow(uniq_genus_tax)){
	if(!is.na(uniq_genus_tax[i,5])){
		if(uniq_genus_tax[i,5] == "f__Pseudanabaenaceae"){
			uniq_genus_tax[i,6] <- "g__Pseudanabaenaceae"
		}
		if(uniq_genus_tax[i,5] == "f__Synechococcaceae"){
			uniq_genus_tax[i,6] <- "g__Synechococcaceae"
		}		
	}
}
uniq_genus_tax <- as.data.frame(uniq_genus_tax)
genus_16s <- uniq_genus_tax$Rank6
genus_16s_nog <- gsub("g__","",genus_16s)

genus_16s_nog[genus_16s_nog %in% micro_name2]
#"Pseudanabaenaceae" %in% genus_16s_nog

#one working
require(gplots)
venn(list(sequencing_16s=genus_16s_nog, Microscopy=micro_name2))

#both
sub <- subset_taxa(physeq, Rank6 == "g__Anabaena"| Rank6 == "g__Chroococcus" | Rank6 == "g__Cylindrospermopsis" | Rank6 == "g__Planktothrix"  | Rank6 == "g__Microcystis" | Rank6 == "g__Pseudanabaena")

glom <- tax_glom(sub, "Rank6")
taxa_sums(glom)
tax_table(glom)




##########  Figure ?  ##########  
physeq <- readRDS("~/Box Sync/2017/1January/water/water64/water-project-object_physeq_normalized.rds")

p = plot_bar(physeq, fill="Rank2")
head(p)
names(p)
#"data"        "layers"      "scales"      "mapping"     "theme"       "coordinates" "facet"       "plot_env"    "labels"  
p$data[1:10,1:10]
dim(p$data)
rnew_table <- data.frame()




melted <- melt(otu_table(physeq))
colnames(melted) <- c("OTU","Sample","value")
merged <- merge(melted, )
head(melted)

	new_table <- data.frame(taxa_sums(physeq), tax_table(physeq)[,2])
	colnames(new_table) = c("abundance", "phylum")
	sum_table = data.frame()
	for (x in unique(new_table$phylum) ){
		temp = subset(new_table, phylum==x)
		su = sum(temp$abundance)
		sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
	}
	colnames(sum_table) = c("abundance", "phylum")
	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
	
	other_table = data.frame()
	other = c()
	for (i in 1:nrow(perc_table)){
		if(perc_table[i,1] > 0.01) {
			other_table = rbind(other_table, perc_table[i,])
		}else{
			other = c(other, perc_table[i,1])
		}
		
	}
	

ggplot(sum_table, aes(x=group,y=abundance, fill=phylum))+geom_bar(stat="identity")+labs(y="relative abundance")

####################################

#two group
sample_data(physeq)$two_group <- sample_data(physeq)$new_group
sample_data(physeq)$two_group <- gsub("PlimHC", "Plim", sample_data(physeq)$two_group)
sample_data(physeq)$two_group <- gsub("PlimLC", "Plim", sample_data(physeq)$two_group)
sample_data(physeq)$two_group <- gsub("NlimHC", "Nlim", sample_data(physeq)$two_group)
sample_data(physeq)$two_group <- gsub("NlimLC", "Nlim", sample_data(physeq)$two_group)

fin_table  = data.frame()
for (group in unique(sample_data(physeq)$two_group)){
	temp_physeq = prune_samples(sample_data(physeq)$two_group == group, physeq)

	new_table <- data.frame(taxa_sums(temp_physeq), tax_table(temp_physeq)[,2])
	colnames(new_table) = c("abundance", "phylum")
	sum_table = data.frame()
	for (x in unique(new_table$phylum) ){
		temp = subset(new_table, phylum==x)
		su = sum(temp$abundance)
		sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
	}
	colnames(sum_table) = c("abundance", "phylum")
	perc_table = sum_table
	for (i in 1:nrow(sum_table)){
			perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
	}
	
	other_table = data.frame()
	other = c()
	for (i in 1:nrow(perc_table)){
		if(perc_table[i,1] > 0.01) {
			other_table = rbind(other_table, perc_table[i,])
		}else{
			other = c(other, perc_table[i,1])
		}
		
	}
	
	sum(other)
	tep = data.frame(sum(other), "other")
	colnames(tep) = c("abundance", "phylum")
	tfin = rbind(other_table, tep)
	ttfin = cbind(tfin,group)
	fin_table = rbind(fin_table, ttfin)
	phy = unique(fin_table$phylum)
}

fin_table$group = factor(fin_table$group, levels = c("Nlim","Plim"))
pdf("Figure4.pdf", width=6, height=6)
ggplot(fin_table, aes(x=group,y=abundance, fill=phylum))+geom_bar(stat="identity")+labs(y="relative abundance")
dev.off()

##########  Figure rarefaction curve ##########  
physeq <- readRDS("~/Box Sync/2017/1January/water/water64/water-project-object_physeq_before_normalization.rds")

#two group
sample_data(physeq)$two_group <- sample_data(physeq)$new_group
sample_data(physeq)$two_group <- gsub("PlimHC", "Plim", sample_data(physeq)$two_group)
sample_data(physeq)$two_group <- gsub("PlimLC", "Plim", sample_data(physeq)$two_group)
sample_data(physeq)$two_group <- gsub("NlimHC", "Nlim", sample_data(physeq)$two_group)
sample_data(physeq)$two_group <- gsub("NlimLC", "Nlim", sample_data(physeq)$two_group)


for (group in unique(sample_data(physeq)$two_group)){
	temp_physeq = prune_samples(sample_data(physeq)$two_group == group, physeq)
	raredata <- t(otu_table(temp_physeq))
	
	filename = paste0(group,'.pdf')
	pdf(filename, width=6, height=6)
	rarecurve(raredata)
	dev.off()

}



# raredata <- t(otu_table(diver))
# rare2 <- raredata[1:2,]

# raremax <- min(rowSums(rare2))
# p=rarecurve(rare2, sample = raremax)



###
not_in_db <- c("Anabaenopsis","Aphanizomenon","Aphanothece","Coelosphaerium","Cuspidothrix.Aphanizomenon.","Cyanodictyon","Gomphosphaeria","Heteroleibleinia","Jaaginema","Limnothrix","Planktolyngbya","Pseudanabaena","Pseudanabaenaceae.family.","Raphidiopsis","Rhabdogloea","Romeria","Snowella","Synechococcaceae.family.","Synechocystis","Woronichinia")
names(optical64) %in% not_in_db
all <- sum(colSums(optical64[,2:ncol(optical64)]))
sort(colSums(optical64[,2:ncol(optical64)])/all )

