# DrosEU-metagentics
The bioinformatics pipeline for https://pubmed.ncbi.nlm.nih.gov/32003146/.
Wang Y, Kapun M, Waidele L, Kuenzel S, Bergland AO, Staubach F. Common structuring principles of the Drosophila melanogaster microbiome on a continental scale and between host and substrate. Environ Microbiol Rep. 2020 Apr;12(2):220-228. doi: 10.1111/1758-2229.12826. Epub 2020 Feb 8. PMID: 32003146.

raw data: https://www.ncbi.nlm.nih.gov/bioproject/?term=prjna515407

## A) mothur script

```mothur
####################mothur script######################
#######quality assessment and preparation of sequences
#pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319, keepdots=F, processors=8)
make.contigs(file=stability.rename.files, processors=8);
summary.seqs(fasta=current);
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275,processors=8);
unique.seqs(fasta=stability.files.trim.contigs.good.fasta);
count.seqs(name=stability.files.trim.contigs.good.names, group=stability.files.contigs.good.groups);
align.seqs(fasta=current,reference=silva.nr_v132.pcr.align);
summary.seqs(fasta=stability.files.trim.contigs.good.unique.align, count=current);
screen.seqs(fasta=stability.files.trim.contigs.good.unique.align, count=current, summary=stability.files.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8);
summary.seqs(fasta=current, count=current);
filter.seqs(fasta=current,vertical=T,trump=.);
unique.seqs(fasta=current,count=current);
pre.cluster(fasta=current,count=current,diffs=2);
chimera.uchime(fasta=current,count=current,dereplicate=t);
remove.seqs(fasta=current,accnos=current);
get.current();

# for stack plot Figure S1 and 
# and calculate wolbachia percentage in Table 
classify.seqs(fasta=current, count=current, reference=silva.nr_v132.pcr.align, taxonomy=silva.nr_v132.tax, cutoff=80,processors=8);
remove.lineage(fasta=stability.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.files.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=stability.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
classify.seqs(fasta=current, count=current, reference=silva.nr_v132.pcr.align, taxonomy=silva.nr_v132.tax, cutoff=80,processors=8)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Wolbachia);
classify.seqs(fasta=current, count=current, reference=silva.nr_v132.pcr.align, taxonomy=silva.nr_v132.tax, cutoff=80);

######assign sequences to OTUs
rename.file(fasta=current, count=current, taxonomy=current, prefix=start);
cluster.split(fasta=start.fasta, count=start.count_table, taxonomy=start.taxonomy, splitmethod=classify, taxlevel=4,method=average, cutoff=0.15,processors=16)
make.shared(list=start.an.unique_list.list, count=start.count_table, label=unique)
rename.file(input=start.an.unique_list.unique.pick.shared,new=start.all.shared)
count.groups(shared=start.all.shared)
#get taxonomy for each otu
classify.otu(list=start.an.unique_list.list, count=start.count_table, taxonomy=start.taxonomy, label=unique)
#get representative sequence for each OTU
get.oturep(list=start.an.unique_list.list, count=start.count_table,fasta=start.fasta,method=abundance,label=unique)

###get DrosEu samples
get.groups(shared=start.all.shared,groups=AT_Mau_14_01-UK_Sou_14_10-CY_Nic_14_11-UK_Mar_14_12-UK_Lut_14_13-DE_Bro_14_14-DE_Bro_14_15-UA_Yal_14_16-UA_Yal_14_17-UA_Yal_14_18-UA_Ode_14_19-AT_Mau_14_02-UA_Ode_14_20-UA_Ode_14_21-UA_Ode_14_22-UA_Kyi_14_23-UA_Kyi_14_24-UA_Var_14_25-UA_Pyr_14_26-UA_Dro_14_27-UA_Cho_14_28-UA_Cho_14_29-TR_Yes_14_03-SE_Lun_14_30-DE_Mun_14_31-DE_Mun_14_32-PT_Rec_14_33-ES_Gim_14_34-ES_Gim_14_35-FI_Aka_14_36-FI_Aka_14_37-FI_Ves_14_38-DK_Kar_14_39-TR_Yes_14_04-DK_Kar_14_40-DK_Kar_14_41-CH_Cha_14_42-CH_Cha_14_43-AT_See_14_44-UA_Kha_14_45-UA_Kha_14_46-UA_Cho_14_47-UA_Cho_14_48-UA_Kyi_14_49-FR_Vil_14_05-UA_Uma_14_50-FR_Vil_14_06-FR_Vil_14_07-FR_Got_14_08-UK_She_14_09)
rename.file(input=start.all.unique.pick.shared,new=start.droseu.shared)
count.groups(shared=start.droseu.shared)
#subsample 5135 sequences
sub.sample(shared=start.droseu.shared,size=5135,persample=true)

#redo tax.summary file for analysis at family level (droseu)
sub.sample(fasta=start.fasta, count=start.count_table,size=5135,persample=true)
classify.seqs(fasta=current, count=current, reference=silva.nr_v132.pcr.align, taxonomy=silva.nr_v132.tax, cutoff=80);
rename.file(fasta=start.subsample.fasta, count=start.subsample.count_table,prefix=5203)
rename.file(input=start.subsample.nr_v132.wang.tax.summary,new=5203start.subsample.nr_v132.wang.tax.summary)

###get fly versus substrate samples
get.groups(shared=start.all.shared,groups=grape_A2_f-grape_A2_s-grape_A5_f-grape_A5_s-apple_A7_f-apple_A7_s-apple_A8_f-apple_A8_s-apple_B7_f-apple_B7_s-cherry_C7_f-cherry_C7_s-plum_D8_f-plum_D8_s-apple_a3_f-apple_a3_s-apple_a4_f-apple_a4_s-apple_Bars_s-apple_Bars_f-cactus_s-cactus_f-lemon_s-lemon_f)
rename.file(input=start.all.unique.pick.shared,new=start.flysubs.shared)
count.groups(shared=start.flysubs.shared)
#sub sample 893 sequences
sub.sample(shared=start.flysubs.shared,size=893,persample=true)

#redo tax.summary file for analysis at family level (substrate vs fly)
sub.sample(fasta=start.fasta, count=start.count_table,size=893,persample=true)
classify.seqs(fasta=current, count=current, reference=silva.nr_v132.pcr.align, taxonomy=silva.nr_v132.tax, cutoff=80);
rename.file(fasta=start.subsample.fasta, count=start.subsample.count_table,prefix=893)
rename.file(input=start.subsample.nr_v132.wang.tax.summary,new=893start.subsample.nr_v132.wang.tax.summary)

###Alpha diversity analysis
summary.single(shared=start.droseu.unique.subsample.shared, calc=nseqs-coverage-sobs-chao-shannon-simpson-invsimpson-shannoneven-simpsoneven)

# assign OTU at 97% for comparing alpha diversity to previous study
make.shared(list=start.an.unique_list.list, count=start.count_table, label=0.03)
rename.file(input=start.an.unique_list.0.03.pick.shared,new=start.all.97.shared)
count.groups(shared=start.all.97.shared)
get.groups(shared=start.all.97.shared,groups=AT_Mau_14_01-UK_Sou_14_10-CY_Nic_14_11-UK_Mar_14_12-UK_Lut_14_13-DE_Bro_14_14-DE_Bro_14_15-UA_Yal_14_16-UA_Yal_14_17-UA_Yal_14_18-UA_Ode_14_19-AT_Mau_14_02-UA_Ode_14_20-UA_Ode_14_21-UA_Ode_14_22-UA_Kyi_14_23-UA_Kyi_14_24-UA_Var_14_25-UA_Pyr_14_26-UA_Dro_14_27-UA_Cho_14_28-UA_Cho_14_29-TR_Yes_14_03-SE_Lun_14_30-DE_Mun_14_31-DE_Mun_14_32-PT_Rec_14_33-ES_Gim_14_34-ES_Gim_14_35-FI_Aka_14_36-FI_Aka_14_37-FI_Ves_14_38-DK_Kar_14_39-TR_Yes_14_04-DK_Kar_14_40-DK_Kar_14_41-CH_Cha_14_42-CH_Cha_14_43-AT_See_14_44-UA_Kha_14_45-UA_Kha_14_46-UA_Cho_14_47-UA_Cho_14_48-UA_Kyi_14_49-FR_Vil_14_05-UA_Uma_14_50-FR_Vil_14_06-FR_Vil_14_07-FR_Got_14_08-UK_She_14_09)
rename.file(input=start.all.97.0.03.pick.shared,new=start.droseu97.shared)
sub.sample(shared=start.droseu97.shared,size=5135,persample=true)
summary.single(shared=start.droseu97.0.03.subsample.shared, calc=nseqs-coverage-sobs-chao-shannon-simpson-invsimpson-shannoneven-simpsoneven)
```

## B) download worldclim climate data
```R
####get climate data from worldclim
library(raster)
library(sp)

###get annual temperature and precipitation
r <- getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")

#load geographic cooridantes
xy=read.csv('DrosEUmetadata.csv',sep=',',header = T)
lats <- xy$lat
lons <- xy$long
coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points)
df <- cbind.data.frame(xy,values)
df

###get tmean data, monthly average temperature
t <- getData("worldclim",var="tmean",res=10)
points <- SpatialPoints(coords, proj4string = t@crs)
values <- extract(t,points)

month=c(7,10,8,7,8,9,10,7,8,7,8,10,10,6,10,6,7,8,7,7,8,10,8,9,8,8,8,9,9,7,6,9,9,10,8,7,8,7,9,9,11,7,10,8,7,9,9,9,10,10)

temp=c()
for (i in 1:nrow(values)) {
  temp=c(temp,as.vector(values[i,m[i]]))
}
temp # average monthly mean temperature (?C * 10)


p <- getData("worldclim",var="prec",res=10)
points <- SpatialPoints(coords, proj4string = p@crs)
values <- extract(p,points)
prec=c()
for (i in 1:nrow(values)) {
  prec=c(prec,as.vector(values[i,m[i]]))
}
prec # average monthly precipitation (mm)
```

## C) RDA & dbMEM analysis
```R
library(adespatial)
library(vegan)

#load meta table
meta<-read.csv("DrosEUmetadata.csv",header=T)
meta<-meta[order(meta$Sample),]
row.names(meta)<-c(1:length(meta[[1]]))

#exclude samples have no host genetic differentiation info
meta[c(-9,-35,-47),]->meta


######load and prepare OTU table (output from mothur analysis)
#100% OTU
shared100<-read.table("start.droseu.unique.subsample.shared",header=T)
shared100=shared100[,-1]
shared100=shared100[,-2]
shared100=data.frame(shared100[,-1],row.names = as.vector(shared100[,1]))
#plot columnsums
colSums(shared100)[1:100]
#take OTUs with at least 1000 reads
shared100EU<-shared100[,colSums(shared100[1:length(shared100)])>1000]

#remove samples for which there is no allel frequecy info
shared100EU[c(-9,-35,-47),]->shared100EU

# Hellinger transform the species table
shared100EU.h <- decostand (shared100EU, "hellinger")


########### dbMEM analysis of autocorrelation
#follow the protocol 
#Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, 2018
#chapter 7
library(SoDA)
droseu.xy=geoXY(meta$lat,meta$long)
droseu.xy=as.data.frame(droseu.xy,row.names = rownames(shared100EU.h))

#Is there a linear trend in the data?
anova(rda(shared100EU.h, droseu.xy))	# yes
# Computation of linearly detrended data
shared100EU.h.det <- resid(lm(as.matrix(shared100EU.h) ~ ., data = droseu.xy))
## Step 1. Construct the matrix of dbMEM variables
dbmem.tmp <- dbmem(droseu.xy, silent = FALSE)
dbmem <- as.data.frame(dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(droseu.xy)))

# Display and count the eigenvalues
attributes(dbmem.tmp)$values
length(attributes(dbmem.tmp)$values)
# Argument silent = FALSE allows the function to display 
# the truncation level.

## Step 2. Run the global dbMEM analysis on the *detrended*
##    Hellinger-transformed  data
dbmem.rda <- rda(shared100EU.h.det ~ ., dbmem)
anova(dbmem.rda)
#not significant globally


#########rda analysis######
##   at 100% OTU level 
#full model
Mfull.rda=rda(shared100EU.h~AF.PC1+AF.PC2+AF.PC3+tempmonth+precmonth+Precyear+Tempyear+substrate,data=meta)
#null model
M1000.rda<-rda(shared100EU.h~1,data=meta)

#run a global test first
anova(Mfull.rda, permutations = how(nperm = 999))  #0.001
anova(Mfull.rda, permutations = how(nperm = 999),by = "axis")   #3

#model selection
Mordi <- ordistep(M1000.rda, scope=formula(Mfull.rda), direction="forward",permutations = how(nperm = 999))
# AF.PC1 + substrate + AF.PC2 + Tempyear+ (tempmonth) were selected

#followed by backward selection
ordistep(Mordi,scope=formula(M1000.rda),perm.max=1000,direction="backward") # no change

#sense the model
extractAIC(Mordi) #12.00000 -19.90485
RsquareAdj(Mordi) #0.3623468  0.1771852

#examine VIF
vif.cca(Mordi)

#####sense the model with or without axes of host genetic differentiation
####Table S3
#best model from ordistep() but without axes of host genetic differentiation
Mordi.withoutAFPC1PC2 <- rda(shared100EU.h ~ substrate + Tempyear, data = meta)
#add either PC1 or PC2 to the model abova
Mordi.withoutAFPC2 <- rda(shared100EU.h ~ substrate + Tempyear + AF.PC1, data = meta)
Mordi.withoutAFPC1 <- rda(shared100EU.h ~ substrate + Tempyear + AF.PC2, data = meta)
#it is significantly different compared with a model without either PC1 or PC2
anova(Mordi, Mordi.withoutAFPC2)    
anova(Mordi, Mordi.withoutAFPC1)
extractAIC(Mordi.withoutAFPC2)
extractAIC(Mordi.withoutAFPC1)


####Table S4
#if we included latitude and longitude in the full model
#PC1 and PC2 were still selected
Mfull.rda=rda(shared100EU.h~AF.PC1+AF.PC2+AF.PC3+lat+long+tempmonth+precmonth+Precyear+Tempyear+substrate,data=meta)
M1000.rda<-rda(shared100EU.h~1,data=meta)
Mordi <- ordistep(M1000.rda, scope=formula(Mfull.rda), direction="forward",permutations = how(nperm = 999))
#AF.PC1 + substrate + lat + AF.PC2 + long + Precyear 
#it is significantly different compared with a model without PC1
Mordi.withoutAFPC1=rda(shared100EU.h ~ substrate + lat + AF.PC2 + long + Precyear, data = meta)
anova(Mordi, Mordi.withoutAFPC1) #0.001


##
#########rda analysis######  at family level  ##
#Table S5
#load tax.summary file from mothur analysis
data5203=read.table(file="5135start.subsample.nr_v132.wang.tax.summary",header=T) 
#pick taxlevel 5
data5203[data5203$taxlevel==5,]->phyla5203
#pick DrosEu samples
droseu=c("taxlevel","rankID","taxon","daughterlevels","total","D1","D10","D11","D12","D13","D14","D15","D16","D18","D19","D2","D20","D21","D22","D23","D24","D25","D26","D27","D28","D29","D3","D30","D31","D32","D33","D34","D35","D36","D37","D38","D39","D4","D41","D42","D43","D44","D45","D46","D47","D48","D49","D5","D50","D7","D8","D9")
dataeu=data5203[,droseu]
phylaeu=phyla5203[,droseu]

#prepare phyla table
phylaeu2=phylaeu[,-c(1,2,4,5)]
#take a look of sums
plot(phylaeu$total)
#take phyla with at least 1000 reads
phylaeu2=phylaeu2[rowSums(phylaeu2[,-1])>1000,]
rownames(phylaeu2)=phylaeu2$taxon
phylaeu2=phylaeu2[,-1]
phylaeu2=as.data.frame(t(phylaeu2))

## Hellinger transform the species table
phylaeu2.h <- decostand (phylaeu2, "hellinger")
#null model 
M1000.rda<-rda(phylaeu2.h~1,data=meta)
#full model
Mfull.rda=rda(phylaeu2.h~AF.PC1+AF.PC2+AF.PC3+tempmonth+precmonth+Precyear+Tempyear+substrate,data=meta)
#run a global test first
anova(Mfull.rda, permutations = how(nperm = 999)) #0.003

#model seletion
Mordi=ordistep(M1000.rda, scope=formula(Mfull.rda), direction="forward", permutations = how(nperm = 999))
#Tempyear + AF.PC1 were selected

#sense the model
extractAIC(Mordi) #3.00000 -77.23731
RsquareAdj(Mordi) #0.1466827  0.1078955

```

## D) R script for figure 1b and fig S3
```R
library(Imap)
library(vegan)

###read geographical matrix function
#function to get geographical distances
#from https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
    GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
    # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
    # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
    # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
    # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
    return(mat.distances)
}
###end of reading function 


#load meta table for latitude,longitude info
meta<-read.csv("DrosEUmetadata.csv",header=T)
meta<-meta[order(meta$Sample),]
row.names(meta)<-c(1:length(meta[[1]]))

#calculate geographic distance
df.cities <- data.frame(name =meta$Sample,  lat  = meta$lat, lon  = meta$long)
geodist=GeoDistanceInMetresMatrix(df.cities) / 1000


#load OTU table and calculate bray curtis matrix
shared100<-read.table("start.droseu.unique.subsample.shared",header=T) #output from mothur analysis
shared100=shared100[,-1]
shared100=shared100[,-2]
shared100=data.frame(shared100[,-1],row.names = as.vector(shared100[,1]))
braymatrix=as.matrix(vegdist(shared100))

#mantel test
mantel(braymatrix,geodist,permutations=9999)$statistic
mantel(braymatrix,geodist,permutations=9999)$signif


##fig 1b
library(ggplot2)
x=as.vector(as.dist(geodist))
y=as.vector(as.dist(braymatrix))
xy=as.data.frame(cbind(x,y))
ggplot(xy, aes(x=x, y=y)) +
  geom_point( shape=16,color=alpha('darkblue',0.35),cex=2.8)+
  #geom_smooth(method=lm, se=FALSE,linetype="dashed",color="red",size=1)+
  labs(x = "geographic distance in km",y = "bacterial community dissimilarity")+
  theme(axis.text=element_text(size=16,face="bold",colour = "black"),axis.title=element_text(size=17.5,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.grid.major =element_line(colour = 'grey'))+
  scale_y_continuous(limits = c(0.2, 1))+
  annotate(geom="text", x=2800, y=0.3, label="\nP = 0.0015\nr = 0.20\n",
           color="black",size=6,fontface="bold")+
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE,size=1)

##fig S3a
#read pariwise fst table
fst=read.csv("DrosEU-FST_Geo_update.csv",header = T)
fst=fst[order(fst$Population1),]
fst=fst[fst$Population1!='s51',]
fst=fst[fst$Population2!='s51',]
droplevels(fst)->fst
fst.mat=as.matrix(xtabs(fst[, 3] ~ fst[, 2] + fst[, 1]))
s1=c(0,rep(0,dim(fst.mat)[2]-1))
names(s1)=colnames(fst.mat)
as.dist(rbind(s1,fst.mat))

#remove samples in meta table and shared file without fst infomation
meta[c(-9,-35,-47),]->meta
df.cities <- data.frame(name =meta$Sample,  lat  = meta$lat, lon  = meta$long)
geodist=GeoDistanceInMetresMatrix(df.cities) / 1000
shared100[c(-9,-35,-47),]->shared100
braymatrix=as.matrix(vegdist(shared100))

#do partial mantel test
mantel.partial(as.dist(braymatrix),as.dist(rbind(s1,fst.mat)),as.dist(geodist),permutations=9999)$statistic 
#0.1855181
mantel.partial(as.dist(braymatrix),as.dist(rbind(s1,fst.mat)),as.dist(geodist),permutations=9999)$signif 
#0.0211

#plot fig S3a
x=as.vector(as.dist(rbind(s1,fst.mat)))
y=as.vector(resid(lm(as.vector(as.dist(braymatrix))~as.vector(as.dist(geodist)))))
xy=as.data.frame(cbind(x,y))
ggplot(xy, aes(x=x, y=y)) +
  geom_point( shape=16,color=alpha('darkblue',0.3),cex=2.8)+
  #geom_smooth(method=lm, se=FALSE,linetype="dashed",color="red",size=1)+
  # labs(x = expression("pairwise F"[ST]),y = "residual community dissimilarity")+
  labs(x ="pairwise Fst",y = "residual community dissimilarity")+
  theme(axis.text=element_text(size=15,face="bold",colour = "black"),axis.title=element_text(size=16,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.grid.major =element_line(colour = 'grey'))+
  # scale_y_continuous(limits = c(0.2, 1))+
  annotate(geom="text", x=0.0375, y=-0.4, label="\nP = 0.02\nr = 0.19\n",
           color="black",size=6,fontface="bold")+
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE,size=1.25,col='black')


######using different distance metrics for
#fig S3bcde
#jaccard index
jmatrix=as.matrix(vegdist(shared100,method="jaccard"))
mantel(jmatrix,geodist,permutations=9999)$statistic  #0.2040113
mantel(jmatrix,geodist,permutations=9999)$signif #8e-04
#partial mantel test
mantel.partial(as.dist(jmatrix),as.dist(rbind(s1,fst.mat)),as.dist(geodist),permutations=9999)$statistic #0.1895445
mantel.partial(as.dist(jmatrix),as.dist(rbind(s1,fst.mat)),as.dist(geodist),permutations=9999)$signif #0.0168

#morisita-horn index
hmatrix=as.matrix(vegdist(shared100,method="horn"))
mantel(hmatrix,geodist,permutations=9999)$statistic  #0.1827412
mantel(hmatrix,geodist,permutations=9999)$signif #0.001
#partial mantel test
mantel.partial(as.dist(hmatrix),as.dist(rbind(s1,fst.mat)),as.dist(geodist),permutations=9999)$statistic #0.1622406
mantel.partial(as.dist(hmatrix),as.dist(rbind(s1,fst.mat)),as.dist(geodist),permutations=9999)$signif #0.0267

x=as.vector(as.dist(geodist))
#figs3b
png("fig.mantel.jaccard.png", width = 550, height = 500,res = 100)
y=as.vector(as.dist(jmatrix))
xy=as.data.frame(cbind(x,y))
ggplot(xy, aes(x=x, y=y)) +
  geom_point( shape=16,color=alpha('darkblue',0.3),cex=2.8)+
  #geom_smooth(method=lm, se=FALSE,linetype="dashed",color="red",size=1)+
  labs(title="jaccard",x = "geographic distance in km",y = "bacterial community dissimilarity")+
  theme(plot.title = element_text(size=17.5,face="bold",colour = "black",hjust = 0.5),
        axis.text=element_text(size=16,face="bold",colour = "black"),
        axis.title=element_text(size=17.5,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.grid.major =element_line(colour = 'grey'))+
  scale_y_continuous(limits = c(0.2, 1))+
  annotate(geom="text", x=2800, y=0.3, label="\nP = 0.0008\nr = 0.20\n",
           color="black",size=6,fontface="bold")+
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE,size=1.25,col='black')
dev.off()
#figure s3d
png("fig.mantel.horn-morisita.png", width = 550, height = 500,res = 100)
y=as.vector(as.dist(hmatrix))
xy=as.data.frame(cbind(x,y))
ggplot(xy, aes(x=x, y=y)) +
  geom_point( shape=16,color=alpha('darkblue',0.3),cex=2.8)+
  #geom_smooth(method=lm, se=FALSE,linetype="dashed",color="red",size=1)+
  labs(title="horn-morisita",x = "geographic distance in km",y = "bacterial community dissimilarity")+
  theme(plot.title = element_text(size=17.5,face="bold",colour = "black",hjust = 0.5),
        axis.text=element_text(size=16,face="bold",colour = "black"),
        axis.title=element_text(size=17.5,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.grid.major =element_line(colour = 'grey'))+
  scale_y_continuous(limits = c(0.2, 1))+
  annotate(geom="text", x=2800, y=0.3, label="\nP = 0.0010\nr = 0.18\n",
           color="black",size=6,fontface="bold")+
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE,size=1.25,col='black')
dev.off()
#residual plot
x=as.vector(as.dist(rbind(s1,fst.mat)))
#figure s3c
png("fig.mantel.jaccard.resid.png", width = 550, height = 500,res = 100)
y=as.vector(resid(lm(as.vector(as.dist(jmatrix))~as.vector(as.dist(geodist)))))
xy=as.data.frame(cbind(x,y))
ggplot(xy, aes(x=x, y=y)) +
  geom_point( shape=16,color=alpha('darkblue',0.3),cex=2.8)+
  #geom_smooth(method=lm, se=FALSE,linetype="dashed",color="red",size=1)+
  # labs(x = expression("pairwise F"[ST]),y = "residual community dissimilarity")+
  labs(title="jaccard",x ="pairwise Fst",y = "residual community dissimilarity")+
  theme(plot.title = element_text(size=17.5,face="bold",colour = "black",hjust = 0.5),
        axis.text=element_text(size=15,face="bold",colour = "black"),
        axis.title=element_text(size=16,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.grid.major =element_line(colour = 'grey'))+
  # scale_y_continuous(limits = c(0.2, 1))+
  annotate(geom="text", x=0.0375, y=-0.4, label="\nP = 0.017\nr = 0.19\n",
           color="black",size=6,fontface="bold")+
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE,size=1.25,col='black')
dev.off()
#figure s3e
png("fig.mantel.horn-morisita.resid.png", width = 550, height = 500,res = 100)
y=as.vector(resid(lm(as.vector(as.dist(hmatrix))~as.vector(as.dist(geodist)))))
xy=as.data.frame(cbind(x,y))
ggplot(xy, aes(x=x, y=y)) +
  geom_point( shape=16,color=alpha('darkblue',0.3),cex=2.8)+
  #geom_smooth(method=lm, se=FALSE,linetype="dashed",color="red",size=1)+
  # labs(x = expression("pairwise F"[ST]),y = "residual community dissimilarity")+
  labs(title="horn-morisita",x ="pairwise Fst",y = "residual community dissimilarity")+
  theme(plot.title = element_text(size=17.5,face="bold",colour = "black",hjust = 0.5),
        axis.text=element_text(size=15,face="bold",colour = "black"),
        axis.title=element_text(size=16,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_rect(fill = 'white', colour = 'red'),
        panel.grid.major =element_line(colour = 'grey'))+
  # scale_y_continuous(limits = c(0.2, 1))+
  annotate(geom="text", x=0.0375, y=-0.5, label="\nP = 0.027\nr = 0.16\n",
           color="black",size=6,fontface="bold")+
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE,size=1.25,col='black')
dev.off()
```

## References

