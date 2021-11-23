library("data.table")
library("epiR")
library("ggplot2")
library("ape")
library("raster")
library("spdep")
library("maps")
library("prevR")
library("ggplotify")
library("gridExtra")
library("rgdal") 
library("maptools")
library("plyr")
library("ggsn")
library("cowplot")

# Read Data ---------------------------------------------------------------

#Read qPCR results.
vir <- read.csv("../Data/Viruses_data.csv")
vir$Mean.Qty <- nafill(vir$Mean.Qty, fill=0)
vir$Mean.Qty[vir$Mean.Qty < 1000] <- 0
vir$Well <- NULL
vir$Task <- NULL
vir$Ct<- NULL
vir$StdDev.Ct<- NULL
vir$Qty<- NULL
vir$StdDev.Qty<- NULL
vir$Tm<- NULL
vir$Filtered<- NULL
colnames(vir) <- c("Sample", "Virus","ViralLoad")


#Read metadata
Metadat <- read.csv("../Data/Genotype_location.csv")
Metadat$Index <- NULL
Metadat$CODA <- NULL
Metadat$Werklijst <- NULL
Metadat$Ugentnummer <- NULL
Metadat$Staalnummer <- NULL
Metadat$FAVV <- NULL
Metadat$Colony.Nr. <- NULL
Metadat$Postcode <- NULL
Virdf_metadata <- merge(vir, Metadat, by="Sample")
write.csv(Virdf_metadata, file="../Data/Virdf_metadata.csv", quote=FALSE)
AMFV <- Virdf_metadata[Virdf_metadata$Virus == 'AMFV',]
BMLV <- Virdf_metadata[Virdf_metadata$Virus == 'bmlv',]
thog <- Virdf_metadata[Virdf_metadata$Virus == 'thog',]
thika <- Virdf_metadata[Virdf_metadata$Virus == 'thika',]
unc <- Virdf_metadata[Virdf_metadata$Virus == 'Unc',]
vdv <- Virdf_metadata[Virdf_metadata$Virus == 'vdv',]
dwv <- Virdf_metadata[Virdf_metadata$Virus == 'dwv',]


# Prevalences -------------------------------------------------------------

#Estimate prevalences:
AMFVprev <- epi.prev(pos = sum(AMFV$'ViralLoad' > 0), tested = length(AMFV$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
BMLVprev <- epi.prev(pos = sum(BMLV$'ViralLoad' > 0), tested = length(BMLV$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
thogprev <- epi.prev(pos = sum(thog$'ViralLoad' > 0), tested = length(thog$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
thikaprev <- epi.prev(pos = sum(thika$'ViralLoad' > 0), tested = length(thika$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
uncrev <- epi.prev(pos = sum(unc$'ViralLoad' > 0), tested = length(unc$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
vdvprev <- epi.prev(pos = sum(vdv$'ViralLoad' > 0), tested = length(vdv$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
dwvprev <- epi.prev(pos = sum(dwv$'ViralLoad' > 0), tested = length(dwv$'ViralLoad'), se=0.95, sp=0.95, method="blaker", conf.level=0.95)
#Throw together results:
estprev <- t(rbind(AMFVprev$tp, BMLVprev$tp,thogprev$tp,thikaprev$tp,uncrev$tp,vdvprev$tp,dwvprev$tp))
colnames(estprev) <- c("Apis Melifera filamentous virus","Bee-macula like virus","Apis Orthomyxovirus 1","Apis Thika-like virus","Apis Uncultured virus","VDV","DWV")
estprev <- t(estprev)
estprev <- as.data.frame(estprev)
#Plot prevalences:
prevalenceplot <- ggplot(estprev, aes(x=reorder(rownames(estprev), -est), y=est)) + 
  geom_errorbar(aes(ymin=estprev$lower, ymax=estprev$upper), colour="black", width=.1) +
  geom_bar(stat = "identity",color="black",fill="lightgray") + 
  labs(x = "Virus",y="Estimated Prevalence (%)") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2))
prevalenceplot
ggsave("Plots/EstPrevalences.pdf",prevalenceplot, dpi=300)


# Moran's I test ----------------------------------------------------------

##Ape
BEL <- raster::getData("GADM", country="BEL", level=4)
FLA <- BEL[BEL$NAME_1 == "Vlaanderen",]
#Calculate all:
##AMFV
##Aggregate per commune:
AMFV_comm <- aggregate(ViralLoad~Commune, data=AMFV, FUN=function(x) c(median=median(x), count=length(x)))
AMFV_FLA <- merge(FLA,AMFV_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
AMFV_distmat <- as.matrix(dist(cbind(coordinates(AMFV_FLA)[,1], coordinates(AMFV_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
AMFV_distmatinv <- 1/AMFV_distmat
diag(AMFV_distmatinv) <- 0
AMFVmoran <- Moran.I(AMFV_FLA$ViralLoad[,1], AMFV_distmatinv, na.rm=TRUE, scaled=TRUE)
AMFVmoranPRES <- Moran.I(ifelse(AMFV_FLA$ViralLoad[,1] > 1000,1,0), AMFV_distmatinv, na.rm=TRUE, scaled=TRUE)
##BMLV
BMLV_comm <- aggregate(ViralLoad~Commune, data=BMLV, FUN=function(x) c(median=median(x), count=length(x)))
BMLV_FLA <- merge(FLA,BMLV_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
BMLV_distmat <- as.matrix(dist(cbind(coordinates(BMLV_FLA)[,1], coordinates(BMLV_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
BMLV_distmatinv <- 1/BMLV_distmat
diag(BMLV_distmatinv) <- 0
BMLVmoran <- Moran.I(BMLV_FLA$ViralLoad[,1], BMLV_distmatinv, na.rm=TRUE, scaled=TRUE)
BMLVmoranPRES <- Moran.I(ifelse(BMLV_FLA$ViralLoad[,1] > 1000, 1,0), BMLV_distmatinv, na.rm=TRUE, scaled=TRUE)

#thog
thog_comm <- aggregate(ViralLoad~Commune, data=thog, FUN=function(x) c(median=median(x), count=length(x)))
thog_FLA <- merge(FLA,thog_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
thog_distmat <- as.matrix(dist(cbind(coordinates(thog_FLA)[,1], coordinates(thog_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
thog_distmatinv <- 1/thog_distmat
diag(thog_distmatinv) <- 0
thogmoran <- Moran.I(thog_FLA$ViralLoad[,1], thog_distmatinv, na.rm=TRUE, scaled=TRUE)
thogmoranPRES <- Moran.I(ifelse(thog_FLA$ViralLoad[,1] > 1000,1,0), thog_distmatinv, na.rm=TRUE, scaled=TRUE)

##thika
thika_comm <- aggregate(ViralLoad~Commune, data=thika, FUN=function(x) c(median=median(x), count=length(x)))
thika_FLA <- merge(FLA,thika_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
thika_distmat <- as.matrix(dist(cbind(coordinates(thika_FLA)[,1], coordinates(thika_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
thika_distmatinv <- 1/thika_distmat
diag(thika_distmatinv) <- 0
thikamoran <- Moran.I(thika_FLA$ViralLoad[,1], thika_distmatinv, na.rm=TRUE, scaled=TRUE)
thikamoranPRES <- Moran.I(ifelse(thika_FLA$ViralLoad[,1] > 1000, 1,0), thika_distmatinv, na.rm=TRUE, scaled=TRUE)


##unc
unc_comm <- aggregate(ViralLoad~Commune, data=unc, FUN=function(x) c(median=median(x), count=length(x)))
unc_FLA <- merge(FLA,unc_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
unc_distmat <- as.matrix(dist(cbind(coordinates(unc_FLA)[,1], coordinates(unc_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
unc_distmatinv <- 1/unc_distmat
diag(unc_distmatinv) <- 0
uncmoran <- Moran.I(unc_FLA$ViralLoad[,1], unc_distmatinv, na.rm=TRUE, scaled=TRUE)
uncmoranPRES <- Moran.I(ifelse(unc_FLA$ViralLoad[,1] > 1000,1,0), unc_distmatinv, na.rm=TRUE, scaled=TRUE)

##vdv
vdv_comm <- aggregate(ViralLoad~Commune, data=vdv, FUN=function(x) c(median=median(x), count=length(x)))
vdv_FLA <- merge(FLA,vdv_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
vdv_distmat <- as.matrix(dist(cbind(coordinates(vdv_FLA)[,1], coordinates(vdv_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
vdv_distmatinv <- 1/vdv_distmat
diag(vdv_distmatinv) <- 0
vdvmoran <- Moran.I(vdv_FLA$ViralLoad[,1], vdv_distmatinv, na.rm=TRUE, scaled=TRUE)
vdvmoranPRES <- Moran.I(ifelse(vdv_FLA$ViralLoad[,1] > 1000,1,0), vdv_distmatinv, na.rm=TRUE, scaled=TRUE)
##dwv
dwv_comm <- aggregate(ViralLoad~Commune, data=dwv, FUN=function(x) c(median=median(x), count=length(x)))
dwv_FLA <- merge(FLA,dwv_comm,by.x="NAME_4",by.y="Commune")
#Create distances matrix
dwv_distmat <- as.matrix(dist(cbind(coordinates(dwv_FLA)[,1], coordinates(dwv_FLA)[,2])))
#Invert distance matrix + put diagonal at 0
dwv_distmatinv <- 1/dwv_distmat
diag(dwv_distmatinv) <- 0
dwvmoran <- Moran.I(dwv_FLA$ViralLoad[,1], dwv_distmatinv, na.rm=TRUE, scaled=TRUE)
dwvmoranPRES <- Moran.I(ifelse(dwv_FLA$ViralLoad[,1] > 1000,1,0), dwv_distmatinv, na.rm=TRUE, scaled=TRUE)
#Wrap results together:
moranAPE <- rbind(c(AMFVmoran$observed[[1]],AMFVmoran$expected[[1]],AMFVmoran$sd[[1]],AMFVmoran$p.value,"APE"),
                  c(BMLVmoran$observed[[1]],BMLVmoran$expected[[1]],BMLVmoran$sd[[1]],BMLVmoran$p.value,"APE"),
                  c(thogmoran$observed[[1]],thogmoran$expected[[1]],thogmoran$sd[[1]],thogmoran$p.value,"APE"),
                  c(thikamoran$observed[[1]],thikamoran$expected[[1]],thikamoran$sd[[1]],thikamoran$p.value,"APE"),
                  c(uncmoran$observed[[1]],uncmoran$expected[[1]],uncmoran$sd[[1]],uncmoran$p.value,"APE"),
                  c(vdvmoran$observed[[1]],vdvmoran$expected[[1]],vdvmoran$sd[[1]],vdvmoran$p.value,"APE"),
                  c(dwvmoran$observed[[1]],dwvmoran$expected[[1]],dwvmoran$sd[[1]],dwvmoran$p.value,"APE"))
moranAPE <- as.data.frame(moranAPE)
moranAPE$virus <- c("Apis Melifera filamentous virus","Bee-macula like virus","Apis Orthomyxovirus 1","Apis Thika-like virus","Apis Uncultured virus","VDV","DWV")
colnames(moranAPE) <- c("observed","expected","sd","p.value","test","virus")
#Wrap results together:
moranPRESAPE <- rbind(c(AMFVmoranPRES$observed[[1]],AMFVmoranPRES$expected[[1]],AMFVmoranPRES$sd[[1]],AMFVmoranPRES$p.value,"APE"),
                      c(BMLVmoranPRES$observed[[1]],BMLVmoranPRES$expected[[1]],BMLVmoranPRES$sd[[1]],BMLVmoranPRES$p.value,"APE"),
                      c(thogmoranPRES$observed[[1]],thogmoranPRES$expected[[1]],thogmoranPRES$sd[[1]],thogmoranPRES$p.value,"APE"),
                      c(thikamoranPRES$observed[[1]],thikamoranPRES$expected[[1]],thikamoranPRES$sd[[1]],thikamoranPRES$p.value,"APE"),
                      c(uncmoranPRES$observed[[1]],uncmoranPRES$expected[[1]],uncmoranPRES$sd[[1]],uncmoranPRES$p.value,"APE"),
                      c(vdvmoranPRES$observed[[1]],vdvmoranPRES$expected[[1]],vdvmoranPRES$sd[[1]],vdvmoranPRES$p.value,"APE"),
                      c(dwvmoranPRES$observed[[1]],dwvmoranPRES$expected[[1]],dwvmoranPRES$sd[[1]],dwvmoranPRES$p.value,"APE"))
moranPRESAPE <- as.data.frame(moranPRESAPE)
moranPRESAPE$virus <- c("Apis Melifera filamentous virus","Bee-macula like virus","Apis Orthomyxovirus 1","Apis Thika-like virus","Apis Uncultured virus","VDV","DWV")
colnames(moranPRESAPE) <- c("observed","expected","sd","p.value","test","virus")

#SPDEP
#Test Moran's I with SPDEP
#use distances rather then number of neighbours.
#Put distances at 10000m (https://naldc.nal.usda.gov/download/IND43968380/PDF)

coordsfla <- coordinates(FLA)
poldists <- dnearneigh(coordsfla, 0, 10000)
lw <- nb2listw(poldists, style="W",zero.policy=T) 
nb <- poly2nb(FLA, queen=TRUE)
lw <- nb2listw(nb, style="S",zero.policy=T)
#Moran's I tests

AMFV_loads <- merge(FLA,AMFV_comm,by.x="NAME_4",by.y="Commune")
AMFV_moranspdep <- moran.test(AMFV_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
AMFV_moranspdepPRES <- moran.test(ifelse(AMFV_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
BMLV_loads <- merge(FLA,BMLV_comm,by.x="NAME_4",by.y="Commune")
BMLV_moranspdep <- moran.test(BMLV_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
BMLV_moranspdepPRES <- moran.test(ifelse(BMLV_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
thog_loads <- merge(FLA,thog_comm,by.x="NAME_4",by.y="Commune")
thog_moranspdep <- moran.test(thog_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
thog_moranspdepPRES <- moran.test(ifelse(thog_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
thika_loads <- merge(FLA,thika_comm,by.x="NAME_4",by.y="Commune")
thika_moranspdep <- moran.test(thika_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
thika_moranspdepPRES <- moran.test(ifelse(thika_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
unc_loads <- merge(FLA,unc_comm,by.x="NAME_4",by.y="Commune")
unc_moranspdep <- moran.test(unc_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
unc_moranspdepPRES <- moran.test(ifelse(unc_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
vdv_loads <- merge(FLA,vdv_comm,by.x="NAME_4",by.y="Commune")
vdv_moranspdep <- moran.test(vdv_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
vdv_moranspdepPRES <- moran.test(ifelse(vdv_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
dwv_loads <- merge(FLA,dwv_comm,by.x="NAME_4",by.y="Commune")
dwv_moranspdep <- moran.test(dwv_loads$ViralLoad[,1],lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)
dwv_moranspdepPRES <- moran.test(ifelse(dwv_loads$ViralLoad[,1] > 1000,1,0),lw, zero.policy=TRUE, na.action=na.exclude,randomisation=FALSE)

moranSPDEP <- rbind(c(AMFV_moranspdep$estimate[[1]],AMFV_moranspdep$estimate[[2]],abs(AMFV_moranspdep$statistic[[1]]),AMFV_moranspdep$p.value,"SPDEP"),
                    c(BMLV_moranspdep$estimate[[1]],BMLV_moranspdep$estimate[[2]],abs(BMLV_moranspdep$statistic[[1]]),BMLV_moranspdep$p.value,"SPDEP"),
                    c(thog_moranspdep$estimate[[1]],thog_moranspdep$estimate[[2]],abs(thog_moranspdep$statistic[[1]]),thog_moranspdep$p.value,"SPDEP"),
                    c(thika_moranspdep$estimate[[1]],thika_moranspdep$estimate[[2]],abs(thika_moranspdep$statistic[[1]]),thika_moranspdep$p.value,"SPDEP"),
                    c(unc_moranspdep$estimate[[1]],unc_moranspdep$estimate[[2]],abs(unc_moranspdep$statistic[[1]]),unc_moranspdep$p.value,"SPDEP"),
                    c(vdv_moranspdep$estimate[[1]],vdv_moranspdep$estimate[[2]],abs(vdv_moranspdep$statistic[[1]]),vdv_moranspdep$p.value,"SPDEP"),
                    c(dwv_moranspdep$estimate[[1]],dwv_moranspdep$estimate[[2]],abs(dwv_moranspdep$statistic[[1]]),dwv_moranspdep$p.value,"SPDEP"))
moranSPDEP <- as.data.frame(moranSPDEP)
colnames(moranSPDEP) <- c("observed","expected","sd","p.value","test")
moranSPDEP$virus <- c("Apis Melifera filamentous virus","Bee-macula like virus","Apis Orthomyxovirus 1","Apis Thika-like virus","Apis Uncultured virus","VDV","DWV")

moranSPDEPPRES <- rbind(c(AMFV_moranspdep$estimate[[1]],AMFV_moranspdep$estimate[[2]],abs(AMFV_moranspdep$statistic[[1]]),AMFV_moranspdep$p.value,"SPDEP"),
                    c(BMLV_moranspdep$estimate[[1]],BMLV_moranspdep$estimate[[2]],abs(BMLV_moranspdep$statistic[[1]]),BMLV_moranspdep$p.value,"SPDEP"),
                    c(thog_moranspdep$estimate[[1]],thog_moranspdep$estimate[[2]],abs(thog_moranspdep$statistic[[1]]),thog_moranspdep$p.value,"SPDEP"),
                    c(thika_moranspdep$estimate[[1]],thika_moranspdep$estimate[[2]],abs(thika_moranspdep$statistic[[1]]),thika_moranspdep$p.value,"SPDEP"),
                    c(unc_moranspdep$estimate[[1]],unc_moranspdep$estimate[[2]],abs(unc_moranspdep$statistic[[1]]),unc_moranspdep$p.value,"SPDEP"),
                    c(vdv_moranspdep$estimate[[1]],vdv_moranspdep$estimate[[2]],abs(vdv_moranspdep$statistic[[1]]),vdv_moranspdep$p.value,"SPDEP"),
                    c(dwv_moranspdep$estimate[[1]],dwv_moranspdep$estimate[[2]],abs(dwv_moranspdep$statistic[[1]]),dwv_moranspdep$p.value,"SPDEP"))
moranSPDEPPRES <- as.data.frame(moranSPDEPPRES)
colnames(moranSPDEPPRES) <- c("observed","expected","sd","p.value","test")
moranSPDEPPRES$virus <- c("Apis Melifera filamentous virus","Bee-macula like virus","Apis Orthomyxovirus 1","Apis Thika-like virus","Apis Uncultured virus","VDV","DWV")

Moran_results <- rbind(moranAPE,moranSPDEP)
Moran_results$observed <- as.numeric(as.character(Moran_results$observed))
Moran_results$expected <- as.numeric(as.character(Moran_results$expected))
Moran_results$sd <- as.numeric(as.character(Moran_results$sd))
Moran_results$p.value <- as.numeric(as.character(Moran_results$p.value))

Moranplot <- ggplot(Moran_results, aes(x=virus, y=observed, fill=test)) + 
  scale_y_continuous(limits = c(-1, 1)) +
  geom_errorbar(position=position_dodge(width=1),alpha=0.3,aes(ymin=observed-sd, ymax=observed+sd)) +
  geom_bar(position=position_dodge(width=1), stat="identity",color="black") +
  geom_text(aes(label=c("-0.008","-0.015","-0.009","-0.027","-0.014",NA,"-0.028","-0.024","-0.091",NA,"-0.017","-0.030",NA,"-0.127")), position=position_dodge(width=1), hjust=1.1, vjust=-0.2) +
  geom_text(aes(label=c(NA,NA,NA,NA,NA,0.008,NA,NA,NA,0.011,NA,NA,0.043,NA)), position=position_dodge(width=1), hjust=-0.2, vjust=-0.2) +
 #geom_text(aes(y=c(-0.057607382,-0.064547090,-0.058624241,-0.076606636,-0.063619480,0.058069378,-0.077840934,-0.073673040,-0.140963243,0.061485294,-0.066587787,-0.080439125,0.093336006,-0.176663343),label=round(observed,3)), position=position_dodge(width=1), hjust=1) +
  labs(x = "Virus",y="Moran's I statistic") +
  scale_fill_manual(values=c("lightgray","black")) +
  coord_flip() +
  theme_minimal()
Moranplot
ggsave("Plots/MoranI.pdf",Moranplot, dpi=300)

# KDE estimations ---------------------------------------------------------
FLAcommcoor <- as.data.frame(cbind(coordinates(FLA),FLA$NAME_4))
colnames(FLAcommcoor) <- c("LON","LAT","Commune")
AMFVprev <- aggregate(ViralLoad~Commune, data=AMFV, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
AMFVprev2 <- as.data.frame(cbind(AMFVprev$ViralLoad[,1], AMFVprev$ViralLoad[,2]))
colnames(AMFVprev2) <- c("pos","num")
AMFVprev2$Commune <- AMFVprev$Commune
AMFVprevdf <- merge(FLAcommcoor,AMFVprev2,by="Commune")

FLAgg <- FLA
FLAgg@data$id <- rownames(FLAgg@data)
FLAgg.points <- fortify(FLAgg, region="id")
FLAgg.df <- join(FLAgg.points, FLAgg@data, by="id")

FLAsamples<- ggplot(FLAgg.df) + 
  aes(long,lat,group=group) + 
  geom_polygon(fill="#D3D3D3") +
  geom_path(color="white", size=0.05) +
  geom_point(data=AMFVprevdf, aes(x=as.numeric(as.character(AMFVprevdf$LON)), as.numeric(as.character(AMFVprevdf$LAT))), inherit.aes = FALSE, size = AMFVprevdf$num/10) + 
  geom_point(aes(x=c(5.9),y=c(51.1)), inherit.aes = FALSE, size=1/10) +
  geom_point(aes(x=c(5.9),y=c(50.9)), inherit.aes = FALSE, size=10/10) +
  #geom_text(aes(x=c(5.97), y=c(51), label="1"), inherit.aes = FALSE) +
  #geom_text(aes(x=c(5.97), y=c(50.9), label="10"), inherit.aes = FALSE) +
  annotate(geom="text", x=6.0, y=51.1, label="1",color="black") + 
  annotate(geom="text", x=6.0, y=50.9, label="10",color="black") + 
  theme_minimal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) +
  scale_colour_manual(values=c("#000000")) + 
  labs(title='Samples Flanders') +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_equal() +
  ggsn::scalebar(FLAgg.df, transform = TRUE, dist_unit = "km", dist=10,height=0.01, st.dist=0.5, st.size=0.1) +
  annotate(geom="text", x=5.75, y=50.55, label="10 km",color="black")
FLAsamples

#plot(FLA, border="white",col="lightgray",main="Samples")
#points(as.numeric(as.character(AMFVprevdf$LON)), as.numeric(as.character(AMFVprevdf$LAT)), pch=16, cex=AMFVprevdf$num/5)
#points(5.9,51,cex=0.1,pch=16)
#points(5.9,50.9,cex=2,pch=16)
#text(5.97,51,"1")
#text(5.97,50.9,"10")
#map.scale(ratio=FALSE, metric=TRUE)
#flasampleplot <- recordPlot()

#Create needed df:
FLAcommcoor <- as.data.frame(cbind(coordinates(FLA),FLA$NAME_4))
colnames(FLAcommcoor) <- c("LON","LAT","Commune")
AMFVprev <- aggregate(ViralLoad~Commune, data=AMFV, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
AMFVprev2 <- as.data.frame(cbind(AMFVprev$ViralLoad[,1], AMFVprev$ViralLoad[,2]))
colnames(AMFVprev2) <- c("pos","num")
AMFVprev2$Commune <- AMFVprev$Commune
AMFVprevdf <- merge(FLAcommcoor,AMFVprev2,by="Commune")
AMFVprevdf$cluster <- rownames(AMFVprevdf)
AMFVprevdf$LON <- as.numeric(as.character(AMFVprevdf$LON))
AMFVprevdf$LAT <- as.numeric(as.character(AMFVprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
AMFV_prevr <- as.prevR(AMFVprevdf, col, FLA)
AMFV_prevr_r <- rings(AMFV_prevr, N = c(68,100,200), progression = FALSE)
AMFV.N100 <- kde(AMFV_prevr_r, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
AMFVres <- as.data.frame(AMFV.N100)
AMFVres <- AMFVres[!is.na(AMFVres$k.prev.N100.RInf), ]
AMFVkdeplot <- ggplot(data = AMFVres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='Apis Melifera filamentous virus',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
AMFVkdeplot
#ggsave("Plots/AMFVkde.pdf",AMFVkdeplot, dpi=300)



BMLVprev <- aggregate(ViralLoad~Commune, data=BMLV, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
BMLVprev2 <- as.data.frame(cbind(BMLVprev$ViralLoad[,1], BMLVprev$ViralLoad[,2]))
colnames(BMLVprev2) <- c("pos","num")
BMLVprev2$Commune <- BMLVprev$Commune
BMLVprevdf <- merge(FLAcommcoor,BMLVprev2,by="Commune")
BMLVprevdf$cluster <- rownames(BMLVprevdf)
BMLVprevdf$LON <- as.numeric(as.character(BMLVprevdf$LON))
BMLVprevdf$LAT <- as.numeric(as.character(BMLVprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
BMLV_prevr <- as.prevR(BMLVprevdf, col, FLA)
BMLV_prevr_R <- rings(BMLV_prevr, N = c(68,100,200), progression = FALSE)
BMLV.N100 <- kde(BMLV_prevr_R, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
BMLVres <- as.data.frame(BMLV.N100)
BMLVres <- BMLVres[!is.na(BMLVres$k.prev.N100.RInf), ]
BMLVkdeplot <- ggplot(data = BMLVres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='Bee-macula like virus',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
BMLVkdeplot
#ggsave("Plots/BMLVkde.pdf",BMLVkdeplot, dpi=300)


thogprev <- aggregate(ViralLoad~Commune, data=thog, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
thogprev2 <- as.data.frame(cbind(thogprev$ViralLoad[,1], thogprev$ViralLoad[,2]))
colnames(thogprev2) <- c("pos","num")
thogprev2$Commune <- thogprev$Commune
thogprevdf <- merge(FLAcommcoor,thogprev2,by="Commune")
thogprevdf$cluster <- rownames(thogprevdf)
thogprevdf$LON <- as.numeric(as.character(thogprevdf$LON))
thogprevdf$LAT <- as.numeric(as.character(thogprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
thog_prevr <- as.prevR(thogprevdf, col, FLA)
thog_prevr_R <- rings(thog_prevr, N = c(68,100,200), progression = FALSE)
thog.N100 <- kde(thog_prevr_R, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
thogres <- as.data.frame(thog.N100)
thogres <- thogres[!is.na(thogres$k.prev.N100.RInf), ]
thogkdeplot <- ggplot(data = thogres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='Apis Orthomyxovirus 1',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
thogkdeplot
#ggsave("Plots/thogkde.pdf",thogkdeplot, dpi=300)

thikaprev <- aggregate(ViralLoad~Commune, data=thika, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
thikaprev2 <- as.data.frame(cbind(thikaprev$ViralLoad[,1], thikaprev$ViralLoad[,2]))
colnames(thikaprev2) <- c("pos","num")
thikaprev2$Commune <- thikaprev$Commune
thikaprevdf <- merge(FLAcommcoor,thikaprev2,by="Commune")
thikaprevdf$cluster <- rownames(thikaprevdf)
thikaprevdf$LON <- as.numeric(as.character(thikaprevdf$LON))
thikaprevdf$LAT <- as.numeric(as.character(thikaprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
thika_prevr <- as.prevR(thikaprevdf, col, FLA)
thika_prevr_R <- rings(thika_prevr, N = c(68,100,200), progression = FALSE)
thika.N100 <- kde(thika_prevr_R, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
thikares <- as.data.frame(thika.N100)
thikares <- thikares[!is.na(thikares$k.prev.N100.RInf), ]
thikakdeplot <- ggplot(data = thikares) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='Apis Thika-like virus',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
thikakdeplot
#ggsave("Plots/thikakde.pdf",thikakdeplot, dpi=300)




uncprev <- aggregate(ViralLoad~Commune, data=unc, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
uncprev2 <- as.data.frame(cbind(uncprev$ViralLoad[,1], uncprev$ViralLoad[,2]))
colnames(uncprev2) <- c("pos","num")
uncprev2$Commune <- uncprev$Commune
uncprevdf <- merge(FLAcommcoor,uncprev2,by="Commune")
uncprevdf$cluster <- rownames(uncprevdf)
uncprevdf$LON <- as.numeric(as.character(uncprevdf$LON))
uncprevdf$LAT <- as.numeric(as.character(uncprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
unc_prevr <- as.prevR(uncprevdf, col, FLA)
unc_prevr_R <- rings(unc_prevr, N = c(68,100,200), progression = FALSE)
unc.N100 <- kde(unc_prevr_R, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
uncres <- as.data.frame(unc.N100)
uncres <- uncres[!is.na(uncres$k.prev.N100.RInf), ]
unckdeplot <- ggplot(data = uncres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='Apis Uncultured virus',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
unckdeplot
#ggsave("Plots/unckde.pdf",unckdeplot, dpi=300)


vdvprev <- aggregate(ViralLoad~Commune, data=vdv, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
vdvprev2 <- as.data.frame(cbind(vdvprev$ViralLoad[,1], vdvprev$ViralLoad[,2]))
colnames(vdvprev2) <- c("pos","num")
vdvprev2$Commune <- vdvprev$Commune
vdvprevdf <- merge(FLAcommcoor,vdvprev2,by="Commune")
vdvprevdf$cluster <- rownames(vdvprevdf)
vdvprevdf$LON <- as.numeric(as.character(vdvprevdf$LON))
vdvprevdf$LAT <- as.numeric(as.character(vdvprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
vdv_prevr <- as.prevR(vdvprevdf, col, FLA)
vdv_prevr_R <- rings(vdv_prevr, N = c(68,100,200), progression = FALSE)
vdv.N100 <- kde(vdv_prevr_R, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
vdvres <- as.data.frame(vdv.N100)
vdvres <- vdvres[!is.na(vdvres$k.prev.N100.RInf), ]
vdvkdeplot <- ggplot(data = vdvres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='VDV',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
vdvkdeplot
#ggsave("Plots/vdvkde.pdf",vdvkdeplot, dpi=300)


dwvprev <- aggregate(ViralLoad~Commune, data=dwv, FUN=function(x) c(pos=sum(x > 1000), count=length(x)))
dwvprev2 <- as.data.frame(cbind(dwvprev$ViralLoad[,1], dwvprev$ViralLoad[,2]))
colnames(dwvprev2) <- c("pos","num")
dwvprev2$Commune <- dwvprev$Commune
dwvprevdf <- merge(FLAcommcoor,dwvprev2,by="Commune")
dwvprevdf$cluster <- rownames(dwvprevdf)
dwvprevdf$LON <- as.numeric(as.character(dwvprevdf$LON))
dwvprevdf$LAT <- as.numeric(as.character(dwvprevdf$LAT))
col <- c(id = "cluster", x = "LON", y = "LAT", n = "num", pos = "pos")
dwv_prevr <- as.prevR(dwvprevdf, col, FLA)
dwv_prevr_R <- rings(dwv_prevr, N = c(68,100,200), progression = FALSE)
dwv.N100 <- kde(dwv_prevr_R, N = 100, nb.cells = 800, progression = FALSE, weighted=FALSE)
dwvres <- as.data.frame(dwv.N100)
dwvres <- dwvres[!is.na(dwvres$k.prev.N100.RInf), ]
dwvkdeplot <- ggplot(data = dwvres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='DWV',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_path(size=0.05,inherit.aes=FALSE, data=FLAgg.df,aes(long,lat,group=group),color="black") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
dwvkdeplot





#ggsave("Plots/dwvkde.pdf",dwvkdeplot, dpi=300)

#Quickly grab a legend:
legendplot <- ggplot(data = dwvres) +
  aes(x = x, y = y, fill = k.prev.N100.RInf) +
  geom_raster() +
  scale_fill_gradientn(colours=prevR.colors.green(101), limits=c(0,100)) +
  coord_fixed() +
  theme_prevR_light() +
  labs(title='DWV',fill = "Prevalence (%)") +
  theme(plot.title = element_text(hjust = 0.5))
legend <- cowplot::get_legend(legendplot)
plot(legend)


lay <- rbind(c(1,1,2,2,NA),
             c(3,3,4,4,9),
             c(5,5,6,6,9),
             c(7,7,8,8,NA))

testgrid <- grid.arrange(FLAsamples,AMFVkdeplot,BMLVkdeplot,thogkdeplot,thikakdeplot,unckdeplot,vdvkdeplot,dwvkdeplot,legend,layout_matrix =lay)
testgrid
ggsave("Plots/SpatialKDE_samplesize.png",testgrid, device = "png", dpi=300, limitsize=FALSE)
