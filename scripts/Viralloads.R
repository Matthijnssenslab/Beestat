virdf <- read.csv("../Data/Virdf_metadata.csv", row.names=1)
virdf$Year <- as.factor(as.character(virdf$Year))
cumulcount <- read.csv("../Data/Virdf_metadata_cumulcount.csv")
cumulcount$X <- NULL
cumulcount$Year <- as.factor(as.character(cumulcount$Year))
virdf <- rbind(virdf, cumulcount)

#Check viral loads + MWU-tests
library(ggplot2)
library(ggpubr)
cols <- c("Virus","Year","ViralLoad","Status")
viralloads <- (virdf[,colnames(virdf) %in% cols])
viralloads$logload <- log10(viralloads$ViralLoad + 1)
levels(viralloads$Virus) <- c('Apis Melifera filamentous virus','Bee-macula like virus','DWV','Apis Thika-like virus','Apis Orthomyxovirus 1','Apis Uncultured virus','VDV','Cumulative')
levels(viralloads$Status) <- c('Healthy','Weak')
levels(viralloads$Year) <- c('2012','2013')
compare_means(ViralLoad ~ Virus, data=viralloads, group.by="Status")
#Get actual numbers:
#DWV:
dwvquick = viralloads[viralloads$Virus == 'DWV',]
vdvquick = viralloads[viralloads$Virus == 'VDV',]
thikaquick = viralloads[viralloads$Virus == 'Apis Thika-like virus',]
amfvquick = viralloads[viralloads$Virus == 'Apis Melifera filamentous virus',]
orthoquick = viralloads[viralloads$Virus == 'Apis Orthomyxovirus 1',]
bmlvquick = viralloads[viralloads$Virus == 'Bee-macula like virus',]

wilcox.test(dwvquick[dwvquick$Status == 'Healthy',]$ViralLoad, dwvquick[dwvquick$Status == 'Weak',]$ViralLoad,)
wilcox.test(vdvquick[vdvquick$Status == 'Healthy',]$ViralLoad, vdvquick[vdvquick$Status == 'Weak',]$ViralLoad,)
wilcox.test(thikaquick[thikaquick$Status == 'Healthy',]$ViralLoad, thikaquick[thikaquick$Status == 'Weak',]$ViralLoad,)
wilcox.test(amfvquick[amfvquick$Status == 'Healthy',]$ViralLoad, amfvquick[amfvquick$Status == 'Weak',]$ViralLoad,)
wilcox.test(orthoquick[orthoquick$Status == 'Healthy',]$ViralLoad, orthoquick[orthoquick$Status == 'Weak',]$ViralLoad,)


logviralload_stat <- ggplot(data=viralloads, aes(x=reorder(viralloads$Virus, viralloads$ViralLoad),y=viralloads$logload,fill=viralloads$Status)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha=0.2) +
  geom_violin(color=NA, alpha=0.4) +
  stat_compare_means(method = "wilcox.test",label.y=10.5, aes(x=Virus, y=ViralLoad,group = Status), label = c("p.signif"),hide.ns=TRUE) +
  theme_minimal() +
  labs(x="Virus", y="Log10(Viralload + 1", fill="Status") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2)) + 
  scale_fill_manual(values=c("#2ca02c", "#d62728"))
logviralload_stat
ggsave("../Plots/logviralload_stat.pdf",logviralload_stat,dpi=300)

#yearly:
logviralload_year <- ggplot(data=viralloads, aes(x=reorder(viralloads$Virus, viralloads$ViralLoad),y=viralloads$logload,fill=viralloads$Year)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha=0.2) +
  geom_violin(color=NA, alpha=0.4) +
  stat_compare_means(method = "wilcox.test",label.y=10.5, aes(x=Virus, y=ViralLoad,group = Year), label = c("p.signif"),hide.ns=TRUE) +
  theme_minimal() +
  labs(x="Virus", y="Log10(Viralload + 1", fill="Year") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.2)) + 
  scale_fill_manual(values=c("#2ca02c", "#d62728"))
logviralload_year
ggsave("../Plots/logviralload_year.pdf",logviralload_year,dpi=300)

#Actual stats for temporal effects:
wilcox.test(dwvquick[dwvquick$Year == '2012',]$ViralLoad, dwvquick[dwvquick$Year == '2013',]$ViralLoad,)
wilcox.test(vdvquick[vdvquick$Year == '2012',]$ViralLoad, vdvquick[vdvquick$Year == '2013',]$ViralLoad,)
wilcox.test(thikaquick[thikaquick$Year == '2012',]$ViralLoad, thikaquick[thikaquick$Year == '2013',]$ViralLoad,)
wilcox.test(bmlvquick[bmlvquick$Year == '2012',]$ViralLoad, bmlvquick[bmlvquick$Year == '2013',]$ViralLoad,)

