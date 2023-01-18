library(ggtree)
library(phangorn)
library(ggplot2)
library(gridExtra)


# Apis Orthomyxovirus 1 ---------------------------------------------------

#Orthomyxo
orthotree <- read.tree("data/trees/orthomyxo.treefile")
#midpoint root
orthotree <- phangorn::midpoint(orthotree)
orthotree <- reorder(orthotree)
cls_thog <- list(
  'NODE' = c('Apisorthomyxovirus1'),
  'thogotolike' = c('AED98371.1','AHB34055.1','AHB34061.1','BAQ22394.1','QBQ64972.1','QCF29600.1','QFR36189.1','YP_009352882.1','YP_009553280.1','YP_009553280.1','YP_145794.1'),
  'unclas' = c('APG77896.1','APG77906.1','APP91612.1','QGA69818.1'),
  'influenza' = c('YP_009449556.1','NP_056657.1','AVL84663.1','AWK48712.1'))
orthotree <- groupOTU(orthotree, cls_thog)
cols_ortho <- c("unclas" = "#7f7f7f","thogotolike" = "#1f77b4", "influenza" = "#ff7f0e","NODE" = "#d62728", "NA" = "green")

orthodraw <- ggtree(orthotree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=1) +
  geom_tiplab(align=TRUE, size=3) +
  scale_colour_manual(values = cols_ortho) + 
  geom_treescale(offset=-0.5) +
  geom_cladelabel(node=23, label="Thogoto-like viruses", align=T, geom='label', col='#1f77b4', offset=0.8) +
  geom_cladelabel(node=35, label="Influenza", align=T, geom='label', col='#ff7f0e', offset=0.8) +
  xlim(0, 5) +
  ggtitle("Apis Orthomyxovirus 1") + 
  theme(plot.title = element_text(size=22)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 7))) +
  theme(legend.position = 'none')
orthodraw


# Aphtili virus -----------------------------------------------------------


#Apthili
Apthilitree <- read.tree("data/trees/Aphtili.treefile")
Apthilitree <- phangorn::midpoint(Apthilitree)
Apthilitree <- reorder(Apthilitree)

cls_aphtili <- list(
  'NODE' = c("Aphtilivirus"),
  'Thika' = c('AKH66829.1', 'AKH66830.1', 'AKH66831.1', 'AKH66832.1', 'AKH66833.1', 'AKH66834.1','AKH66835.1', 'AKH66836.1', 'AKH66837.1', 'YP_009140561.1'),
  'Old port virus' = c('QIJ25865.1', 'QIJ25866.1', 'QIJ25867.1', 'QIJ25868.1', 'QIJ25869.1'),
  'Acyrthosiphon pisum virus' = c('ABB89048.1', 'NP_620557.1', 'QAA78863.1', 'QAA78869.1', 'QAA78871.1'),
  'crawfish virus' = c('APG77989.1', 'APG78044.1', 'YP_009336630.1'),
  'Unclassified' = c('APG77457.1', 'APG78030.1', 'ARU76991.1', 'AVD69111.1', 'AWK77859.1', 'AWS06670.1','AXV43882.1', 'AYQ66683.1', 'QAY29242.1', 'QBL75890.1', 'QCI31816.1', 'QGL51720.1','QGL51724.1', 'QIJ25870.1', 'YP_009140560.1', 'YP_009336543.1', 'YP_009337038.1','YP_009337062.1', 'YP_009337100.1', 'YP_009337160.1', 'YP_009337269.1','YP_009337395.1', 'YP_009337405.1', 'YP_009337411.1')
)
cols_aphtili <- c("Unclassified" = "#7f7f7f","crawfish virus" = "#1f77b4", "Acyrthosiphon pisum virus" = "#ff7f0e","NODE" = "#d62728", "Old port virus" = "#2ca02c","Thika"="#9467bd")

Apthilitree <- groupOTU(Apthilitree, cls_aphtili)

aphtilidraw <- ggtree(Apthilitree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=1) +
  geom_tiplab(align=TRUE, size=3) +
  scale_colour_manual(values = cols_aphtili) + 
  geom_treescale(offset=-1.5) +
  geom_cladelabel(node=67, label="Thika virus", align=T, geom='label', col='#9467bd', offset=0.55) +
  geom_cladelabel(node=90, label="Acyrthosiphon pisum virus", align=T, geom='label', col='#ff7f0e', offset=0.50) +
  geom_cladelabel(node=89, label="Changjiang crawfish virus", align=T, geom='label', col='#1f77b4', offset=0.55) +
  geom_cladelabel(node=59, label="Old port virus", align=T, geom='label', col='#2ca02c', offset=0.45) +
  xlim(0, 4) +
  ggtitle("Aphtilivirus") + 
  theme(plot.title = element_text(size=22)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none')
aphtilidraw


# Apparli virus -----------------------------------------------------------

Apparlitree <- read.tree("data/trees/Apparli.treefile")
Apparlitree <- phangorn::midpoint(Apparlitree)
Apparlitree <- reorder(Apparlitree)
cls_apparli <- list(
  "NODE" = c("Apparlivirus"),
  #"RNA_virosphere" = c('APG78155.1',  'APG78162.1',  'APG78183.1',  'APG78199.1',  'APG78200.1',  'APG78203.1',  'APG78213.1',  'APG78215.1',  'APG78216.1',  'APG78217.1',  'APG78218.1',  'APG78222.1',  'APG78224.1',  'APG78227.1',  'APG78231.1',  'APG78233.1',  'APG78238.1',  'APG78243.1',  'APG78244.1',  'APG78245.1',  'APG78248.1',  'APG78249.1',  'APG78251.1',  'APG78252.1',  'APG78253.1',  'APG78254.1',  'APG78256.1',  'APG78257.1',  'APG78260.1',  'APG78261.1',  'APG78266.1',  'APG78275.1',  'APG78276.1',  'APG78277.1',  'APG78278.1',  'APG78281.1',  'APG78282.1',  'APG78283.1',  'APG78293.1',  'APG78294.1',  'APG78296.1',  'APG78301.1',  'APG78306.1',  'APG78307.1',  'APG78310.1',  'APG78313.1',  'APG78316.1',  'APG78317.1',  'APG78322.1',  'APG78330.1',  'APG78331.1',  'APG78332.1',  'APG78333.1',  'APG78342.1',  'APG78345.1',  'APG78364.1',  'BBQ05104.1',  'BBQ05105.1',  'YP_009329869.1',  'YP_009329875.1',  'YP_009329882.1',  'YP_009333350.1',  'YP_009333370.1',  'YP_009337885.1',  'YP_009342308.1',  'YP_009345133.1'),
  #"Insect_derived" = c('AEA40853.1',  'AGW51759.1',  'ASV45859.1',  'ATY36109.1',  'AXQ04877.1',  'BBE15516.1',  'KFB47171.1',  'QBP37031.1',  'QFR15908.1',  'QGA70940.1',  'QGL51732.1',  'QHA33703.1',  'QHA33832.1',  'QHA33899.1',  'QHA33901.1',  'QHA33902.1'),
  "RNAvirosph_insect" = c('APG78155.1',  'APG78162.1',  'APG78183.1',  'APG78199.1',  'APG78200.1',  'APG78203.1',  'APG78213.1',  'APG78215.1',  'APG78216.1',  'APG78217.1',  'APG78218.1',  'APG78222.1',  'APG78224.1',  'APG78227.1',  'APG78231.1',  'APG78233.1',  'APG78238.1',  'APG78243.1',  'APG78244.1',  'APG78245.1',  'APG78248.1',  'APG78249.1',  'APG78251.1',  'APG78252.1',  'APG78253.1',  'APG78254.1',  'APG78256.1',  'APG78257.1',  'APG78260.1',  'APG78261.1',  'APG78266.1',  'APG78275.1',  'APG78276.1',  'APG78277.1',  'APG78278.1',  'APG78281.1',  'APG78282.1',  'APG78283.1',  'APG78293.1',  'APG78294.1',  'APG78296.1',  'APG78301.1',  'APG78306.1',  'APG78307.1',  'APG78310.1',  'APG78313.1',  'APG78316.1',  'APG78317.1',  'APG78322.1',  'APG78330.1',  'APG78331.1',  'APG78332.1',  'APG78333.1',  'APG78342.1',  'APG78345.1',  'APG78364.1',  'BBQ05104.1',  'BBQ05105.1',  'YP_009329869.1',  'YP_009329875.1',  'YP_009329882.1',  'YP_009333350.1',  'YP_009333370.1',  'YP_009337885.1',  'YP_009342308.1',  'YP_009345133.1', 'AEA40853.1',  'AGW51759.1',  'ASV45859.1',  'ATY36109.1',  'AXQ04877.1',  'BBE15516.1',  'KFB47171.1',  'QBP37031.1',  'QFR15908.1',  'QGA70940.1',  'QGL51732.1',  'QHA33703.1',  'QHA33832.1',  'QHA33899.1',  'QHA33901.1',  'QHA33902.1'),
  "Partiti" = c('AOR51388.1', 'AWD38959.1', 'QGZ98411.1'),
  "Animal" = c('QIJ70076.1', 'QIJ70086.1', 'QIJ70087.1', 'QIJ70088.1', 'QIJ70089.1')
)
cols_apparli <-c("Animal" = "#1f77b4", "Partiti" = "#ff7f0e","NODE" = "#d62728","RNAvirosph_insect"="#9467bd")

Apparlitree <- groupOTU(Apparlitree, cls_apparli)

apparlidraw <- ggtree(Apparlitree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=1) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_apparli) + 
  geom_treescale(offset=-1.5) +
  xlim(0, 3) +
  theme(legend.position="right") +
  ggtitle("Apparli virus") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=22)) +
  scale_color_manual(name="Viral group",values=c("#1f77b4","#d62728","#ff7f0e","#9467bd"), labels=c("Animal-derived Partitiviruses", "Apis Uncultured virus","Partiti(-like) viruses","RNA Virosphere partiti-like viruses"))
apparlidraw


# BMLV --------------------------------------------------------------------
bmlvtree <- read.tree("data/trees/BMLV.treefile")
bmlvtree <- phangorn::midpoint(bmlvtree)
bmlvtree <- reorder(bmlvtree)
cls_bmlv <- list(
  "NODE" = c("Bee-macula-likevirus_Belgium"),
  "Tymovirus" = c('AAB92649.1',  'AFC95826.1',  'AIE44664.1',  'AIE44667.1',  'AIG05827.1',  'ALA65431.1',  'ALA65434.1',  'AMC38501.1',  'AMH40122.1',  'AMH40125.1',  'AMH40128.1',  'AMH40131.1',  'AMH40134.1',  'AMH40137.1',  'AMH40140.1',  'ARK20051.1',  'AVW89220.1',  'AYC35256.1',  'CBH31042.1',  'NP_040968.1',  'NP_041257.1',  'NP_047920.1',  'NP_619756.1',  'NP_663297.1',  'P10358.2',  'P20128.1',  'P28477.1',  'YP_001285472.1',  'YP_002048673.1',  'YP_002308439.1',  'YP_002308442.1',  'YP_002308445.1',  'YP_002308578.1',  'YP_003620401.1',  'YP_004464924.1',  'YP_007517180.1',  'YP_007517183.1',  'YP_008318042.1'),
  "Marafivirus" = c('AAK52838.2',  'AAW33732.1',  'ABA54133.1',  'ACV83739.1',  'AFV34757.1',  'AIY22517.1',  'AKQ08147.1',  'AKQ08148.1',  'AKQ08149.1',  'AKZ17743.1',  'AKZ17760.1',  'ALX72769.1',  'ALX72770.1',  'ANP94321.1',  'AOX49247.1',  'AQV11970.1',  'ARH56444.1',  'ARO38248.1',  'ARO38249.1',  'ATB20100.1',  'ATB20101.1',  'ATJ00054.1',  'AXM42964.1',  'AXM42965.1',  'AYC35261.1',  'NP_115454.1',  'QBZ78634.1',  'QCC30252.1',  'QCC30253.1',  'QID59002.1',  'YP_002756536.1',  'YP_004464920.1',  'YP_009222597.1',  'YP_009315883.1',  'YP_009345914.1',  'YP_009351862.1',  'YP_009505639.1',  'YP_009551972.1',  'YP_224218.1'),
  "Unclassified" = c('AEP40395.1',  'AKQ48574.1',  'ARI47198.1',  'ARI47200.1',  'AUR53413.1',  'AYV61007.1',  'AZF99024.1',  'BAJ14665.1',  'BAJ14667.1',  'BAJ14669.1',  'BAJ14673.1',  'CAL80776.1',  'NP_037647.1',  'QBP37021.1',  'YP_009159826.1',  'YP_009160324.1',  'YP_009505642.1'),
  "Maculavirus" = c('AHX22590.1',  'ANV22069.1',  'ANV22095.1',  'AOX49245.1',  'YP_004464930.1',  'YP_009268923.1')
)

cols_bmlv <-c("Tymovirus" = "#1f77b4","Unclassified"= "#7f7f7f", "Marafivirus" = "#ff7f0e","NODE" = "#d62728","Maculavirus"="#9467bd")

bmlvtree <- groupOTU(bmlvtree, cls_bmlv)

bmlvdraw <- ggtree(bmlvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=1) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_bmlv) + 
  geom_cladelabel(node=149, label="Marafivirus", align=T, geom='label', col='#ff7f0e', offset=0.35) +
  geom_cladelabel(node=154, label="Tymovirus", align=T, geom='label', col='#1f77b4', offset=0.35) +
  geom_cladelabel(node=110, label="Marafivirus", align=T, geom='label', col='#ff7f0e', offset=0.35) +
  geom_cladelabel(node=150, label="Maculavirus", align=T, geom='label', col='#9467bd', offset=0.35) +
  geom_cladelabel(node=198, label="Maculavirus", align=T, geom='label', col='#9467bd', offset=0.35) +
  geom_treescale(offset=-1.5) +
  ggtitle("Bee-macula like virus") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=22)) +
  xlim(c(0,2.5)) +
  theme(legend.position = 'none')
bmlvdraw


# AMFV --------------------------------------------------------------------
amfvtree <- read.tree("data/trees/AMFV.treefile")
amfvtree <- phangorn::midpoint(amfvtree)
amfvtree <- reorder(amfvtree)
cls_amfv <- list(
  "NODE" = c("ApisMeliferaFilamentousvirus_Belgium"),
  "Poxviridae" = c('AAA60776.1',  'AAA69333.1',  'AAA69374.1',  'AAF33894.1',  'AAK43560.1',  'AAN02745.1',  'AAS49745.1',  'AAW23440.1',  'ABD97395.1',  'ABH08147.1',  'ABI99013.1',  'ABQ43492.1',  'ACB28637.1',  'ADZ29378.1',  'ADZ29593.1',  'ADZ29806.1',  'ADZ30447.1',  'ADZ30659.1',  'AEY72851.1',  'AFU76948.1',  'AGJ91195.1',  'AGJ91739.1',  'AGU99698.1',  'AGY97303.1',  'AGY97736.1',  'AGZ95335.1',  'AMB18349.1',  'ARB50277.1',  'ARR29609.2',  'ARR30008.1',  'ARR30425.1',  'ARR30628.1',  'ARR31017.1',  'AUI80588.1',  'AXN74832.1',  'AXN75271.1',  'AYP74130.1',  'NP_051729.1',  'NP_051904.1',  'NP_073405.1',  'NP_150454.1',  'NP_570429.1',  'NP_619839.1',  'NP_659593.1',  'NP_938276.1',  'O57175.1',  'P0DSS8.1',  'P20493.1',  'QAV34328.1',  'QAV42102.1',  'QDJ94996.1',  'QED21168.1',  'QEJ79019.1',  'QEJ79168.1',  'QEJ79318.1',  'SNB49884.1',  'SNB53815.1',  'YP_001293211.1',  'YP_004821373.1',  'YP_005296225.1',  'YP_009143350.1',  'YP_009268721.1',  'YP_009281788.1',  'YP_009282734.1',  'YP_009329641.1',  'YP_009407970.1',  'YP_009408218.1',  'YP_009408420.1',  'YP_227403.1',  'YP_232925.1',  'YP_717351.1'),
  "Unclassified" = c('ADD74390.1', 'YP_009165967.1'),
  "Phycodnaviridae" = c('AEO97836.1', 'AEO98243.1', 'CAZ69365.1', 'YP_293780.1'),
  "Baculoviridae" = c('AHC69623.1',  'AIX47958.1',  'AJR20393.1',  'AMO27616.1',  'AMO27788.1',  'AMO27974.1',  'ANS71004.1',  'AOW42794.1',  'AQQ80135.1',  'AWJ76733.1',  'NP_047757.1',  'NP_258331.1',  'QCQ67370.1',  'QCQ67530.1',  'QCQ67689.1',  'QDH05963.1',  'YP_003517877.1',  'YP_009505866.1'),
  "Other_families" = c('AUS94194.1','AXN91020.1','BAV31392.1','QCQ67809.1','YP_001883390.1')
)

cols_amfv <-c("Poxviridae" = "#1f77b4","Unclassified"= "#7f7f7f", "Phycodnaviridae" = "#ff7f0e","NODE" = "#d62728","Baculoviridae"="#9467bd","Other_families"="#7f7f7f", "0"="#7f7f7f")
amfvtree <- groupOTU(amfvtree, cls_amfv)

amfvdraw <- ggtree(amfvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=1) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_amfv) + 
  geom_treescale(offset=-1.5) +
  geom_cladelabel(node=108, label="Poxviridae", align=T, geom='label', col='#1f77b4', offset=0.24) +
  geom_cladelabel(node=180, label="Baculoviridae", align=T, geom='label', col='#9467bd', offset=0.24) +
  geom_cladelabel(node=198, label="Phycodnaviridae", align=T, geom='label', col='#ff7f0e', offset=0.19) +
  ggtitle("Apis Melifera filamentous virus") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=22)) +
  xlim(0, 1.6 ) +
  theme(legend.position = 'none')
amfvdraw


# dwv - vdv ---------------------------------------------------------------
vdvtree <- read.tree("data/trees/vdv_dwv.treefile")
#midpoint root
vdvtree <- phangorn::midpoint(vdvtree)
vdvtree <- reorder(vdvtree)
cls_vdv <- list('NODE'=c("DeformedWingVirus_Belgium","VarroaDestructorVirus_Belgium"),
                'Iflaviridae' = c('YP_003622540.1',
                                  'YP_008130309.1',
                                  'YP_008130310.1',
                                  'YP_008888537.1',
                                  'YP_009002581.1',
                                  'YP_009010941.1',
                                  'YP_009026409.1',
                                  'YP_009047245.1',
                                  'YP_009116875.1',
                                  'YP_009129265.1',
                                  'YP_009162630.1',
                                  'YP_009305421.1',
                                  'YP_009315906.1',
                                  'YP_009328891.1',
                                  'YP_009344960.1',
                                  'YP_009345906.1',
                                  'YP_009351892.1',
                                  'YP_009444707.1',
                                  'YP_009505599.1',
                                  'YP_009552080.1',
                                  'YP_009552119.1',
                                  'YP_009553259.1',
                                  'YP_145791.1',
                                  'NP_853560.2'),
                'Unclassified' = c('YP_009329861.1',
                                   'YP_009330050.1',
                                   'YP_009330055.1',
                                   'YP_009333199.1',
                                   'YP_009336533.1',
                                   'YP_009336552.1',
                                   'YP_009336575.1',
                                   'YP_009336939.1',
                                   'YP_009337003.1',
                                   'YP_009337046.1',
                                   'YP_009337127.1',
                                   'YP_009337161.1',
                                   'YP_009337271.1',
                                   'YP_009337284.1',
                                   'YP_009337760.1',
                                   'YP_009342053.1',
                                   'YP_009342321.1'))
cols_vdv <-c("Iflaviridae" = "#1f77b4","Unclassified"= "#7f7f7f", "NODE" = "#d62728")
vdvtree <- groupOTU(vdvtree, cls_vdv)

vdvdraw <- ggtree(vdvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=1) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_vdv)+
  geom_treescale(offset=-1.5) +
  ggtitle("DWV/VDV") + 
  xlim(0, 4)+
  theme(legend.position=c(0.8,0.5)) +
  theme(plot.title = element_text(size=22)) +
  ggtitle("VDV/DWV") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Viral group",values=c("#1f77b4","#d62728","#7f7f7f"), labels=c("Iflaviridae","Iflavirus Belgium","Unclassified"))
vdvdraw






lay <- rbind(c(1,1,2,2,3,3),
             c(1,1,2,2,3,3),
             c(4,4,5,5,6,6))



testgrid <- grid.arrange(amfvdraw,bmlvdraw,apparlidraw,vdvdraw, aphtilidraw ,orthodraw,layout_matrix =lay)

ggsave("output/fig1.png",testgrid, device = "png", height=14,width=23.5, limitsize=FALSE)
