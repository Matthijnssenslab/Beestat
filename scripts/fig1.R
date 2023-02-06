library(ggtree)
library(phangorn)
library(ggplot2)
library(gridExtra)

# Apis Orthomyxovirus 1 ---------------------------------------------------

# acc 2 name mappings
aov.acc2name <- read.delim("data/trees/orthomyxo.acc2name", header=FALSE)
aov.a2n <- c(paste(aov.acc2name$V1, aov.acc2name$V2, sep=' - '))
names(aov.a2n) <- aov.acc2name$V1

#Orthomyxo
orthotree <- read.tree("data/trees/orthomyxo.treefile")
orthotree$tip.label <- as.character(aov.a2n[match(names(aov.a2n), orthotree$tip.label)])
orthotree$tip.label[1] <- "Apis orthomyxovirus 1"

#midpoint root
orthotree <- phangorn::midpoint(orthotree)
orthotree <- reorder(orthotree)
cls_thog <- list(
  'NODE' = c(
    'Apis orthomyxovirus 1'
  ),
  'thogotolike' = c(
    "AED98371.1 - Jos virus","AHB34055.1 - Upolu virus",
    "AHB34061.1 - Aransas Bay virus","BAQ22394.1 - Thogotovirus thogotoense",
    "YP_145794.1 - Thogotovirus thogotoense","QBQ64972.1 - Thogotovirus dhoriense",
    "YP_009352882.1 - Thogotovirus dhoriense","QCF29600.1 - Bourbon virus",
    "YP_009553280.1 - Oz virus","QFR36189.1 - Thailand tick thogotovirus"
  ),
  'unclas' = c(
    "APG77906.1 - Hubei orthoptera virus 6",
    "APG77896.1 - Hubei orthomyxo-like virus 4",
    "QGA69818.1 - Varroa orthomyxovirus-1",
    "APP91612.1 - Sinu virus"
  ),
  'influenza' = c(
    'AWK48712.1 - Influenza C virus',
    'AVL84663.1 - Influenza A virus',
    'YP_009449556.1 - Influenza D virus (D/swine/Oklahoma/1334/2011)',
    "NP_056657.1 - Influenza B virus (B/Lee/1940)"
  ))
orthotree <- groupOTU(orthotree, cls_thog)
cols_ortho <- c("unclas" = "#35274A","thogotolike" = "#0B775E", "influenza" = "#EABE94","NODE" = "#F2300F")

orthodraw <- ggtree(orthotree, aes(color=group)) + 
  geom_tiplab(align = TRUE, aes(subset = isTip & label %in% c(
    'Apis orthomyxovirus 1',
    "QGA69818.1 - Varroa orthomyxovirus-1",
    "APG77896.1 - Hubei orthomyxo-like virus 4",
    "APG77906.1 - Hubei orthoptera virus 6",
    "APP91612.1 - Sinu virus"
  ))) +
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  #geom_tiplab(align=TRUE, size=3) +
  scale_colour_manual(values = cols_ortho) + 
  geom_treescale(2,5) +
  geom_cladelabel(node=23, label="Thogoto-like viruses", align=T, geom='label', col='#0B775E') +
  geom_cladelabel(node=35, label="Influenza", align=T, geom='label', col='#EABE94') +
  ggtitle("Apis Orthomyxovirus 1") + 
  theme(plot.title = element_text(size=22)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(0,5)) +
  theme(legend.position = 'none')
orthodraw <- scaleClade(orthodraw, node=23, scale=.1) 
orthodraw <- scaleClade(orthodraw, node=35, scale=.1) 

orthodraw_full <- ggtree(orthotree, aes(color=group)) + 
  geom_tiplab(align=TRUE) +
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  scale_colour_manual(values = cols_ortho) + 
  geom_treescale(offset=-0.5) +
  geom_cladelabel(node=23, label="Thogoto-like viruses", align=T, geom='label', col='#0B775E', offset=1.4) +
  geom_cladelabel(node=35, label="Influenza", align=T, geom='label', col='#EABE94', offset=2.2) +
  ggtitle("Apis Orthomyxovirus 1") + 
  theme(plot.title = element_text(size=22)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(0,5)) +
  theme(legend.position = 'none')
ggsave('output/sup_fig7.png', orthodraw_full, dpi=300, width=12, height=12)

# apthili virus -----------------------------------------------------------

# acc 2 name mappings
apthili.acc2name <- read.delim("data/trees/apthili.acc2name", header=FALSE)
apthili.acc2name$V3 <- paste(apthili.acc2name$V1, apthili.acc2name$V2, sep=' - ')
apthili.a2n <- c(paste(apthili.acc2name$V1, apthili.acc2name$V2, sep=' - '))
names(apthili.a2n) <- apthili.acc2name$V1

#apthili
apthilitree <- read.tree("data/trees/apthili.treefile")
apthilitree$tip.label <- as.character(apthili.a2n[match(names(apthili.a2n), apthilitree$tip.label)])
apthilitree$tip.label[1] <- "Apthili virus"
apthilitree <- phangorn::midpoint(apthilitree)
apthilitree <- reorder(apthilitree)
cls_aphtili <- list(
  'NODE' = c("Apthili virus"),
  'Thika' = c(
    "YP_009140561.1 - Thika virus",
    "AKH66831.1 - Thika virus",
    "AKH66832.1 - Thika virus",
    "AKH66829.1 - Thika virus",
    "AKH66830.1 - Thika virus",
    "AKH66835.1 - Thika virus",
    "AKH66836.1 - Thika virus",
    "AKH66837.1 - Thika virus",
    "AKH66834.1 - Thika virus",
    "AKH66833.1 - Thika virus"
  ),
  'Old port virus' = c(
    "QIJ25868.1 - Old Port virus",
    "QIJ25867.1 - Old Port virus",
    "QIJ25869.1 - Old Port virus",
    "QIJ25865.1 - Old Port virus",
    "QIJ25866.1 - Old Port virus"
  ),
  'Acyrthosiphon pisum virus' = c(
    "NP_620557.1 - Acyrthosiphon pisum virus",
    "QAA78863.1 - Acyrthosiphon pisum virus",
    "QAA78871.1 - Acyrthosiphon pisum virus",
    "QAA78869.1 - Acyrthosiphon pisum virus",
    "ABB89048.1 - Rosy apple aphid virus"
  ),
  'crawfish virus' = c(
    "APG77989.1 - Changjiang crawfish virus 6",
    "APG78044.1 - Changjiang crawfish virus 6",
    "YP_009336630.1 - Changjiang crawfish virus 6"
  ),
  'Unclassified' = c(
    "YP_009337411.1 - Hubei picorna-like virus 58",
    "APG77457.1 - Hubei picorna-like virus 60",
    "YP_009337395.1 - Hubei picorna-like virus 60",
    "AXV43882.1 - Yongsan picorna-like virus 2",
    "YP_009337405.1 - Hubei tetragnatha maxillosa virus 5",
    "YP_009337038.1 - Hubei picorna-like virus 63",
    "YP_009337160.1 - Hubei picorna-like virus 62",
    "QIJ25870.1 - Warroolaba Creek virus 3",
    "YP_009337062.1 - Hubei picorna-like virus 61",
    "YP_009337100.1 - Hubei picorna-like virus 59",
    "YP_009140560.1 - Kilifi Virus",
    "AYQ66683.1 - Kilifi Virus",
    "QBL75890.1 - Solenopsis invicta virus 7",
    "AWK77859.1 - Bundaberg bee virus 8",
    "ARU76991.1 - Milolii virus",
    "AVD69111.1 - HVAC-associated RNA virus 1",
    "QGL51724.1 - Vespa velutina associated acypi-like virus",
    "QGL51720.1 - Vespa velutina associated acypi-like virus",
    "APG78030.1 - Hubei picorna-like virus 57",
    "QAY29242.1 - Boghill Burn virus",
    "YP_009336543.1 - Hubei picorna-like virus 56",
    "YP_009337269.1 - Hubei picorna-like virus 55",
    "AWS06670.1 - Hubei picorna-like virus 55",
    "QCI31816.1 - Sitobion miscanthi virus 1"
  )
)
cols_aphtili <- c(
  "Unclassified" = "#35274A",
  "crawfish virus" = "#0B775E",
  "Acyrthosiphon pisum virus" = "#EABE94",
  "NODE" = "#F2300F",
  "Old port virus" = "#7FC0C6",
    "Thika"="#02401B"
)
apthilitree <- groupOTU(apthilitree, cls_aphtili)
apthilidraw <- ggtree(apthilitree, aes(color=group)) + 
  geom_tiplab(align=TRUE,aes(subset = isTip & label %in% c(
    "YP_009337411.1 - Hubei picorna-like virus 58",
    "APG77457.1 - Hubei picorna-like virus 60",
    "YP_009337395.1 - Hubei picorna-like virus 60",
    "AXV43882.1 - Yongsan picorna-like virus 2",
    "YP_009337405.1 - Hubei tetragnatha maxillosa virus 5",
    "YP_009337038.1 - Hubei picorna-like virus 63",
    "YP_009337160.1 - Hubei picorna-like virus 62",
    "QIJ25870.1 - Warroolaba Creek virus 3",
    "YP_009337062.1 - Hubei picorna-like virus 61",
    "YP_009337100.1 - Hubei picorna-like virus 59",
    "YP_009140560.1 - Kilifi Virus",
    "AYQ66683.1 - Kilifi Virus",
    "QBL75890.1 - Solenopsis invicta virus 7",
    "AWK77859.1 - Bundaberg bee virus 8",
    "ARU76991.1 - Milolii virus",
    "AVD69111.1 - HVAC-associated RNA virus 1",
    "QGL51724.1 - Vespa velutina associated acypi-like virus",
    "QGL51720.1 - Vespa velutina associated acypi-like virus",
    "APG78030.1 - Hubei picorna-like virus 57",
    "QAY29242.1 - Boghill Burn virus",
    "YP_009336543.1 - Hubei picorna-like virus 56",
    "YP_009337269.1 - Hubei picorna-like virus 55",
    "AWS06670.1 - Hubei picorna-like virus 55",
    'Apthili virus'
  ))) +
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  #geom_tiplab(align=TRUE, size=3) +
  scale_colour_manual(values = cols_aphtili) + 
  geom_treescale(2,4.5) +
  geom_cladelabel(node=67, label="Thika virus", align=T, geom='label', col='#02401B') +
  geom_cladelabel(node=90, label="Acyrthosiphon pisum virus", align=T, geom='label', col='#EABE94') +
  geom_cladelabel(node=89, label="Changjiang crawfish virus", align=T, geom='label', col='#0B775E') +
  geom_cladelabel(node=59, label="Old port virus", align=T, geom='label', col='#7FC0C6') +
  xlim(0, 6) +
  ggtitle("Apthilivirus") + 
  theme(plot.title = element_text(size=22)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none')
apthilidraw <- scaleClade(apthilidraw, node=67, scale=.1) 
apthilidraw <- scaleClade(apthilidraw, node=90, scale=.1) 
apthilidraw <- scaleClade(apthilidraw, node=89, scale=.1) 
apthilidraw <- scaleClade(apthilidraw, node=59, scale=.1) 
apthilidraw

apthilidraw_full <- ggtree(apthilitree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(align=TRUE, size=3) +
  scale_colour_manual(values = cols_aphtili) + 
  geom_treescale(offset=-1.5) +
  geom_cladelabel(node=67, label="Thika virus", align=T, geom='label', col='#02401B', offset=0.6) +
  geom_cladelabel(node=90, label="Acyrthosiphon pisum virus", align=T, geom='label', col='#EABE94', offset=0.9) +
  geom_cladelabel(node=89, label="Changjiang crawfish virus", align=T, geom='label', col='#0B775E', offset=0.9) +
  geom_cladelabel(node=59, label="Old port virus", align=T, geom='label', col='#7FC0C6', offset=0.55) +
  xlim(0, 4) +
  ggtitle("Apthili virus") + 
  theme(plot.title = element_text(size=22)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none')
apthilidraw_full
ggsave('output/sup_fig6.png', apthilidraw_full, dpi=300, width=12, height=12)

# Apparli virus -----------------------------------------------------------
# acc 2 name mappings
apparli.acc2name <- read.delim("data/trees/apparli.acc2name", header=FALSE)
apparli.acc2name$V3 <- paste(apparli.acc2name$V1, apparli.acc2name$V2, sep=' - ')
apparli.a2n <- c(paste(apparli.acc2name$V1, apparli.acc2name$V2, sep=' - '))
names(apparli.a2n) <- apparli.acc2name$V1
#apparli
apparlitree <- read.tree("data/trees/apparli.treefile")
apparlitree$tip.label <- as.character(apparli.a2n[match(names(apparli.a2n), apparlitree$tip.label)])
apparlitree$tip.label[1] <- "Apparli virus"
apparlitree <- phangorn::midpoint(apparlitree)
apparlitree <- reorder(apparlitree)

cls_apparli <- list(
  "Apparli virus Belgium" = c("Apparli virus"),
  "RNA Virosphere" = c(
    "APG78322.1 - Hubei partiti-like virus 34",
    "AEA40853.1 - uncultured virus",
    "ATY36109.1 - Vespa velutina partiti-like virus 1",
    "APG78222.1 - Hubei partiti-like virus 33",
    "BBE15516.1 - Osugoroshi virus 1",
    "QHA33899.1 - Atrato Partiti-like virus 3",
    "QHA33901.1 - Atrato Partiti-like virus 3",
    "APG78254.1 - Hubei odonate virus 13",
    "APG78233.1 - Hubei partiti-like virus 39",
    "APG78200.1 - Wuhan insect virus 25",
    "APG78199.1 - Wuhan insect virus 24",
    "APG78249.1 - Hubei partiti-like virus 36",
    "APG78333.1 - Hubei partiti-like virus 35",
    "YP_009329882.1 - Wuhan insect virus 23",
    "APG78216.1 - Wuhan insect virus 23",
    "APG78231.1 - Hubei partiti-like virus 37",
    "KFB47171.1 - Anopheles sinensis",
    "APG78261.1 - Hubei partiti-like virus 38",
    "APG78203.1 - Wuhan fly virus 6",
    "APG78332.1 - Wuhan house centipede virus 8",
    "APG78213.1 - Shuangao partiti-like virus 1",
    "APG78215.1 - Shuangao partiti-like virus 1",
    "APG78243.1 - Hubei partiti-like virus 46",
    "APG78331.1 - Wuhan house centipede virus 7",
    "APG78330.1 - Hubei partiti-like virus 45",
    "APG78313.1 - Hubei partiti-like virus 43",
    "APG78364.1 - Hubei partiti-like virus 43",
    "AGW51759.1 - uncultured virus",
    "APG78310.1 - Hubei partiti-like virus 40",
    "YP_009329869.1 - Beihai barnacle virus 13",
    "APG78238.1 - Hubei partiti-like virus 41",
    "APG78248.1 - Hubei partiti-like virus 44",
    "ASV45859.1 - Araticum virus",
    "APG78281.1 - Hubei partiti-like virus 42",
    "APG78345.1 - Wenling partiti-like virus 3",
    "APG78245.1 - Hubei diptera virus 19",
    "APG78266.1 - Hubei diptera virus 19",
    "QFR15908.1 - Spodoptera exempta insect virus 1",
    "QBP37031.1 - Lampyris noctiluca partitivirus-like virus 1",
    "APG78277.1 - Hubei partiti-like virus 31",
    "APG78244.1 - Hubei partiti-like virus 29",
    "APG78342.1 - Wenling partiti-like virus 2",
    "APG78251.1 - Hubei partiti-like virus 32",
    "APG78296.1 - Sanxia water strider virus 18",
    "APG78293.1 - Sanxia water strider virus 18",
    "YP_009333350.1 - Beihai partiti-like virus 2",
    "APG78282.1 - Hubei partiti-like virus 15",
    "APG78316.1 - Hubei partiti-like virus 16",
    "APG78301.1 - Hubei partiti-like virus 16",
    "AXQ04877.1 - Partitivirus-like Culex mosquito virus",
    "QGA70940.1 - Sonnbo virus",
    "APG78278.1 - Hubei partiti-like virus 14",
    "APG78275.1 - Hubei partiti-like virus 13",
    "YP_009329875.1 - Hubei partiti-like virus 11",
    "APG78227.1 - Hubei partiti-like virus 10",
    "APG78307.1 - Hubei partiti-like virus 8",
    "YP_009333370.1 - Beihai barnacle virus 12",
    "APG78183.1 - Beihai partiti-like virus 1",
    "YP_009342308.1 - Wuhan Millipede virus 4",
    "APG78317.1 - Wuhan Millipede virus 4",
    "YP_009345133.1 - Wuhan cricket virus 2",
    "APG78276.1 - Hubei partiti-like virus 6",
    "APG78252.1 - Hubei partiti-like virus 5",
    "APG78306.1 - Hubei partiti-like virus 1",
    "YP_009337885.1 - Hubei tetragnatha maxillosa virus 8",
    "APG78155.1 - Hubei partiti-like virus 2",
    "QHA33832.1 - Atrato Partiti-like virus 1",
    "APG78162.1 - Hubei partiti-like virus 3",
    "APG78224.1 - Hubei partiti-like virus 4",
    "QHA33902.1 - Atrato Partiti-like virus 2",
    "APG78257.1 - Hubei partiti-like virus 12",
    "APG78253.1 - Hubei partiti-like virus 20",
    "APG78260.1 - Hubei partiti-like virus 19",
    "APG78256.1 - Hubei partiti-like virus 17",
    "APG78283.1 - Hubei partiti-like virus 22",
    "BBQ05104.1 - Hubei partiti-like virus 22",
    "BBQ05105.1 - Hubei partiti-like virus 22",
    "APG78217.1 - Hubei partiti-like virus 22",
    "QGL51732.1 - Vespa velutina associated partiti-like virus 3",
    "APG78294.1 - Sanxia partiti-like virus 1",
    "APG78218.1 - Hubei partiti-like virus 48",
    "QHA33703.1 - Atrato Partiti-like virus 4"
  ),
  "Partiti (plant/fungus)" = c(
    "AWD38959.1 - Maize associated partiti-like virus",
    "AOR51388.1 - Partitivirus-like 1",
    "QGZ98411.1 - Plasmopara viticola lesion associated Partiti-like 1"
  ),
  "Animal-derived" = c(
    "QIJ70087.1 - Varnsen partiti-like virus",
    "QIJ70088.1 - Pennypacker partiti-like virus",
    "QIJ70089.1 - Pennypacker partiti-like virus",
    "QIJ70076.1 - Costanza partiti-like virus",
    "QIJ70086.1 - Vandelay partiti-like virus"
  )
)
cols_apparli <-c("Animal-derived" = "#0B775E", "Partiti (plant/fungus)" = "#7FC0C6","Apparli virus Belgium" = "#F2300F","RNA Virosphere"="#35274A")

apparlitree <- groupOTU(apparlitree, cls_apparli)

apparlidraw <- ggtree(apparlitree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(align = TRUE, aes(subset = isTip & label %in% c(
    'Apparli virus', 'AEA40853.1 - uncultured virus', "APG78322.1 - Hubei partiti-like virus 34",
    "ATY36109.1 - Vespa velutina partiti-like virus 1", "BBE15516.1 - Osugoroshi virus 1",
    "APG78222.1 - Hubei partiti-like virus 33", 
    "QHA33899.1 - Atrato Partiti-like virus 3",
    "QHA33901.1 - Atrato Partiti-like virus 3",
    "QHA33703.1 - Atrato Partiti-like virus 4",
    "APG78254.1 - Hubei odonate virus 13",
    "APG78345.1 - Wenling partiti-like virus 3",
    "QGL51732.1 - Vespa velutina associated partiti-like virus 3",
    "APG78294.1 - Sanxia partiti-like virus 1",
    "APG78218.1 - Hubei partiti-like virus 48"
  ))) +
  geom_cladelabel(node=103, label="RNA virosphere / Animal-derived partiti-viruses", align=T, geom='label', col='#35274A') +
  geom_cladelabel(node=134, label="RNA virosphere / Animal-derived partiti-viruses", align=T, geom='label', col='#35274A', offset=-0.6) +
  geom_cladelabel(node=145, label="RNA virosphere / Plant&Fungus partiti viruses", align=T, geom='label', col='#35274A') +
  scale_colour_manual(values = cols_apparli) + 
  geom_treescale(2,13.5) +
  xlim(0, 10) +
  theme(legend.position="left") +
  ggtitle("Apparli virus") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=22)) +
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 6)))

apparlidraw <- scaleClade(apparlidraw, node=103, scale=.1) 
apparlidraw <- scaleClade(apparlidraw, node=134, scale=.1) 
apparlidraw <- scaleClade(apparlidraw, node=145, scale=.1) 
apparlidraw

apparlidraw_full <- ggtree(apparlitree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(align=TRUE) + 
  geom_cladelabel(node=103, label="RNA virosphere / Animal-derived partiti-viruses", align=T, geom='label', col='#35274A', offset=4) +
  geom_cladelabel(node=134, label="RNA virosphere / Animal-derived partiti-viruses", align=T, geom='label', col='#35274A', offset=4.5) +
  geom_cladelabel(node=145, label="RNA virosphere / Plant&Fungus partiti viruses", align=T, geom='label', col='#35274A', offset=4.9) +
  scale_colour_manual(values = cols_apparli) + 
  geom_treescale(offset=-1.5) +
  xlim(0, 10.2) +
  theme(legend.position="left") +
  ggtitle("Apparli virus") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=22)) +
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 6)))
apparlidraw_full
ggsave('output/sup_fig5.png', apparlidraw_full, dpi=300, width=12, height=12)


# BMLV --------------------------------------------------------------------

bmlv.acc2name <- read.delim("data/trees/bmlv.acc2name", header=FALSE)
bmlv.acc2name$V3 <- paste(bmlv.acc2name$V1, bmlv.acc2name$V2, sep=' - ')
bmlv.a2n <- c(paste(bmlv.acc2name$V1, bmlv.acc2name$V2, sep=' - '))
names(bmlv.a2n) <- bmlv.acc2name$V1

#bmlv
bmlvtree <- read.tree("data/trees/bmlv.treefile")
bmlvtree$tip.label <- as.character(bmlv.a2n[match(names(bmlv.a2n), bmlvtree$tip.label)])
bmlvtree$tip.label[1] <- "Bee Macula-like virus Belgium"
bmlvtree <- phangorn::midpoint(bmlvtree)
bmlvtree <- reorder(bmlvtree)

cls_bmlv <- list(
  "NODE" = c("Bee Macula-like virus Belgium"),
  "Tymovirus" = c(
    "YP_002308439.1 - Scrophularia mottle virus",
    "YP_002308578.1 - Anagyris vein yellowing virus",
    "YP_002308445.1 - Plantago mottle virus",
    "NP_041257.1 - Ononis yellow mosaic virus",
    "ARK20051.1 - Diascia yellow mottle virus",
    "YP_002048673.1 - Diascia yellow mottle virus",
    "YP_002308442.1 - Nemesia ring necrosis virus",
    "CBH31042.1 - Chiltepin yellow mosaic virus",
    "YP_003620401.1 - Chiltepin yellow mosaic virus",
    "AYC35256.1 - Naranjilla mild mosaic virus",
    "AIE44664.1 - Eggplant mosaic virus",
    "NP_040968.1 - Eggplant mosaic virus",
    "YP_007517183.1 - Andean potato mild mosaic virus",
    "YP_008318042.1 - Tomato blistering mosaic virus",
    "AIG05827.1 - Tomato blistering mosaic virus",
    "AMC38501.1 - Tomato blistering mosaic virus",
    "YP_007517180.1 - Andean potato latent virus",
    "NP_619756.1 - Physalis mottle virus",
    "AVW89220.1 - Passion fruit yellow mosaic virus",
    "YP_004464924.1 - Asclepias asymptomatic virus",
    "YP_001285472.1 - Okra mosaic virus",
    "AMH40128.1 - Turnip yellow mosaic virus",
    "AMH40131.1 - Turnip yellow mosaic virus",
    "AMH40125.1 - Turnip yellow mosaic virus",
    "AMH40122.1 - Turnip yellow mosaic virus",
    "AMH40137.1 - Turnip yellow mosaic virus",
    "AMH40140.1 - Turnip yellow mosaic virus",
    "AMH40134.1 - Turnip yellow mosaic virus",
    "P28477.1 - Turnip yellow mosaic virus (isolate TYMC)",
    "P10358.2 - Turnip yellow mosaic virus",
    "NP_663297.1 - Turnip yellow mosaic virus",
    "ALA65431.1 - Turnip yellow mosaic virus",
    "ALA65434.1 - Turnip yellow mosaic virus",
    "AIE44667.1 - Turnip yellow mosaic virus",
    "AAB92649.1 - Turnip yellow mosaic virus",
    "P20128.1 - Turnip yellow mosaic virus (Australian isolate)",
    "AFC95826.1 - Watercress white vein virus",
    "NP_047920.1 - Erysimum latent virus"
  ),
  "Marafivirus" = c(
    "AOX49247.1 - Grapevine rupestris vein feathering virus",
    "YP_009351862.1 - Grapevine rupestris vein feathering virus",
    "ATB20101.1 - Grapevine rupestris vein feathering virus",
    "AQV11970.1 - Grapevine rupestris vein feathering virus",
    "ATB20100.1 - Grapevine rupestris vein feathering virus",
    "AAW33732.1 - Grapevine rupestris vein feathering virus",
    "YP_002756536.1 - Grapevine Syrah virus 1",
    "AFV34757.1 - Grapevine Syrah virus 1",
    "ACV83739.1 - Grapevine virus Q",
    "AKZ17743.1 - Grapevine Syrah virus 1",
    "AKQ08149.1 - Grapevine Syrah virus 1",
    "ANP94321.1 - Grapevine Syrah virus 1",
    "AKQ08148.1 - Grapevine Syrah virus 1",
    "AKQ08147.1 - Grapevine Syrah virus 1",
    "AKZ17760.1 - Grapevine Syrah virus 1",
    "YP_004464920.1 - Switchgrass mosaic virus",
    "AIY22517.1 - Maize rayado fino virus",
    "NP_115454.1 - Maize rayado fino virus",
    "AAK52838.2 - Maize rayado fino virus",
    "QID59002.1 - Camellia-associated marafivirus",
    "AXM42964.1 - Grapevine asteroid mosaic associated virus",
    "ARH56444.1 - Grapevine asteroid mosaic associated virus",
    "AXM42965.1 - Grapevine asteroid mosaic associated virus",
    "QBZ78634.1 - Grapevine asteroid mosaic associated virus",
    "YP_009315883.1 - Grapevine asteroid mosaic associated virus",
    "ABA54133.1 - Citrus sudden death-associated virus",
    "ARO38249.1 - Citrus sudden death-associated virus",
    "YP_224218.1 - Citrus sudden death-associated virus",
    "AYC35261.1 - Citrus sudden death-associated virus",
    "ARO38248.1 - Citrus sudden death-associated virus",
    "YP_009222597.1 - Nectarine marafivirus M",
    "ALX72769.1 - Nectarine marafivirus M",
    "ALX72770.1 - Nectarine marafivirus M",
    "YP_009505639.1 - Blackberry virus S",
    "QCC30252.1 - Peach marafivirus D",
    "QCC30253.1 - Peach marafivirus D",
    "YP_009345914.1 - Peach virus D",
    "ATJ00054.1 - Medicago sativa marafivirus 1",
    "YP_009551972.1 - Alfalfa virus F"
  ),
  "Unclassified" = c(
    "YP_009160324.1 - Bee Macula-like virus",
    "AKQ48574.1 - Bee Macula-like virus",
    "AYV61007.1 - Guarapuava tymovirus-like 3",
    "ARI47198.1 - Peach virus T",
    "ARI47200.1 - Peach virus T",
    "AUR53413.1 - Naranjilla chlorotic mosaic virus",
    "AZF99024.1 - Ullucus tymovirus 2",
    "AEP40395.1 - Tomato yellow blotch virus",
    "BAJ14669.1 - Poinsettia mosaic virus",
    "BAJ14667.1 - Poinsettia mosaic virus",
    "BAJ14673.1 - Poinsettia mosaic virus",
    "BAJ14665.1 - Poinsettia mosaic virus",
    "NP_037647.1 - Poinsettia mosaic virus",
    "CAL80776.1 - Poinsettia mosaic virus",
    "QBP37021.1 - Lampyris noctiluca tymovirus-like virus 1",
    "YP_009505642.1 - Bombyx mori latent virus",
    "YP_009159826.1 - Varroa Tymo-like virus"
  ),
  "Maculavirus" = c(
    "ANV22069.1 - Grapevine Red Globe virus",
    "YP_009268923.1 - Grapevine Red Globe virus",
    "ANV22095.1 - Grapevine Red Globe virus",
    "AOX49245.1 - Grapevine Red Globe virus",
    "YP_004464930.1 - Bombyx mori Macula-like virus",
    "AHX22590.1 - Bombyx mori Macula-like virus"
  )
)

cols_bmlv <-c("Tymovirus" = "#0B775E","Unclassified"= "#35274A", "Marafivirus" = "#EABE94","NODE" = "#F2300F","Maculavirus"="#02401B")

bmlvtree <- groupOTU(bmlvtree, cls_bmlv)

bmlvdraw <- ggtree(bmlvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(aes(subset = isTip & label %in% 
                    c(
                      "Bee Macula-like virus Belgium",
                      "YP_009160324.1 - Bee Macula-like virus",
                      "AKQ48574.1 - Bee Macula-like virus",
                      "AYV61007.1 - Guarapuava tymovirus-like 3",
                      "BAJ14669.1 - Poinsettia mosaic virus",
                      "BAJ14667.1 - Poinsettia mosaic virus",
                      "BAJ14673.1 - Poinsettia mosaic virus",
                      "BAJ14665.1 - Poinsettia mosaic virus",
                      "NP_037647.1 - Poinsettia mosaic virus",
                      "CAL80776.1 - Poinsettia mosaic virus",
                      "QBP37021.1 - Lampyris noctiluca tymovirus-like virus 1",
                      "YP_009159826.1 - Varroa Tymo-like virus"
                    )
                    )) +
  #geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_bmlv) + 
  geom_cladelabel(node=149, label="Marafivirus", align=T, geom='label', col='#EABE94', offset = -0.7) +
  geom_cladelabel(node=154, label="Tymovirus", align=T, geom='label', col='#0B775E') +
  geom_cladelabel(node=110, label="Marafivirus", align=T, geom='label', col='#EABE94', offset = -0.5) +
  geom_cladelabel(node=150, label="Maculavirus", align=T, geom='label', col='#02401B', offset=-0.6) +
  geom_cladelabel(node=199, label="Maculavirus", align=T, geom='label', col='#02401B', offset=-0.65) +
  geom_treescale(2,14.5) +
  ggtitle("Bee-macula like virus") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=22)) +
  xlim(c(0,8)) +
  theme(legend.position = 'none')
bmlvdraw <- scaleClade(bmlvdraw, node=154, scale=.1) 
bmlvdraw <- scaleClade(bmlvdraw, node=110, scale=.1)
bmlvdraw <- scaleClade(bmlvdraw, node=149, scale=.1)
bmlvdraw

bmlvdraw_full <- ggtree(bmlvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_bmlv) + 
  geom_cladelabel(node=149, label="Marafivirus", align=T, geom='label', col='#EABE94', offset = 1.2) +
  geom_cladelabel(node=154, label="Tymovirus", align=T, geom='label', col='#0B775E', offset = 1.2) +
  geom_cladelabel(node=110, label="Marafivirus", align=T, geom='label', col='#EABE94', offset = 1.2) +
  geom_cladelabel(node=150, label="Maculavirus", align=T, geom='label', col='#02401B', offset=1.2) +
  geom_cladelabel(node=199, label="Maculavirus", align=T, geom='label', col='#02401B', offset=1.2) +
  geom_treescale(offset=-1.5) +
  ggtitle("Bee-macula like virus") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=22)) +
  xlim(c(0,3.5)) +
  theme(legend.position = 'none')
bmlvdraw_full
ggsave('output/sup_fig4.png', bmlvdraw_full, dpi=300, width=12, height=12)


# AMFV --------------------------------------------------------------------

amfv.acc2name <- read.delim("data/trees/amfv.acc2name", header=FALSE)
amfv.acc2name$V3 <- paste(amfv.acc2name$V1, amfv.acc2name$V2, sep=' - ')
amfv.a2n <- c(paste(amfv.acc2name$V1, amfv.acc2name$V2, sep=' - '))
names(amfv.a2n) <- amfv.acc2name$V1

#amfv
amfvtree <- read.tree("data/trees/amfv.treefile")
amfvtree$tip.label <- as.character(amfv.a2n[match(names(amfv.a2n), amfvtree$tip.label)])
amfvtree$tip.label[1] <- "Apis mellifera filamentous virus Belgium"
amfvtree <- phangorn::midpoint(amfvtree)
amfvtree <- reorder(amfvtree)

cls_amfv <- list(
  "NODE" = c("Apis mellifera filamentous virus Belgium"),
  "Poxviridae" = c(
    "QDJ94996.1 - Hypsugopox virus",
    "YP_009407970.1 - Eptesipox virus",
    "YP_004821373.1 - Yokapox virus",
    "YP_009408420.1 - NY_014 poxvirus",
    "YP_009408218.1 - Murmansk poxvirus",
    "YP_009143350.1 - Raccoonpox virus",
    "YP_009282734.1 - Skunkpox virus",
    "YP_009281788.1 - Volepox virus",
    "NP_619839.1 - Cowpox virus",
    "ARR30425.1 - Cowpox virus",
    "ARR30628.1 - Cowpox virus",
    "ARR31017.1 - Cowpox virus",
    "ADZ30447.1 - Cowpox virus",
    "ADZ30659.1 - Cowpox virus",
    "ARR30008.1 - Cowpox virus",
    "SNB53815.1 - Cowpox virus",
    "ARR29609.2 - Cowpox virus",
    "ADZ29378.1 - Cowpox virus",
    "AGY97303.1 - Cowpox virus",
    "AAA60776.1 - Variola major virus",
    "P0DSS8.1 - Variola virus",
    "AAA69333.1 - Variola virus",
    "AAA69374.1 - Variola virus",
    "YP_717351.1 - Taterapox virus",
    "NP_570429.1 - Camelpox virus",
    "SNB49884.1 - Cowpox virus",
    "ADZ29593.1 - Cowpox virus",
    "QED21168.1 - Alaskapox virus",
    "AGY97736.1 - Cowpox virus",
    "O57175.1 - Vaccinia virus Ankara",
    "ABD97395.1 - Cowpox virus",
    "AAF33894.1 - Vaccinia virus Tian Tan",
    "ADZ29806.1 - Cowpox virus",
    "YP_232925.1 - Vaccinia virus",
    "AGJ91195.1 - Vaccinia virus",
    "AEY72851.1 - Vaccinia virus",
    "AAS49745.1 - Rabbitpox virus",
    "AAW23440.1 - Vaccinia virus",
    "AGJ91739.1 - Vaccinia virus",
    "ARB50277.1 - Cowpox virus",
    "ABH08147.1 - Horsepox virus",
    "P20493.1 - Vaccinia virus Copenhagen",
    "AXN75271.1 - Akhmeta virus",
    "AXN74832.1 - Akhmeta virus",
    "QEJ79168.1 - Goatpox virus",
    "QEJ79318.1 - Goatpox virus",
    "QEJ79019.1 - Goatpox virus",
    "YP_001293211.1 - Goatpox virus Pellor",
    "AAN02745.1 - Lumpy skin disease virus",
    "AAK43560.1 - Lumpy skin disease virus",
    "AGZ95335.1 - Goatpox virus FZ",
    "NP_659593.1 - Sheeppox virus",
    "NP_150454.1 - Lumpy skin disease virus NI-2490",
    "NP_073405.1 - Yaba-like disease virus",
    "ABQ43492.1 - Tanapox virus",
    "NP_938276.1 - Yaba monkey tumor virus",
    "YP_005296225.1 - Cotia virus SPAn232",
    "YP_009329641.1 - BeAn 58058 virus",
    "AUI80588.1 - White-tailed deer poxvirus",
    "ABI99013.1 - Deerpox virus W-1170-84",
    "YP_227403.1 - Deerpox virus W-848-83",
    "ACB28637.1 - Myxoma virus",
    "NP_051729.1 - Myxoma virus",
    "QAV34328.1 - Myxoma virus",
    "QAV42102.1 - Myxoma virus",
    "AFU76948.1 - Myxoma virus",
    "AMB18349.1 - Myxoma virus",
    "AGU99698.1 - Myxoma virus",
    "NP_051904.1 - Rabbit fibroma virus",
    "AYP74130.1 - Fowlpox virus",
    "YP_009268721.1 - Pteropox virus"
  ),
  "Unclassified" = c("YP_009165967.1 - Apis mellifera filamentous virus",
                     "ADD74390.1 - uncultured virus"),
  "Phycodnaviridae" = c(
    "YP_293780.1 - Emiliania huxleyi virus 86",
    "AEO97836.1 - Emiliania huxleyi virus 84",
    "CAZ69365.1 - Emiliania huxleyi virus 99B1",
    "AEO98243.1 - Emiliania huxleyi virus 203"
  ),
  "Baculoviridae" = c(
    "NP_258331.1 - Spodoptera litura nucleopolyhedrovirus",
    "YP_009505866.1 - Spodoptera littoralis nucleopolyhedrovirus",
    "AHC69623.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AJR20393.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "QCQ67370.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AIX47958.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "QDH05963.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AQQ80135.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "QCQ67530.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AMO27788.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AWJ76733.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AMO27616.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "NP_047757.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AMO27974.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "AOW42794.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "QCQ67689.1 - Lymantria dispar multiple nucleopolyhedrovirus",
    "YP_003517877.1 - Lymantria xylina nucleopolyhedrovirus",
    "ANS71004.1 - Lymantria dispar multiple nucleopolyhedrovirus"
  ),
  "Other_families" = c(
    "AXN91020.1 - Namao virus",
    "QCQ67809.1 - European chub iridovirus",
    "AUS94194.1 - Trichoplusia ni ascovirus 6b",
    "YP_001883390.1 - Musca domestica salivary gland hypertrophy virus",
    "BAV31392.1 - Tenacibaculum phage pT24"
  )
)

cols_amfv <-c("Poxviridae" = "#0B775E","Unclassified"= "#35274A", "Phycodnaviridae" = "#EABE94","NODE" = "#F2300F","Baculoviridae"="#02401B","Other_families"="#7f7f7f")
amfvtree <- groupOTU(amfvtree, cls_amfv)

amfvdraw <- ggtree(amfvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(aes(subset = isTip & label %in% 
                    c(
                      "Apis mellifera filamentous virus Belgium",
                      "YP_001883390.1 - Musca domestica salivary gland hypertrophy virus",
                      "BAV31392.1 - Tenacibaculum phage pT24",
                      "YP_009165967.1 - Apis mellifera filamentous virus",
                      "ADD74390.1 - uncultured virus"
                    )
  )) +
  #geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  #geom_text(aes(label=node)) +
  scale_colour_manual(values = cols_amfv) + 
  geom_treescale(2,9.5) +
  geom_cladelabel(node=108, label="Poxviridae", align=T, geom='label', col='#0B775E') +
  geom_cladelabel(node=180, label="Baculoviridae", align=T, geom='label', col='#02401B', offset = -0.2) +
  geom_cladelabel(node=198, label="Phycodnaviridae", align=T, geom='label', col='#EABE94', offset=-0.35) +
  ggtitle("Apis Melifera filamentous virus") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=22)) +
  xlim(0, 6 ) +
  theme(legend.position = 'none')
amfvdraw <- scaleClade(amfvdraw, node=108, scale=.1) 
amfvdraw <- scaleClade(amfvdraw, node=180, scale=.1)
amfvdraw <- scaleClade(amfvdraw, node=198, scale=.1)
amfvdraw

amfvdraw_full <- ggtree(amfvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  #geom_text(aes(label=node)) +
  scale_colour_manual(values = cols_amfv) + 
  geom_treescale(offset=-1.5) +
  geom_cladelabel(node=108, label="Poxviridae", align=T, geom='label', col='#0B775E', offset=0.5) +
  geom_cladelabel(node=180, label="Baculoviridae", align=T, geom='label', col='#02401B', offset = 0.65) +
  geom_cladelabel(node=198, label="Phycodnaviridae", align=T, geom='label', col='#EABE94', offset=0.45) +
  ggtitle("Apis Melifera filamentous virus") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=22)) +
  xlim(0, 2 ) +
  theme(legend.position = 'none')
amfvdraw_full
ggsave('output/sup_fig3.png', amfvdraw_full, dpi=300, width=12, height=12)

# dwv - vdv ---------------------------------------------------------------

vdvdwv.acc2name <- read.delim("data/trees/vdv_dwv.acc2name", header=FALSE)
vdvdwv.acc2name$V3 <- paste(vdvdwv.acc2name$V1, vdvdwv.acc2name$V2, sep=' - ')
vdvdwv.a2n <- c(paste(vdvdwv.acc2name$V1, vdvdwv.acc2name$V2, sep=' - '))
names(vdvdwv.a2n) <- vdvdwv.acc2name$V1

#vdv
vdvtree <- read.tree("data/trees/vdv_dwv.treefile")
vdvtree$tip.label <- as.character(vdvdwv.a2n[match(names(vdvdwv.a2n), vdvtree$tip.label)])
vdvtree$tip.label[1] <- "Deformed wing virus Belgium"
vdvtree$tip.label[2] <- "Varroa destructor virus 1 Belgium"
vdvtree <- phangorn::midpoint(vdvtree)
vdvtree <- reorder(vdvtree)

cls_vdv <- list('NODE'=c("Deformed wing virus Belgium", "Varroa destructor virus 1 Belgium" ),
                'Iflaviridae' = c("YP_145791.1 - Varroa destructor virus 1",
                                  "YP_009162630.1 - Bombyx mori iflavirus",
                                  "YP_009002581.1 - Antheraea pernyi iflavirus",
                                  "YP_009116875.1 - Thaumetopoea pityocampa iflavirus 1",
                                  "YP_009047245.1 - Lymantria dispar iflavirus 1",
                                  "YP_009026409.1 - Heliconius erato iflavirus",
                                  "YP_009344960.1 - Helicoverpa armigera iflavirus",
                                  "YP_009315906.1 - King virus",
                                  "YP_009345906.1 - Bat iflavirus",
                                  "YP_009305421.1 - Moku virus",
                                  "YP_003622540.1 - Slow bee paralysis virus",
                                  "YP_009328891.1 - Euscelidius variegatus virus 1",
                                  "YP_009129265.1 - Graminella nigrifrons virus 1",
                                  "YP_009553259.1 - Psammotettix alienus iflavirus 1",
                                  "YP_009351892.1 - Pityohyphantes rubrofasciatus iflavirus",
                                  "YP_009552080.1 - Yongsan iflavirus 1",
                                  "YP_009552119.1 - Varroa destructor virus 2",
                                  "YP_009444707.1 - Chequa iflavirus",
                                  "YP_009010941.1 - Laodelphax striatella honeydew virus 1",
                                  "YP_009505599.1 - Nilaparvata lugens honeydew virus 1",
                                  "YP_008130309.1 - Nilaparvata lugens honeydew virus-2",
                                  "YP_008130310.1 - Nilaparvata lugens honeydew virus-3",
                                  "YP_008888537.1 - Formica exsecta virus 2",
                                  "NP_853560.2 - Deformed wing virus"),
                'Unclassified' = c("YP_009337003.1 - Hubei picorna-like virus 26",
                                   "YP_009337760.1 - Hubei odonate virus 4",
                                   "YP_009337271.1 - Hubei tetragnatha maxillosa virus 2",
                                   "YP_009337284.1 - Hubei picorna-like virus 28",
                                   "YP_009329861.1 - Wuhan spider virus 2",
                                   "YP_009330055.1 - Hubei picorna-like virus 31",
                                   "YP_009337046.1 - Shuangao insect virus 12",
                                   "YP_009336533.1 - Hubei tick virus 3",
                                   "YP_009336552.1 - Hubei tick virus 1",
                                   "YP_009330050.1 - Hubei myriapoda virus 1",
                                   "YP_009342053.1 - Wuhan coneheads virus 1",
                                   "YP_009342321.1 - Wuhan insect virus 13",
                                   "YP_009333199.1 - Shahe heteroptera virus 2",
                                   "YP_009336939.1 - Shahe heteroptera virus 1",
                                   "YP_009336575.1 - Hubei picorna-like virus 29",
                                   "YP_009337127.1 - Hubei coleoptera virus 1",
                                   "YP_009337161.1 - Hubei picorna-like virus 27"))

cols_vdv <-c("Iflaviridae" = "#EABE94","Unclassified"= "#35274A", "NODE" = "#F2300F")
vdvtree <- groupOTU(vdvtree, cls_vdv)

vdvdraw <- ggtree(vdvtree, aes(color=group)) + 
  geom_tiplab(aes(subset = isTip & label %in% c(
    "YP_145791.1 - Varroa destructor virus 1",
    "Deformed wing virus Belgium", "Varroa destructor virus 1 Belgium",
    "YP_008130310.1 - Nilaparvata lugens honeydew virus-3",
    "YP_008888537.1 - Formica exsecta virus 2",
    "YP_009328891.1 - Euscelidius variegatus virus 1",
    "YP_009129265.1 - Graminella nigrifrons virus 1",
    "YP_009553259.1 - Psammotettix alienus iflavirus 1",
    "YP_009337271.1 - Hubei tetragnatha maxillosa virus 2",
    "YP_009337284.1 - Hubei picorna-like virus 28",
    "YP_009345906.1 - Bat iflavirus",
    "YP_009305421.1 - Moku virus",
    "YP_003622540.1 - Slow bee paralysis virus",
    "YP_009351892.1 - Pityohyphantes rubrofasciatus iflavirus",
    "YP_009552080.1 - Yongsan iflavirus 1",
    "YP_009337161.1 - Hubei picorna-like virus 27",
    "YP_009337760.1 - Hubei odonate virus 4",
    "YP_009337127.1 - Hubei coleoptera virus 1",
    "NP_853560.2 - Deformed wing virus"))) +
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  #geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_vdv)+
  geom_treescale(2,2.5) +
  ggtitle("DWV/VDV") + 
  xlim(0, 6)+
  theme(legend.position=c(0.8,0.5)) +
  theme(plot.title = element_text(size=22)) +
  theme(legend.position="left") +
  ggtitle("VDV/DWV") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Viral group",values=c("#EABE94","#F2300F","#35274A"), labels=c("Iflaviridae","Iflavirus Belgium","Unclassified")) +
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 6)))
vdvdraw <- scaleClade(vdvdraw, node=69, scale=.1) 
vdvdraw <- scaleClade(vdvdraw, node=51, scale=.1) 
vdvdraw

vdvdraw_full <- ggtree(vdvtree, aes(color=group)) + 
  geom_nodepoint(aes(subset = as.numeric(label) > 70), color='black',size=2) +
  geom_tiplab(align=TRUE, show.legend=FALSE, size=3)+
  scale_colour_manual(values = cols_vdv)+
  geom_treescale(offset=-1.5) +
  ggtitle("DWV/VDV") + 
  xlim(0, 4)+
  theme(legend.position=c(0.8,0.5)) +
  theme(plot.title = element_text(size=22)) +
  theme(legend.position="left") +
  ggtitle("VDV/DWV") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="Viral group",values=c("#EABE94","#F2300F","#35274A"), labels=c("Iflaviridae","Iflavirus Belgium","Unclassified")) +
  guides(color = guide_legend(override.aes = list(label = "\u25A0", size = 6)))
vdvdraw_full
ggsave('output/sup_fig2.png', vdvdraw_full, dpi=300, width=12, height=12)



# Combine -----------------------------------------------------------------

testgrid <- grid.arrange(amfvdraw,bmlvdraw,apparlidraw,vdvdraw, apthilidraw ,orthodraw,nrow=3, ncol=2)
ggsave("output/fig1.png",testgrid, device = "png", height=15,width=20, limitsize=FALSE)
