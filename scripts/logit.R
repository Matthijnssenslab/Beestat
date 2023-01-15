suppressMessages(library('glmnet'))

qpcr <- read.delim("data/qpcr_abs.tsv", row.names=1)
metadf <- read.delim("data/metadata.tsv", row.names=1)
metadf$colony_id <- as.factor(metadf$colony_id)
metadf$status <- as.factor(metadf$status)
metadf$status_year <- as.factor(metadf$status_year)
metadf$year <- as.factor(as.character(metadf$year))
metadf$sample <- rownames(metadf)

qpcr$sample <- rownames(qpcr)


df <- merge(qpcr, metadf, by='sample', all.x = TRUE)

colnames(df) <- c(
  'sample',
  'dwv','thog', 'thika', 'bmlv', 'unc', 'amfv', 'vdv',
  'colony_id', 'status', 'lat', 'lon', 'year', 'status_year'
)
rownames(df) <- df$sample
df$sample <- NULL

df$dwv <- log10(df$dwv + 1)
df$thog <- log10(df$thog + 1)
df$thika <- log10(df$thika + 1)
df$bmlv <- log10(df$bmlv + 1)
df$unc <- log10(df$unc + 1)
df$amfv <- log10(df$amfv + 1)
df$vdv <- log10(df$vdv + 1)


# 2012
set.seed(123)
x <- model.matrix(
  ~ year + vdv + dwv + thog + thika + bmlv + unc + amfv,
  df
)
xint <- model.matrix(
  ~ (year + vdv + dwv + thog + thika + bmlv + unc + amfv)^2,
  df
)

y <- ifelse(df$status == 'wl', 1, 0)

# 1 = dead, 0 = healthy
set.seed(123)
cv.mod <- cv.glmnet(x, y,family='binomial', type.measure='auc', keep=TRUE, alpha=00, lambda=NULL)
cv.mod2 <- cv.glmnet(xint, y,family='binomial', type.measure='auc', keep=TRUE, alpha=00, lambda=NULL)

coef(cv.mod)
plot(cv.mod)

rocs <- roc.glmnet(cv.mod2$fit.preval, newy = y)
best <- cv.mod$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "blue")

# 2013
set.seed(123)
x <- model.matrix(
  ~ (dwv + thog + thika + bmlv + unc + amfv)^2,
  df[df$year == '13',]
)
y <- ifelse(df[df$year == '13',]$status == 'wl', 1, 0)

# 1 = dead, 0 = healthy
set.seed(123)
cv.mod <- cv.glmnet(x, y,family='binomial', type.measure='auc', keep=TRUE, alpha=0.8, lambda=NULL)
plot(cv.mod)

rocs <- roc.glmnet(cv.mod$fit.preval, newy = y)
best <- cv.mod$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")
