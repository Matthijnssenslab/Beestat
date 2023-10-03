suppressMessages(library('caret'))
suppressMessages(library('ROCR'))
suppressMessages(library('tidyr'))
# Prep data
qpcr <- read.delim("data/qpcr_abs.tsv", row.names=1)
metadf <- read.delim("data/metadata.tsv", row.names=1)
metadf$colony_id <- as.factor(metadf$colony_id)
metadf$status <- as.factor(metadf$status)
metadf$status_year <- as.factor(metadf$status_year)
metadf <- metadf %>% separate(status_year, '_', into=c('status', 'year'))
metadf$year <- as.factor(as.character(metadf$year))
metadf$sample <- rownames(metadf)
qpcr$sample <- rownames(qpcr)
df <- merge(qpcr, metadf, by='sample', all.x = TRUE)

colnames(df) <- c(
  'sample',
  'dwv','aov', 'apthili', 'bmlv', 'aparli', 'amfv', 'vdv',
  'colony_id', 'lat', 'lon', 'commune', 'status', 'year'
)

rownames(df) <- df$sample
df$sample <- NULL

df$dwv <- log10(df$dwv + 1)
df$aov <- log10(df$aov + 1)
df$apthili <- log10(df$apthili + 1)
df$bmlv <- log10(df$bmlv + 1)
df$aparli <- log10(df$aparli + 1)
df$amfv <- log10(df$amfv + 1)
df$vdv <- log10(df$vdv + 1)
df$statusbool <- factor(ifelse(df$status == 'wl', 1, 0))
df$commune <- factor(df$commune)
df$year <- factor(df$year)

logdf <- df[, c(
  'statusbool', 'year', 'dwv', 'aov', 'apthili', 'bmlv', 'aparli', 'amfv' ,'vdv', 'commune'
)]

# caret
cv_full <- train(
  statusbool ~ year + (dwv + aov + apthili + bmlv + aparli + amfv + vdv)^2, 
  data = logdf, 
  method = "glm",
  family = "binomial",
  trControl = trainControl(method = "cv", number = 10)
)

cv_full_sum <- data.frame(summary(cv_full)$coefficients)
cv_full_or <- data.frame(cbind(exp(coef(cv_full$finalModel)), exp(confint(cv_full$finalModel))))
cv_full_or_coef <- merge(cv_full_sum, cv_full_or, by.x = 0, by.y = 0)
rownames(cv_full_or_coef) <- cv_full_or_coef$Row.names
cv_full_or_coef$Row.names <- NULL
cv_full_or_coef <- cv_full_or_coef[-1,]
colnames(cv_full_or_coef) <- c('Estimate', 'se', 'z.value', 'pval', 'OR', 'OR_lower_95', 'OR_upper_95')
write.table(cv_full_or_coef, 'data/out_cvfull_or.tsv', sep='\t')

cv_yearassoc <- train(
  statusbool ~  year+ (dwv + apthili + vdv)^2, 
  data = logdf, 
  method = "glm",
  family = "binomial",
  trControl = trainControl(method = "cv", number = 10)
)

cv_yearassoc_sum <- data.frame(summary(cv_yearassoc)$coefficients)
cv_yearassoc_or <- data.frame(cbind(exp(coef(cv_yearassoc$finalModel)), exp(confint(cv_yearassoc$finalModel))))
cv_yearassoc_or_coef <- merge(cv_yearassoc_sum, cv_yearassoc_or, by.x = 0, by.y = 0)
rownames(cv_yearassoc_or_coef) <- cv_yearassoc_or_coef$Row.names
cv_yearassoc_or_coef$Row.names <- NULL
cv_yearassoc_or_coef <- cv_yearassoc_or_coef[-1,]
colnames(cv_yearassoc_or_coef) <- c('Estimate', 'se', 'z.value', 'pval', 'OR', 'OR_lower_95', 'OR_upper_95')
write.table(cv_yearassoc_or_coef, 'data/out_cvyearassoc_or.tsv', sep='\t')

cv_simple <- train(
  statusbool ~ year + dwv + aov + apthili + bmlv + aparli + amfv + vdv,
  data=logdf,
  method='glm',
  family='binomial',
  trControl = trainControl(method = 'cv', number=10)
)

cv_simple_sum <- data.frame(summary(cv_simple)$coefficients)
cv_simple_or <- data.frame(cbind(exp(coef(cv_simple$finalModel)), exp(confint(cv_simple$finalModel))))
cv_simple_or_coef <- merge(cv_simple_sum, cv_simple_or, by.x = 0, by.y = 0)
rownames(cv_simple_or_coef) <- cv_simple_or_coef$Row.names
cv_simple_or_coef$Row.names <- NULL
cv_simple_or_coef <- cv_simple_or_coef[-1,]
colnames(cv_simple_or_coef) <- c('Estimate', 'se', 'z.value', 'pval', 'OR', 'OR_lower_95', 'OR_upper_95')
write.table(cv_simple_or_coef, 'data/out_cvsimple_or.tsv', sep='\t')


pred_full <- predict(cv_full, logdf)
pred_year_assoc <- predict(cv_yearassoc, logdf)
pred_simple <- predict(cv_simple, logdf)

full_mat <- confusionMatrix(
  data = relevel(pred_full, ref = '0'), 
  reference = relevel(logdf$statusbool, ref = '0')
)
yearassoc_mat <- confusionMatrix(
  data = relevel(pred_year_assoc, ref = '0'), 
  reference = relevel(logdf$statusbool, ref = '0')
)
simple_mat <- confusionMatrix(
  data = relevel(pred_simple, ref = '0'), 
  reference = relevel(logdf$statusbool, ref = '0')
)

# Write results.
confdf <- cbind(
  rbind(
    as.matrix(full_mat, what='classes'),
    as.matrix(full_mat, what='overall')
  ),
  rbind(
    as.matrix(yearassoc_mat, what='classes'),
    as.matrix(yearassoc_mat, what='overall')
  ),
  rbind(
    as.matrix(simple_mat, what='classes'),
    as.matrix(simple_mat, what='overall')
  )
)
colnames(confdf) <- c('full', 'assoc', 'simple')
write.table(confdf, 'data/out_conf.tsv', sep='\t')


prob_full <- predict(cv_full, logdf, type='prob')$'1'
prob_year_assoc <- predict(cv_yearassoc, logdf, type='prob')$'1'
prob_simple <- predict(cv_simple, logdf, type='prob')$'1'


perf_full <- prediction(prob_full, logdf$statusbool) %>%
  performance(measure = "tpr", x.measure = "fpr")
perf_yearassoc <- prediction(prob_year_assoc, logdf$statusbool) %>%
  performance(measure = "tpr", x.measure = "fpr")
perf_simple <- prediction(prob_simple, logdf$statusbool) %>%
  performance(measure = "tpr", x.measure = "fpr")

rocdf <- cbind(
  perf_full@x.values[[1]],
  perf_full@y.values[[1]],
  perf_yearassoc@x.values[[1]],
  perf_yearassoc@y.values[[1]],
  perf_simple@x.values[[1]],
  perf_simple@y.values[[1]]
)
colnames(rocdf) <- c('fullfpr', 'fulltpr', 'assocfpr', 'assoctpr', 'simplefpr', 'simpletpr')
write.table(rocdf, 'data/out_roc.tsv', sep='\t')

