library("Boruta")
library("caret")
library("rattle")
library("rpart")
library("factoextra")
# Reading_data ------------------------------------------------------------

#Stats:
vircast <- read.csv("../Data/Virdf_metadata_cast.csv", row.names=1)
#Weak = 1, Healthy= 0
vircast$status <- as.factor(vircast$status)
vircast$Numvir <- as.numeric(vircast$Numvir)
# Quick PCA ---------------------------------------------------------------
Status <- vircast_important$status
princomp <- prcomp(as.matrix(vircast[,-1]), scale = TRUE, center = TRUE)
get_eigenvalue(princomp)
PCAplot <- ggplot() + 
  geom_point(aes(x=princomp$x[,1], y=princomp$x[,2], col=Status)) +
  theme_minimal() +
  labs(x="PC1 (10.0%)", y="PC2 (6.72%)", fill="Status") +
  scale_color_manual(values=c("#2ca02c", "#d62728"))
PCAplot
ggsave("../Plots/PCAplot.pdf",PCAplot,dpi=300)
# Boruta ------------------------------------------------------------------

#Feature selection:
Borutaresults <- Boruta(status~., data=vircast)

#Confirmed:
Borutaresults$finalDecision[Borutaresults$finalDecision == 'Confirmed']
#Tentative:
Borutaresults$finalDecision[Borutaresults$finalDecision == 'Tentative']

keepcols <- c("status","mites","WIN_None","WIN_OA","AMFVload","thikaload","vdvload","dwvload","PRE_None","WIN_Chlorophenyl","Numvir")
vircast_important <- vircast[,colnames(vircast) %in% keepcols]


# Running all predictors --------------------------------------------------


#logit:
set.seed(200)
Train <- createDataPartition(vircast$status, p=0.8, list=FALSE)
training <- vircast[ Train, ]
testing <- vircast[ -Train, ]
#Exclude variables with 0 variance
training <- training[,-nearZeroVar(training)]
testing <- testing[,-nearZeroVar(training)]

#default logit
set.seed(200)
default_glm = train(
  form = status ~ .,
  data = training,
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  method = "glm",
  family = "binomial",
  tuneLength=10
)
#randomforest
set.seed(200)
default_rf = train(
  form = status ~ .,
  data = training,
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  method = "rf",
  family = "binomial",
  tuneLength=10
)
#knn
set.seed(200)
default_knn = train(
  status ~ .,
  data = training,
  method = "knn",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#rpart
set.seed(200)
default_rpart <- train(
  status ~.,
  data = training, 
  method = "rpart",
  parms = list(split = "information"),
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 10
)
#lvq
set.seed(200)
default_lvq <- train(
  status ~.,
  data = training, 
  method = "lvq",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 10
)
#gbm
set.seed(200)
default_gbm<- train(
  status ~.,
  data = training, 
  method = "gbm",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 10
)

#svm
set.seed(200)
default_svm <- train(
  status ~.,
  data = training,
  method = "svmRadial",
  preProc=c("center","scale"),
  trControl=trainControl(method = "repeatedcv",number = 10, repeats=5),
  tuneLength = 10
)
#knn
set.seed(200)
default_kknn = train(
  status ~ .,
  data = training,
  method = "kknn",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#pda
set.seed(200)
default_pda = train(
  status ~ .,
  data = training,
  method = "pda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#sda
set.seed(200)
default_sda = train(
  status ~ .,
  data = training,
  method = "sda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#slda
set.seed(200)
default_slda = train(
  status ~ .,
  data = training,
  method = "slda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#hdda
set.seed(200)
default_hdda = train(
  status ~ .,
  data = training,
  method = "hdda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#pam
set.seed(200)
default_pam = train(
  status ~ .,
  data = training,
  method = "pam",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#SIMCA
set.seed(200)
default_simca = train(
  status ~ .,
  data = training,
  method = "CSimca",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#pls
set.seed(200)
default_pls = train(
  status ~ .,
  data = training,
  method = "pls",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#glmnet
set.seed(200)
default_glmnet = train(
  status ~ .,
  data = training,
  method = "glmnet",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#nnet
set.seed(200)
default_nnet = train(
  form = status ~ .,
  data = training,
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  method = "nnet",
  tuneLength=10
)
#gbm
set.seed(200)
default_gbm<- train(
  status ~.,
  data = training, 
  method = "gbm",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 100
)


# Grid optimize GBM_ALL ---------------------------------------------------


#Since GBM performs best, try to optimize:
gbmoptgrid <-  expand.grid(interaction.depth = c(1, 5, 9, 13, 17, 21, 25, 29), 
                           n.trees = (1:100)*50, 
                           shrinkage = c(0.0001,0.001,0.01,0.1),
                           n.minobsinnode = c(1,5,10,20,50,100))
set.seed(200)
default_gbm_tuned<- train(
  status ~.,
  data = training, 
  method = "gbm",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneGrid = gbmoptgrid
)


# Saving dotplot ----------------------------------------------------------


results <- resamples(list(GBM_tuned=default_gbm_tuned,NNET=default_nnet,GLMNET=default_glmnet,PLS=default_pls,SIMCA=default_simca,PAM=default_pam,HDDA=default_hdda,SLDA=default_slda,SDA=default_sda,PDA=default_pda,KKNN=default_kknn,SVM=default_svm,GBM=default_gbm,LVQ=default_lvq,RPART=default_rpart,KNN=default_knn,RF=default_rf,GLM=default_glm))
dotplot_allpred <- dotplot(results, main="All predictors")
trellis.device(device="png", filename="../Plots/Dotplot_allpred.png")
print(dotplot_allpred)
dev.off()


# Running featsel ---------------------------------------------------------


#1 more on the selected features:
#logit:
set.seed(200)
Train <- createDataPartition(vircast_important$status, p=0.8, list=FALSE)
training <- vircast_important[ Train, ]
testing <- vircast_important[ -Train, ]
#Exclude variables with 0 variance
training <- training[,-nearZeroVar(training)]
testing <- testing[,-nearZeroVar(training)]
#default logit
set.seed(200)
default_glm_imp = train(
  form = status ~ .,
  data = training,
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  method = "glm",
  family = "binomial",
  tuneLength=10
)
#randomforest
set.seed(200)
default_rf_imp = train(
  form = status ~ .,
  data = training,
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  method = "rf",
  family = "binomial",
  tuneLength=10
)
#knn
set.seed(200)
default_knn_imp = train(
  status ~ .,
  data = training,
  method = "knn",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#rpart
set.seed(200)
default_rpart_imp <- train(
  status ~.,
  data = training, 
  method = "rpart",
  parms = list(split = "information"),
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 10
)
#lvq
set.seed(200)
default_lvq_imp <- train(
  status ~.,
  data = training, 
  method = "lvq",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 10
)
#gbm
set.seed(200)
default_gbm_imp<- train(
  status ~.,
  data = training, 
  method = "gbm",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 10
)

#svm
set.seed(200)
default_svm_imp <- train(
  status ~.,
  data = training,
  method = "svmRadial",
  preProc=c("center","scale"),
  trControl=trainControl(method = "repeatedcv",number = 10, repeats=5),
  tuneLength = 10
)
#knn
set.seed(200)
default_kknn_imp = train(
  status ~ .,
  data = training,
  method = "kknn",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#pda
set.seed(200)
default_pda_imp = train(
  status ~ .,
  data = training,
  method = "pda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#sda
set.seed(200)
default_sda_imp = train(
  status ~ .,
  data = training,
  method = "sda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#slda
set.seed(200)
default_slda_imp = train(
  status ~ .,
  data = training,
  method = "slda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#hdda
set.seed(200)
default_hdda_imp = train(
  status ~ .,
  data = training,
  method = "hdda",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#pam
set.seed(200)
default_pam_imp = train(
  status ~ .,
  data = training,
  method = "pam",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#SIMCA
set.seed(200)
default_simca_imp = train(
  status ~ .,
  data = training,
  method = "CSimca",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#pls
set.seed(200)
default_pls_imp = train(
  status ~ .,
  data = training,
  method = "pls",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#glmnet
set.seed(200)
default_glmnet_imp = train(
  status ~ .,
  data = training,
  method = "glmnet",
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength=10
)
#nnet
set.seed(200)
default_nnet_imp = train(
  form = status ~ .,
  data = training,
  preProc=c("center", "scale"),
  trControl = trainControl(method = "repeatedcv", number = 10, repeats=5),
  method = "nnet",
  tuneLength=10
)
#gbm
set.seed(200)
default_gbm_imp <- train(
  status ~.,
  data = training, 
  method = "gbm",
  preProc=c("center", "scale"),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5),
  tuneLength = 100
)



# Saving dotplot featsel --------------------------------------------------

results_imp <- resamples(list(RPART=default_rpart_imp,NNET=default_nnet_imp,GLMNET=default_glmnet_imp,PLS=default_pls_imp,SIMCA=default_simca_imp,PAM=default_pam_imp,HDDA=default_hdda_imp,SLDA=default_slda_imp,SDA=default_sda_imp,PDA=default_pda_imp,KKNN=default_kknn_imp,SVM=default_svm_imp,GBM=default_gbm_imp,LVQ=default_lvq_imp,KNN=default_knn_imp,RF=default_rf_imp,GLM=default_glm_imp))
dotplot(results_imp, main='Feature selection')
dotplot_featuresel <- dotplot(results_imp, main='Feature selection')
trellis.device(device="png", filename="../Plots/Dotplot_featsel.png", dpi=300)
print(dotplot_featuresel)
dev.off()


# Get metrics from Rpart model --------------------------------------------
# First quickly rerun, to make sure we have standardized healthy and weak.
levels(training$status) <- c("Healthy","Weak")
#rpart
set.seed(200)
default_rpart_imp_levels <- train(
  status ~.,
  data = training, 
  method = "rpart",
  parms = list(split = "information"),
  #preProc=c("center", "scale"),
  tuneGrid = expand.grid(cp = seq(from = 0.001, to = 0.1, length = 200)),
  trControl=trainControl(method = "repeatedcv", number = 10, repeats=5, classProbs = TRUE, savePredictions = "final"),
  tuneLength = 10
)

fancyRpartPlot(default_rpart_imp_levels$finalModel, caption=NA)
plot(default_rpart_imp_levels)


preds <- predict(default_rpart_imp_levels, testing)

confusionMatrix(table(preds,testing$status))

