#loading required libraries
library(stabs)
library(factoextra)
library(NbClust)
library(FunCluster)
library(ggfortify)
library(glmnet)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(randomForest)
library(data.table)
library(mlr)
library(h2o)
library(caret)
library(plsVarSel)
library(pROC)
library(jtools)
library(dplyr)
library(cluster)
library("flexclust")
library(corrplot)
library(tsne)
library(clusterCrit)

setwd("C:/Users/gbeng/OneDrive/Documents")

#Importing data
gene_proteins <- read.csv("PAM50_proteins.csv")
clinical <- read.csv("clinical_data_breast_cancer.csv")
proteomes <- read.csv("77_cancer_proteomes_CPTAC_itraq.csv")


#transposing proteome matrix to make observations into rows rather than columns

#save rownames
n <- proteomes$RefSeq_accession_number

#transpose all but the first 3 column 
proteomes <- as.data.frame(t(proteomes[,4:86]))
colnames(proteomes) <- n

#rownames to first column
proteomes <- cbind(rownames(proteomes), data.frame(proteomes, row.names=NULL))
colnames(proteomes)[1] <- "Complete.TCGA.ID"

#reformatting Complete.TCGA.ID as clinical format to allow joining of data sets

#defining formula to restructure:
get.clinical.id <- function(proteome.id) {
  x = substr(proteome.id, 4, 7)
  y = substr(proteome.id, 0, 2)
  paste("TCGA",y,x,sep="-")
}

#sapply to id column in proteomes
proteomes$Complete.TCGA.ID <- sapply(proteomes$Complete.TCGA.ID, get.clinical.id)
proteomes_all <- proteomes

#looking for proteomes with many NAs
naCounts <- colSums(is.na(proteomes)) / nrow(proteomes)

#plotting missing data proportions

plot(sort(naCounts, decreasing = TRUE), col ="red", type = 'h', xlab = "index of proteome", ylab="proportion of missing data", main = "Propotion of missing data for each proteome") 

#how many have more than 25% missing data
length(naCounts[naCounts>0.25])

#remove variables with >25% missing data
proteomes <- proteomes[ , colSums(is.na(proteomes))  / nrow(proteomes) < 0.25] 

#loop to impute means for remaining missing data
for (i in which(sapply(proteomes, is.numeric))) {
  proteomes[is.na(proteomes[, i]), i] <- mean(proteomes[, i],  na.rm = TRUE)
}

#inner join on data to create full data set
library(dplyr)

data <-  inner_join(clinical, proteomes, by = "Complete.TCGA.ID")
## Warning: Column `Complete.TCGA.ID` joining factor and character vector,
## coercing into character vector
#replacing lengthy col name
colnames(data)[3] <- "diag_age"

#Barplot of subtypes
library(ggplot2)

ggplot(data, aes(PAM50.mRNA, col = PAM50.mRNA, fill = PAM50.mRNA, alpha=0.7)) + geom_bar() + ggtitle("Proportion of patients with each cancer subtype") 


#creating test/train split index
set.seed(1)
library(caret)
samp <- createDataPartition(data$PAM50.mRNA, p = 0.7, list = FALSE)

options(warn = -1) #turn off warnings

#Stability analyses
library(glmnet)

#creating a function to repeat lasso regression and return the selected model variables
LassoSub=function(k=1, Xdata, Ydata){
  set.seed(k)
  s=sample(nrow(data), size=0.8*nrow(data))
  Xsub=Xdata[s, ]
  Ysub=Ydata[s]
  model.sub=cv.glmnet(x=Xsub, y=Ysub, alpha=1, family="multinomial") #cross validated lasso
  coef.sub=coef(model.sub, s='lambda.1se')[-1] #using lambda +1se hyperparameter value for parsimony
  return(coef.sub)
}

#Run model 100 times and save results
niter=100
lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=as.matrix(data[,31:ncol(data)]), Ydata=as.matrix(data[,21]))

#create a matrix of all predictor variables
stability_matrix <- matrix(nrow=length(lasso.stab[[1]]),ncol=length(lasso.stab))
rownames(stability_matrix) <- rownames(lasso.stab[[1]])

#loop through to put list contents into matrix
for (i in 1:300){
  temp.data.frame <- as.matrix(lasso.stab[[i]])
  stability_matrix[,i] <- temp.data.frame
}

stability_matrix <- ifelse(stability_matrix != 0, 1, 0) #Replacing beta values with binary 1/0 (selected/not selected)
stability_matrix <- stability_matrix[2:nrow(stability_matrix),] #remove intercept value
stable_variables <- as.data.frame(rowSums(stability_matrix)) #create data frame with count of how many times each variable is selected for a model
stable_variables$protein <- rownames(stable_variables) #create column of variable names

colnames(stable_variables)[1] <- "times_selected" #assign appropriate column name
stable_variables <- stable_variables[!is.na(stable_variables$times_selected),]  #remove NAs
stable_variables <- stable_variables[stable_variables$times_selected != 0,] #remove all variables that were never selected

stable_variables <- stable_variables[order(-stable_variables$times_selected),] #ordering by number of times selected


#plotting stable variables
ggplot(stable_variables[1:30,], aes(x=reorder(as.factor(protein),-abs(times_selected),mean), y=times_selected, col =reorder(as.factor(protein),-abs(times_selected),mean), fill =reorder(as.factor(protein),-abs(times_selected),mean))) + geom_col(show.legend = FALSE, alpha = 0.6) + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=10)) + xlab("Protein") + ylab("Times selected") 

stable_variables <- stable_variables[1:30,]

STABVARS <- stable_variables$protein

STABVARS.ind <- which(colnames(data) %in% STABVARS)
library(gridExtra)

#for (i in 1:length(STABVARS[1:5])){
#  print(ggplot(data, aes_string("PAM50.mRNA", STABVARS[i], col="PAM50.mRNA", fill="PAM50.mRNA")) + geom_boxplot(alpha=0.3) + ggtitle(STABVARS[i]))
#}

set.seed(1)

#Settng up train control for cross validation and calibration of hyperparameter
train_control <- trainControl(method="repeatedcv", number=3, repeats=10, savePredictions = TRUE, summaryFunction = multiClassSummary) 

#Tunegrid for different values of C
grid <- expand.grid(C = seq(0.000001,0.15,0.002))

set.seed(1)
#Model training
svm.lin.mod <- train(PAM50.mRNA ~ ., data=data[samp, c(21, STABVARS.ind)], trControl=train_control, method="svmLinear", preProcess = c("center","scale"), tuneGrid =grid, tuneLength = 10)

#Creating predictions on test set
svm.predicts <- predict(svm.lin.mod, newdata = data[-samp, c(21, STABVARS.ind)])

#viewing confusion matrix
cm_svm <- confusionMatrix(svm.predicts, as.factor(data$PAM50.mRNA[-samp]))


set.seed(1)
#Model training
knn.lin.mod <- train(PAM50.mRNA ~ ., data=data[samp, c(21, STABVARS.ind)], trControl=train_control, method="knn", preProcess = c("center","scale"), tuneLength = 10)


#Creating predictions on test set
knn.predicts <- predict(knn.lin.mod, newdata = data[-samp, c(21, STABVARS.ind)])

#viewing confusion matrix
cm_knn <- confusionMatrix(knn.predicts, as.factor(data$PAM50.mRNA[-samp]))



set.seed(1)
#Model training
rf.lin.mod <- train(PAM50.mRNA ~ ., data=data[samp, c(21, STABVARS.ind)], trControl=train_control, method="rf", preProcess = c("center","scale"), tuneLength = 10)


#Creating predictions on test set
rf.predicts <- predict(rf.lin.mod, newdata = data[-samp, c(21, STABVARS.ind)])

#viewing confusion matrix
cm_rf <- confusionMatrix(rf.predicts, as.factor(data$PAM50.mRNA[-samp]))

library(yardstick)
library(ggplot2)

autoplot(cm_rf$table, type = "heatmap") +
           scale_fill_gradient(low="#D6EAF8",high = "#2E86C1")

library(ggfortify)
df <- data
colnames_remove <- colnames(data[1:30])
df <- df[,!names(df) %in% colnames_remove]
pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res, data = data, colour = 'PAM50.mRNA', size=5) #PCA plot colored by breast cancer subtypes
autoplot(pca_res, data = data, colour = 'Stage', size=5) #PCA plot colored by breast cancer subtypes

chisq.test(data$Node, data$PAM50.mRNA)

ggplot(data, aes(Stage, col = Stage, fill = AJCC.Stage, alpha=0.7)) + geom_bar() + ggtitle("Distribution of patients with each cancer stage") 


 library(dplyr)
 data <- data %>% mutate(AJCC.Stage=recode(AJCC.Stage, 
                                                                 'Stage IIA'="2",
                                                                 'Stage IIB'="2",
                                                                 'Stage II'="2",
                                                                 'Stage I'="1",
                                                                 'Stage IA'="1",
                                                                 'Stage IB'="1",
                                                                 'Stage IIIA'="3",
                                                                 'Stage IIIB'="3",
                                                                 'Stage III'="3",
                                                                 'Stage IIIC'="3",
                                                                 'Stage IV'="4"))
 
 

 