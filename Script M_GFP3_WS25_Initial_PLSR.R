

##############################################
###  PLS base model ###
##############################################

install.packages("prospectr")
install.packages("resemble")
install.packages("caret")
install.packages("pls")
install.packages("randomForest")
install.packages("ranger")

#### Package administration: ####
library(prospectr)
library(resemble)
library(caret)
library(pls)
library(randomForest)
library(ranger)

#### Import of datasets: ####
setwd("C:/Users/camah/OneDrive/Dokumente/Master/WiSe25_26/Spectroscopy/dataset/dataset")

lab <- read.csv("lab.csv")

VNIR<-read.csv("VNIR.csv")
MIR<-read.csv("MIR.csv")
VNIR<-VNIR[,-1]
MIR<-MIR[,-1]

wl <- read.csv("wl.csv")$wl
wn <- read.csv("wn.csv")$wn

# Soil parent material
dev<-seq(1,45,1)
lia<-seq(46,90,1)
mus<-seq(91,135,1)
rot<-seq(136,180,1)

##################################################################################
### Spectral preprocessing: ###
###############################

# VNIR reflectance
nXr00<-as.matrix(VNIR) # reflectance without further preprocessing
nXr01<-sweep(nXr00, MARGIN=1, apply(nXr00, 1, function(x) (mean(x))), FUN="/")# mean centering
nXr02<-standardNormalVariate(nXr00) # SNV
nXr03<-detrend(nXr00, wl) # Detrend, requires additional information on wavelengths/wavenumbers
nXr04<-savitzkyGolay(nXr00, m=1, p=2, w=5) # 1st derivative (Savitzky-Golay)

# MIR reflectance
mXr00<-as.matrix(MIR)
mXr01<-sweep(mXr00, MARGIN=1, apply(mXr00, 1, function(x) (mean(x))), FUN="/")
mXr02<-standardNormalVariate(mXr00)
mXr03<-detrend(mXr00, wn)
mXr04<-savitzkyGolay(mXr00, m=1, p=2, w=15)


###################################################

# optional: outlier analysis - not necessary for our datasets

####################################################

# Calibration sampling (e.g. k-Means, Kennard-Stone, cLHS)

# In our session from 12.12.24 we decided to define our 'independent' test set using 
# the Kennard-Stone algorithm with the 40 "most extreme" samples to be the test set for 
# reasons of robustness (models have to extrapolate a little bit, which simulates a more 
# realistic scenario of new, unknown data)

# definition of calibration and test set based on raw VNIR reflectance (without any preprocessing)
nXdat<-nXr00

# size of test set:
nt<-40

# Kennard-Stone algorithm
split<-kenStone(X = nXdat, 
                k = nt, 
                metric="mahal",
                pc=0.95,   # based on PC feature space explaining 95% of overall variance of the spectra
                .center = TRUE,
                .scale = FALSE)


# which samples are defined as test set
test<-split$model[-21]   # removal of an outlier

# which samples are defined as calibration set
train<-split$test

# Plot in PC feature space
pc<-prcomp(nXdat, center=T, scale. = F, rank. = 5)

par(pty="s")      # pty="s" ---> squared plot, regardless of the plot window size
matplot(pc$x[,1],   # scores of first PC
        pc$x[,2],   # scores of second PC
        pch=16)     # filled points

# adding test set
matpoints(pc$x[test,1], # only test set samples, 1st PC
          pc$x[test,2], # only test set samples, 2nd PC
          pch=16, 
          col="red")

# show the outlier
matpoints(pc$x[21,1], # outlier, 1st PC
          pc$x[21,2], # outlier, 2nd PC
          pch=16, 
          col="blue")

legend("topright", c("model", "test"), pch=16, col=c("black", "red"))

##############################################################################################

###########################
### baseline PLS model: ###
###########################

# in the first step, we calculated a baseline predicion model using PLS regression 
# based on raw VNIR reflectance

# target variable: amount of clay minerals (in %)

# plsr() function requires input data as matrix
nXdat<-nXr00

Ytrain <- as.matrix(lab$Clay[train])   # calibration set
Xtrain <- as.matrix(nXdat[train,])     # calibration set

# Model without cross validation
modf <- plsr(Ytrain~Xtrain,       # variables
             ncomp=20,  # max. number of latent variables to compute
             method="oscorespls",    # oscorespls = NIPALS = standard PLS
             validation="none")      # without any validation
summary(modf)

# Model with 10-fold cross validation
set.seed(123)
modcv <- plsr(Ytrain~Xtrain,        # variables
              ncomp=20,   # max. number of latent variables to compute
              method="oscorespls",     # oscorespls = NIPALS 
              validation="CV",   # with n-fold cross-validation
              segments=10)       # n of segments = 10 --> 10-fold CV
summary(modcv)

# RMSE of model fit
err_fit<-RMSEP(modf)

# RMSE of cross validation
err_modcv<-RMSEP(modcv)

# plot them together
plot(err_fit)       # RMSEcv of model fit
plot(err_modcv, add=T)# RMSEcv of cross validation - find optimal number of latent variables
# with increasing number of latent variables, PLS is able to explain an increasing amount of 
# variance of Y

# with cross validation, we can define the point where the model is overfitting, i.e. where the model
# uses sample specific noise to explain Y, which results in an increase of RMSE in the test set


################
# set optimal number of latent variables here, defined as minimum RMSE:
onc<-which.min(err_modcv$val[1,,])-1      #  - 1 to compensate to remove the RMSE of the intercept

# Inspection of regression coefficients
par(pty="m")          # pty="m" --> maximizing plot to plot window
matplot(wl, modcv$coefficients[,,onc], type="l", lty=1,
        xlab="wavelengths (nm)", ylab="PLS coefficients")
axis(1, seq(500,2500,100), labels = NA)
abline(h=0, lty=2, col="gray")

# we clearly see the most important wavelengths e.g. around 1420, 1905, 2220 nm,
# which are caused by water and hydrogen absorption and can be indirectly linked 
# to clay mineral absorption in soil samples

#######
# evaluation of cross validated results

# scatter plot of observed vs predicted values
par(pty="s")
matplot(Ytrain,
        modcv$validation$pred[,,onc],
        xlab="observed values",
        ylab="predicted values",
        pch=16)
# adding a 1:1 line
abline(0,1, lty=2, col="darkgrey")

# validation metrics:
Yobs<-Ytrain         # observed model values
np<-nrow(Yobs)     # number of values
Ycv <- modcv$validation$pred[,,onc]


# calculation of R²:
SST<-sum((Yobs-mean(Yobs))^2)        
SSE<-sum((Ycv-Yobs)^2)
R2<-round(1-(SSE/SST), digits=2)

# calculation of RMSE
RMSE<-round(sqrt(sum((Ycv-Yobs)^2)/np), digits=2)

# Calculation of RPD
RPD<-round(sd(Yobs)/RMSE, digits=2)

print(c(R2, RMSE, RPD))

# Text für die Legende
legend_text <- c(paste("R²_cv =", R2),
                 paste("RMSE_cv =", RMSE),
                 paste("RPD_cv =", RPD))

# Legende oben links hinzufügen
legend("topleft", legend=legend_text, bty="n", cex=1.2)

#################################################################################
###############################
### Prediction on test set: ###
###############################

Ynew <- as.matrix(lab$Clay[test])  # test set
Xnew <- as.matrix(nXdat[test,])    # test set

# actual prediction:
Ypred <- predict(object=modcv,    # prediction based on model
                 ncomp=onc,       # optimal number of latent variables
                 newdata = Xnew)  # new, 'unknown' data

# Scatter plot:
par(pty="s")
matplot(Ynew, 
        Ypred[,1,1], 
        pch=16)
abline(0,1, lty=2)   # 1:1 line

# validation metrics:
np<-nrow(Ynew)     # number of predictions

# calculation of R²:
SST<-sum((Ynew-mean(Ynew))^2)        
SSE<-sum((Ynew-Ypred[,1,1])^2)
R2<-round(1-(SSE/SST), digits=2)
R2

# calculation of RMSE
RMSE<-round(sqrt(sum((Ypred[,1,1]-Ynew)^2)/np), digits=2)
RMSE

# Calculation of RPD
RPD<-round(sd(Ynew)/RMSE, digits=2)
RPD

print(c(R2, RMSE, RPD))

# text for the legend
legend_text <- c(paste("R² =", R2),
                 paste("RMSE =", RMSE),
                 paste("RPD =", RPD))

# legend
legend("topleft", legend=legend_text, bty="n", cex=1.2)


#################################################################################

#################################################################################
#### Alternative: train() function in the caret package: ####
#############################################################

# In the following, we are using the train() function of the caret package to repeat 
# the PLS model
# The advantage is that the train() function provides a unified and systematic 
# framework for model training, validation, and optimization, including multiple 
# regression and classification methods (e.g. linear regression, random forest, SVM, etc.).

# plsr() function requires input data as matrix
nXdat<-nXr00

# dependent variable and predictors
Ytrain <- as.matrix(lab$Clay[train])   # training set
Xtrain <- as.matrix(nXdat[train,])     # training set 

# combine to a data frame
traindat <- data.frame(target = Ytrain, Xtrain)

# Set up a control object for cross-validation (part of the train function)
set.seed(123)
control <- trainControl(method = "cv", number = 10, savePredictions = "final")

# Set up a tuning object for grid search of optimal model parameters (part of the train function)
tuneGrid = expand.grid(ncomp = 1:20) # Specify a grid of latent variables to tune

# Partial Least Squares Regression (PLSR)
set.seed(123)
plsr_model <- train(target~ .,
                    data = traindat,
                    method = "pls",
                    preProcess = c("center"),
                    trControl = control,
                    tuneGrid = tuneGrid)

print(plsr_model)
summary(plsr_model)

# variable importance:
help(varImp)

# Partial Least Squares: the variable importance measure here is based on 
# weighted sums of the absolute regression coefficients. 
# The weights are a function of the reduction of the sums of squares across 
# the number of PLS components and are computed separately for each outcome. 
# Therefore, the contribution of the coefficients are weighted proportionally 
# to the reduction in the sums of squares.

VI<-varImp(plsr_model)
matplot(wl, VI$importance, type="l", lty=1)
axis(1, seq(500,2500,100), labels = NA)

# regression coefficients from train object:
par(pty="m")          # pty="m" --> maximizing plot to plot window
matplot(wl, coef(plsr_model$finalModel, ncomp=plsr_model$bestTune$ncomp), 
        type="l", lty=1,
        xlab="wavelengths (nm)", ylab="PLS coefficients")
axis(1, seq(500,2500,100), labels = NA)
abline(h=0, lty=2, col="gray")


#############################################
# prediction of test set:

Ynew <- as.matrix(lab$Clay[test])  # test set
Xnew <- as.matrix(nXdat[test,])    # test set

# this way, the optimal model (best fit) is used for prediction
Ypred<- predict(plsr_model,
                newdata=Xnew)

# adjusting the optimal model is possible:

# Apply same preprocessing to the new data
#Xnew_prep <- predict(plsr_model$preProcess, newdata=Xnew)

# specifiying the required model configuration:
# Ypred <- predict(object=plsr_model$finalModel,    # prediction based on model
#                  ncomp=12,
#                  newdata  = Xnew_prep)  # new, 'unknown' data

# Scatter plot:
par(pty="s")
matplot(Ynew, 
        Ypred, 
        pch=16)
abline(0,1, lty=2)   # 1:1 line

# validation metrics:
np<-nrow(Ynew)     # number of predictions

# calculation of R²:
SST<-sum((Ynew-mean(Ynew))^2)        
SSE<-sum((Ynew-Ypred)^2)
R2<-round(1-(SSE/SST), digits=2)
R2

# calculation of RMSE
RMSE<-round(sqrt(sum((Ypred-Ynew)^2)/np), digits=2)
RMSE

# Calculation of RPD
RPD<-round(sd(Ynew)/RMSE, digits=2)
RPD

print(c(R2, RMSE, RPD))

# text for the legend
legend_text <- c(paste("R² =", R2),
                 paste("RMSE =", RMSE),
                 paste("RPD =", RPD))

# legend
legend("topleft", legend=legend_text, bty="n", cex=1.2)



#################################################################################
#################################################
#### Exercise: Optimizing the Baseline Model ####
#################################################

# The aim of this exercise is to optimize the baseline PLS model by exploring, e.g., 
# alternative preprocessing techniques, spectral ranges and modeling algorithms. 
# You will evaluate the optimized models using key performance metrics (R², RMSE, RPD)
# and compare them to the baseline model.

# Properties of the baseline model were:
# spectral range: VNIR
# spectral preprocessing: none
# regression algorithm: PLS regression
# training set: n=140
# test set: n=40 (defined by Kennard-Stone algorithm)

# Model quality metrics of the baseline model applied to the test set were:
# R² = 0.85
# RMSE = 3.59
# RPD = 2.64

# Your task is to improve the baseline model, e.g. by performing one of the following:
# 1) Modify the preprocessing steps
# 2) Change the spectral range (VNIR vs. MIR)
# 3) Change the model algorithm

# After  each change in preprocessing or modeling:

# 1) Train the new model on the training dataset.
# 2) Apply the trained model to the independent test dataset.
# 3) Calculate and record the performance metrics (R², RMSE, RPD)
# 4) Compare these metrics with the baseline model to determine whether 
#    your modifications improved the model.



# Model with PCA instead of PLSR
modf_pca <- pcr(
  Ytrain ~ Xtrain,
  ncomp = 20,            # number of principal components
  method = "svdpc",      # PCA via SVD (default)
  validation = "none"    # no cross-validation
)

summary(modf_pca)


# RMSE of model fit
err_fit_pca<-RMSEP(modf_pca)



# plot them together
plot(err_fit_pca)       # RMSEcv of model fit

