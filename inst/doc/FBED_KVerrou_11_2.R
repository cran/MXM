## ---- warning =  FALSE, message = FALSE------------------------------------

### ~ ~ ~ Load Packages ~ ~ ~ ###
library(MXM) 
library(dplyr)

## --------------------------------------------------------------------------

### ~ ~ ~ Load The Dataset ~ ~ ~ ###
wine.url <- "ftp://ftp.ics.uci.edu/pub/machine-learning-databases/wine/wine.data"
wine <- read.csv(wine.url,
                 check.names = FALSE,
                 header = FALSE) 
head(wine)
str(wine)
colnames(wine) <- c('Type', 'Alcohol', 'Malic', 'Ash', 
                    'Alcalinity', 'Magnesium', 'Phenols', 
                    'Flavanoids', 'Nonflavanoids',
                    'Proanthocyanins', 'Color', 'Hue', 
                    'Dilution', 'Proline')

## --------------------------------------------------------------------------
### ~ ~ ~ Removing The Categorical ('Type') and The Target ('Nonflavanoids') Variables ~ ~ ~ ###

wine_dataset <- dplyr::select(wine,
                              -contains("Type"),
                              -contains("Nonflavanoids")) 
head(wine_dataset)

## --------------------------------------------------------------------------
wine_target <- wine$Nonflavanoids
head(wine_target)

## --------------------------------------------------------------------------
### ~ ~ ~ Running FBED For First Time ~ ~ ~ ###
fbed_default_1st <- MXM::fbed.reg(target       = wine_target,
                                 dataset        = wine_dataset, 
                                 test     = "testIndFisher", 
                                 threshold    = 0.05,
                                 wei      = NULL,
                                 K        = 10,
                                 method   = "eBIC",
                                 gam      = NULL,
                                 backward = TRUE)

## --------------------------------------------------------------------------
fbed_default_1st$res
SelectedVars_names<-colnames(wine_dataset[fbed_default_1st$res[,1]])
SelectedVars_names

## --------------------------------------------------------------------------
fbed_default_1st$res[,2]

## --------------------------------------------------------------------------
fbed_default_1st$info

## --------------------------------------------------------------------------
fbed_default_1st$back.rem 

## --------------------------------------------------------------------------
fbed_default_1st$back.n.tests 

## --------------------------------------------------------------------------
fbed_default_1st$runtime 

## --------------------------------------------------------------------------
### ~ ~ ~ Taking The Whole Dataset ~ ~ ~ ###
wine_dataset <- wine
head(wine_dataset)

## ---- message=FALSE--------------------------------------------------------
### ~ ~ ~ Running FBED For Categorical Variable ~ ~ ~ ###
wine[, 1] <- as.factor(wine[, 1])
fbed_default_2nd <- MXM::fbed.reg(target       = wine_target,
                                 dataset        = wine_dataset, 
                                 test     = "testIndReg", 
                                 threshold    = 0.05,
                                 wei      = NULL,
                                 K        = 10,
                                 method   = "eBIC",
                                 gam      = NULL,
                                 backward = TRUE) 


## --------------------------------------------------------------------------
fbed_default_2nd$res
SelectedVars_names<-colnames(wine_dataset[fbed_default_2nd$res[,1]])
SelectedVars_names

## --------------------------------------------------------------------------
fbed_default_2nd$res[,2]

## --------------------------------------------------------------------------
fbed_default_2nd$info

## --------------------------------------------------------------------------
fbed_default_2nd$back.rem 
fbed_default_2nd$back.n.tests 

## --------------------------------------------------------------------------
fbed_default_2nd$runtime 

## --------------------------------------------------------------------------
sessionInfo()

