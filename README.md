# A multivariate geostatistical framework for combining multiple indices of abundance for disease vectors and reservoirs: A case study of *rattiness* in a low-income urban Brazilian community

## Article
The manuscript has been submitted for publication. A link to a pre-print will be made available here soon.

## Data
### 1. Rat index data
File: [1-PdLRattinessData.csv](https://github.com/maxeyre/Rattiness-1/blob/master/Data/1-PdLRattinessData.csv) 

The cleaned dataset for the cross-sectional survey for rat signs, traps and track plates in Pau da Lima community, Salvador, Brazil, which was collected between October and December 2014. It is georeferenced and includes all relevant covariate values used in the model described in the article. A codebook for variables can be found [here]().

This dataset was collected by and belongs to the [Instituto de Saúde Coletiva, Universidade Federal da Bahia (ISC-UFBA)](http://www.isc.ufba.br/). If you use this dataset please could you cite ISC-UFBA and this published article.

### 2. Prediction grid
File: [2-PdLPredictionGrid.csv](https://github.com/maxeyre/Rattiness-1/blob/master/Data/2-PdLPredictionGrid.csv)

A 5m by 5m prediction grid with covariate values within the Pau da Lima study area.

## Software
The statistical software developed for this project and used for the *rattiness* analysis is divided into four main scripts and six additional scripts for fitting two-indices models and making predictions. Each script includes a description of its function within the R script. All scripts are contained [here](https://github.com/maxeyre/Rattiness-1/tree/master/Scripts).

### 1. Exploratory analysis - exploring covariate relationships
File: [1-ExploreCovariates.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/1-ExploreCovariates.R)

### 2. Exploratory analysis - Fit non-spatial models with covariates, LRT tests for each index, test for residual spatial correlation
File: [2-LRTTestSpatialCorr.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/2-LRTTestSpatialCorr.R)

### 3. Full model - Fit multivariate geostatistical model
File: [3-FitFullModel.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/3-FitFullModel.R)

### 4. Prediction - Make *rattiness* predictions at unobserved locations
File: [4-Prediction.R](https://github.com/maxeyre/Rattiness-1/blob/master/Scripts/4-Prediction.R)

### 5. Two-indices models - Fitting and prediction
Folder: [Two-indices-models](https://github.com/maxeyre/Rattiness-1/tree/master/Scripts/Two-indices-models)

[License](https://github.com/maxeyre/Rattiness-1/blob/master/LICENSE)

## Contact
If you have any questions about the project, data or software please contact max.eyre@lstmed.ac.uk.
