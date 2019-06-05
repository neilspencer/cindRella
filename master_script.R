library(cindRella)
library(plyr)
library(matrixStats)
library(MASS)
library(ggplot2)
library(RColorBrewer)

# show some contact surfaces and accidentals
# from the example data in the package
plotContactSurface(1, show_accidentals = T)
plotContactSurface(2, show_accidentals = T)
plotContactSurface(3, show_accidentals = T)
plotContactSurface(4, show_accidentals = T)
plotContactSurface(5, show_accidentals = T)
plotContactSurface(6, show_accidentals = T)

# the same 6 contact surfaces repeat 5 times with different accidentals

source("data_prep.R") # construct necessary objects to fit models
source("data_run.R") # fit the model, competitor models, and variants
source("data_evaluate.R") # evaluate posterior predictive on held out data
source("data_plot.R") # create plots demonstrating the results

#comparing performance of various competitor models on held out data
competitorsplot

# comparing the performance of variant models on held out data
variantsplot

# an illustration of the posterior predictive distribution for a held out shoe
predictplot


