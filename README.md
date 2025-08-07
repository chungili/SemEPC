# Introduction  
This R package implements estimation uncertainty using the exceedance probability criterion to ensure the desired in-control performance, based on extreme value theory.

# Install the SemEPC package from GitHub
1. Install the `devtools` package
install.packages("devtools")
2. Load the `devtools` package
library(devtools)
3. Install the `SemEPC` package
install_github("chungili/SemEPC")

# Install required packages 
install.packages("DescTools")
install.packages("dfphase1")

# Load the required packages
library(SemEPC)
library(DescTools)
library(dfphase1)

# Example: ucisemcom Dataset
data("ucisemcom")
x = prepro(ucisecom$V2)
n = length(x)

# Construct EPC control limits using Pickands' method
Pickands(x, pn=Pn_Pick(n))

# Construct EPC control limits using the Moment method
Moment(x, pn=Pn_Pick(n))
