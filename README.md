# Introduction  
This R package implements estimation uncertainty using the exceedance probability criterion to ensure the desired in-control performance, based on extreme value theory.

# Install the SemEPC package from GitHub
1. Install the `devtools` package
```{r}
install.packages("devtools")
```
2. Load the `devtools` package
```{r}
library(devtools)
```
3. Install the `SemEPC` package
```{r}
install_github("chungili/SemEPC")
```
# Install required packages 
```{r}
install.packages("DescTools")
install.packages("dfphase1")
```

# Load the required packages
```{r}
library(SemEPC)
library(DescTools)
library(dfphase1)
```

# Example: ucisecom Dataset
```{r}
data("ucisecom")
# The 2nd variable is selected and then named as x1
x1 = prepro(ucisecom$V2)[1:61]
n = length(x1)
```

# Construct EPC control limits using Pickands' method
```{r}
Pickands(x1, pn=Pn_Pick(n))
```

# Construct EPC control limits using the Moment method
```{r}
Moment(x1, pn=Pn_Mom(n))
```
