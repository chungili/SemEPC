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

# Example: ucisemcom Dataset
```{r}
data("ucisecom")
x = prepro(ucisecom$V2)
n = length(x)
```

# Construct EPC control limits using Pickands' method
```{r}
Pickands(x, pn=Pn_Pick(n))
```

# Construct EPC control limits using the Moment method
```{r}
Moment(x, pn=Pn_Pick(n))
```
