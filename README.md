# bayesMRM

R package bayesMRM for Bayesian Multivariate Receptor Modeling 

## How to install (last update: Nov 04, 2020)

### Note:

You may get a message about updating the packages needed by bayesMRM during installation, and updating the packages may cause an error.  
Restart R or Rstudio and skip updating will resolve the problem. 

### Installation steps

Step 1. Install the latest version of R (currently 4.0.3) and Rtools from https://cran.r-project.org

Step 2. Install the latest version of JAGS(currently4.3.0) from http://mcmc-jags.sourceforge.net.

Step 3. Follow either Step 3A or Step 3B


   * Step 3A (download from github): In R, execute 
    
```{r}
   install.packages("devtools")  #this is needed only once in R
   library(devtools)
   install_github("mansukoh/bayesMRM",dependencies=TRUE)
   library(bayesMRM)
```


   * Step 3B (use tar.gz file): Download bayesMRM_currentversion.tar.gz(currently, bayesMRM_0.1.1.tar.gz) from 
https://github.com/mansukoh/Rpackages and save it in a directioy (e.g.,  "c:/temp").  In R, execute 
 
```{r}    
   install.packages(c("rjags","coda","ggplot2","tibble","grDevices","graphics",
                       "gridExtra","rgl","shiny","shinythemes","crayon"))
   install.packages("c:/temp/bayesMRM_0.1.1.tar.gz",repos=NULL,type="source")
   library(bayesMRM)
```

