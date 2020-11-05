# bayesMRM

R package bayesMRM for Bayesian Multivariate Receptor Modeling 

# How to install bayesMRM (updated Nov 04, 2020)

## Note:

Updating the exisiting packages needed for our package may cause an error 
when installing bayesMRM.  In that case, restart Rstudio and select "no update" to the message of update options.  

## Steps

Step 1. Install the latest version of R (currently 4.0.3) and Rtools from https://cran.r-project.org

Step 2. Install the latest version of JAGS(currently4.3.0) from http://mcmc-jags.sourceforge.net.

Step 3. Follow either Step 3A or Step 3B

```{}
 Step 3A (download from github)
   
 In R, execute 

   install.packages("devtools")  #this is needed only once in R
   library(devtools)
   install_github("mansukoh/bayesMRM",dependencies=TRUE)
   library(bayesMRM)
```

 ```{}
  Step 3B (use tar.gz file)

  (1)    Download bayesMRM_currentversion.tar.gz( currently, bayesMRM_0.1.1.tar.gz) from
https://github.com/mansukoh/Rpackages and save this file in the directioy  "C:/temp"
 (or any directory you want). 

  (2) In R, execute 
    
   install.packages(c("rjags","coda","ggplot2","tibble","grDevices","graphics",
                       "gridExtra","rgl","shiny","shinythemes","crayon"))
   install.packages("C:/temp/bayesMRM_0.1.1.tar.gz",repos=NULL,type="source")
   library(bayesMRM)
```

