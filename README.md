# bayesMRM

R package bayesMRM for Bayesian Multivariate Receptor Modeling 

# How to install bayesMRM (updated Feb 15, 2021)

Note: Currently, it seems to be in a transition period from R version 3.x.x to
version 4.x.x. We find out that updating the exisiting packages needed for our package causes an error when installing our package.  We recommend to skip updates until we give a new guideline about installation of the package 'bayesMRM'. 


## install the required packages

Step 1. Install the latest version of R (currently 4.0.3) and Rtools from https://cran.r-project.org

Step 2. Install Rstudio from https://cran.rstudio.com. (not required, but recommanded)

Step 3. Install the latest version of JAGS(currently4.3.0) from http://mcmc-jags.sourceforge.net.

Step 4. Follow either Step 4A or Step 4B

  * Step 4A (download from github)
   
   **In Rstudio, execute the commands: 

```{}
   install.packages("devtools")  #this is needed only once in R
   library(devtools)
   ## Select 'None' or click [return] to 'updating packages' message when running  'install_github'. 
   install_github("mansukoh/bayesMRM",dependencies=TRUE)
   library(bayesMRM)
```
**

 
 * Step 4B (use tar.gz file)

  ** Step 4B.1.    Download bayesMRM_currentversion.tar.gz( currently, bayesMRM_0.1.1.tar.gz) from
https://github.com/mansukoh/Rpackages and save this file in the directioy  "C:/TEMP"
 (or any directory you want). 
**(we need to save tar.gz and zip file for bayesMRM in github.com/mansukoh/Rpackages)**

  Step 4B.2. In Rstudio, execute the commands: 

```
   # select 'No' to the question of updates when running the command below. 
   install.packages(c("rjags","coda","ggplot2","tibble","grDevices","graphics",
                       "gridExtra","rgl","shiny","shinythemes","crayon"))
   install.packages("C://temp/bayesMRM_0.1.1.tar.gz",repos=NULL,type="source")
   library(bayesMRM)
```
**

