##### BRIERS-LOUW ET AL. 2024 #####

#Single-session spatial capture-recapture models: 2021 survey only 

##### PREPARE DATA #####

### PREP WORKSPACE
options(timeout = max(1000, getOption("timeout")))
install.packages("tidyverse", dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("rgdal", dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("secr", dependencies = TRUE, INSTALL_opts = '--no-lock')
install.packages("devtools", dependencies = TRUE, INSTALL_opts = '--no-lock')

require(tidyverse)
require(devtools)
require(rgdal)
require(secr)

setwd(".")

#READ IN DATA

#Read in traps. Coordinates are in UTM
traps= read.traps ("Trapfile2021_hyaenas_Final.txt", detector="proximity", 
                    covnames = c("HumanIndex","HabitatFineScale","PreyIndex"), binary.usage = T)

#Read in captures
captures <- read.table("Capturefile2021_hyaenas_Final.txt",
                       header = F, 
                       col.names = c("Session", "ID", "Occasion", "Detector"))

### MAKE CAPTURE HISTORY OBJECT
#combine captures and traps into capthist
hyaenas21 <- make.capthist(captures, 
                           traps,noccasions = 140)

summary(hyaenas21)

hyaenas21

#test for closure
closure.test(hyaenas21)

##### MAKE MASKS #####

### Suggest buffer to get minimum starting point for buffer
#Buffer should be > 4*sigma #RB = Relative Bias

suggest.buffer(hyaenas21,
               detectfn = 0,
               detectpar = list(g0 = 0.02, sigma = 5100),
               noccasions = 140,
               ncores = 4,
               RBtarget = 1e-5)

mask25.1 <- make.mask(traps, buffer = 25000, 
                      spacing = 1000, 
                      type ="trapbuffer", poly = NULL)
mask25.2 <- make.mask(traps, buffer = 25000, 
                      spacing = 2000, 
                      type ="trapbuffer", poly = NULL)

mask30.1 <- make.mask(traps, buffer = 30000, 
                      spacing = 1000, 
                      type ="trapbuffer", poly = NULL)
mask30.2 <- make.mask(traps, buffer = 30000, 
                      spacing = 2000, 
                      type ="trapbuffer", poly = NULL)

##### Fit preliminary models #####
# run models to check buffer width, buffer spacing, pmix values, & convergence issues

### Basic model ###
# 25km buffer with 1 km spacing and no finite mixture
m0.25.1 <- secr.fit(capthist = hyaenas21, 
                    model = list(D ~ 1,
                                 g0 ~ 1, 
                                 sigma ~ 1), 
                    mask = mask25.1, 
                    detectfn = "HN", #Half-normal
                    binomN = 1, #size determined by usage
                    details = list(fastproximity = TRUE), 
                    method = "Nelder-Mead", 
                    trace = F, 
                    ncores = 4, 
                    control = list(maxit = 19999))

### Check convergence ###
# check m0.25.1 converged properly
m0.25.1_b <- secr.fit(capthist = hyaenas21, 
                      model = list(D ~ 1,
                                   g0 ~ 1, 
                                   sigma ~ 1), 
                      mask = mask25.1, 
                      detectfn = "HN", #Half-normal
                      binomN = 1, #size determined by usage
                      details = list(fastproximity = TRUE), 
                      method = "Nelder-Mead", 
                      trace = F, 
                      start = m0.25.1,
                      ncores = 4, 
                      control = list(maxit = 19999))

# compare coefficient estimates side by side
cbind(coef(m0.25.1)[,1],
      coef(m0.25.1_b)[,1])

#### Check Masks ####
m0.30.1 <- secr.fit(capthist = hyaenas21, 
                    model = list(D ~ 1,
                                 g0 ~ 1, 
                                 sigma ~ 1), 
                    mask = mask30.1, 
                    detectfn = "HN", #Half-normal
                    binomN = 1, #size determined by usage
                    details = list(fastproximity = TRUE), 
                    method = "Nelder-Mead", 
                    trace = F, 
                    start = m0.25.1,
                    ncores = 4, 
                    control = list(maxit = 19999))

# compare coefficient estimates side by side
cbind(coef(m0.25.1)[,1],
      coef(m0.30.1)[,1])

# difference as percentage of relSE
(coef(m0.30.1)[1,1] - coef(m0.25.1)[1,1])/ coef(m0.25.1)[1,2] * 100

### Check coarser mask ###
m0.25.2 <- secr.fit(capthist = hyaenas21, 
                    model = list(D ~ 1,
                                 g0 ~ 1, 
                                 sigma ~ 1), 
                    mask = mask25.2, 
                    detectfn = "HN", #Half-normal
                    binomN = 1, #size determined by usage
                    details = list(fastproximity = TRUE), 
                    method = "Nelder-Mead", 
                    trace = F, 
                    start = m0.25.1,
                    ncores = 4, 
                    control = list(maxit = 19999))

# compare parameter estimates
cbind(coef(m0.25.1)[,1],
      coef(m0.25.2)[,1])

#### Check finite mixture ####
m0.h2 <- secr.fit(capthist = hyaenas21, 
                  model = list(D ~ 1,
                               g0 ~ h2, 
                               sigma ~ h2), 
                  mask = mask25.1, 
                  detectfn = "HN", #Half-normal
                  binomN = 1, #size determined by usage
                  details = list(fastproximity = TRUE), 
                  method = "Nelder-Mead", 
                  trace = T, 
                  start = m0.25.1,
                  ncores = 4, 
                  control = list(maxit = 19999))

predict(m0.h2)
AIC(m0.h2, m0.25.1)

### Explore pmix profile ###
# run pmixprofileLL and plot results
# use coarser mask to make it faster
pmvals <- seq(0.02,0.98,0.01)

profileLL <- pmixProfileLL(hyaenas21, 
                           model = list(g0~h2, sigma~h2), 
                           pmi = 6,
                           pmvals = pmvals, 
                           mask = mask25.2, 
                           CL = FALSE, 
                           details = list(fastproximity = TRUE), 
                           method = "Nelder-Mead",
                           trace = T,
                           ncores = 4, 
                           control = list(maxit = 19999))

par(mfrow = c(1,1))

plot(pmvals, profileLL, xlim = c(0,1),
     xlab = 'Fixed pmix', ylab = 'Profile log-likelihood')

pmix_tbl <- tibble(pmix = pmvals,
                   llk = profileLL) %>%
  arrange(desc(llk)) 

View(pmix_tbl)

### Rerun with a bit more resolution within the window of maximum interest
pmvals2 <- seq(0.85, 0.98, 0.01)

profileLL2 <- pmixProfileLL(hyaenas21, 
                            model = list(g0~h2, sigma~h2), 
                            pmi = 6,
                            pmvals = pmvals2, 
                            mask = mask25.2, 
                            CL = FALSE, 
                            details = list(fastproximity = TRUE), 
                            method = "Nelder-Mead",
                            trace = T,
                            ncores = 4, 
                            control = list(maxit = 19999))
par(mfrow = c(1,1))

plot(pmvals2, profileLL2, xlim = c(0,1),
     xlab = 'Fixed pmix', ylab = 'Profile log-likelihood')

pmix_tbl2 <- tibble(pmix = pmvals2,
                    llk = profileLL2) %>%
  arrange(desc(llk)) 

View(pmix_tbl2)

##### EXTRACT COVARIATE VALUES #####

###EXTRACT RASTER VALUES
#Load spatial packages
library(raster)
library(sf)
library(nngeo)
library(sp)
library(ggmap)

###Read in covariate layers
#CRS: WGS 84 UTM zone 36S

#Land cover tifs
grass <- raster("./secr - shapefiles/Land Cover/GrasslandUTM.tif")
plot(grass)
grass@crs

other <- raster("./secr - shapefiles/Land Cover/OtherUTM.tif")
other@crs

shrubs <- raster("./secr - shapefiles/Land Cover/ShrublandUTM.tif")
shrubs@crs

trees <- raster("./secr - shapefiles/Land Cover/Trees_coverUTM.tif")
trees@crs

#Community distance
community <- st_read("./secr - shapefiles/Community 2021/CommunityUTM.shp")
st_crs(community)

###extract covariate values
msk.covs <- list() #empty list to populate with spatial points of mask with covariate values

#convert mask (25km) to spatial points
#Tibbles: convert into dataframe before input into secr
msk_pts <- mask25.1 %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs =32736, remove = F)

#sample covariates
msk.covs <- msk_pts %>% mutate(Grass = raster::extract(grass, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Trees = raster::extract(trees, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Shrubs = raster::extract(shrubs, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Other = raster::extract(other, msk_pts, buffer=7000, fun="mean",na.rm =TRUE),
                               Community = simplify(st_nn(msk_pts, community, k=1, returnDist = T)[[2]]),
                               Comm_log = log(replace(Community, Community < 1, 1)))
msk.covs %>%
  write_csv("./outputs/msk.covs.csv")

###scale covariates
#data frame with continuous predictors scaled to have mean 0 and std dev 1

mask.covs_st <- 
  msk.covs %>% 
  mutate_at(vars(c("Grass", "Trees", "Shrubs", "Other", "Comm_log")), function(x) scale(x)[,1])

###check correlations

cor_df <- mask.covs_st %>% 
  dplyr::select(Grass, Trees, Shrubs, Other, Comm_log)
cor_df <- st_drop_geometry(cor_df)
cor_tbl <- round(cor(cor_df), digits = 2)
cor_tbl2 <- cor(cor_df, use="pairwise.complete.obs")
write.csv(cor_tbl2, "covariate_correlations.csv")

###add covariates to mask object
covariates(mask25.1) <- mask.covs_st%>%
  as.data.frame()
summary(covariates(mask25.1))

##### MODEL FITTING #####

### SPECIFY MODELS
#specify competing models

secr1=secr.fit(hyaenas21, model=list(D~1, g0~h2, sigma~h2),
                mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, 
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr2=secr.fit(hyaenas21, model=list(D~1, g0~(h2+HumanIndex), sigma~h2),
                mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr3=secr.fit(hyaenas21, model=list(D~1, g0~(h2+PreyIndex), sigma~h2),
                mask=mask25.1,detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr4=secr.fit(hyaenas21, model=list(D~1, g0~(h2+HabitatFineScale), sigma~h2),
                mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, start = secr1, 
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr5=secr.fit(hyaenas21, model=list(D~1, g0~(h2+HabitatFineScale+HumanIndex), sigma~h2),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr6=secr.fit(hyaenas21, model=list(D~Comm_log, g0~h2, sigma~h2),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr7=secr.fit(hyaenas21, model=list(D~Trees, g0~h2, sigma~h2),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr8=secr.fit(hyaenas21, model=list(D~Trees + Comm_log, g0~h2, sigma~h2),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr9=secr.fit(hyaenas21, model=list(D~Trees*Comm_log, g0~h2, sigma~h2),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr1,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

### EVALUATE RESULTS
secrmods <- secrlist(secr1, secr2, secr3, secr4, secr5, secr6, secr7, secr8, secr9)

#check convergence
print(secrmods$fit$convergence)

#check compatibility
AICcompatible(secrmods)

#check support
AIC(secrmods, criterion = "AICc", dmax = 7)

#save AIC table
AIC_tblL <- as_tibble(AIC(secrmods, criterion = "AICc", dmax = 7), rownames = "Hyp") %>% 
  mutate(model = str_sub(model, 1, -46),
         dAICc = round(dAICc, 2), Weight = round(AICcwt, 2)) %>% 
  dplyr::select(-detectfn, -logLik, -AICc) %>%
  write_csv("./outputs/AICc_table.csv")

saveRDS(secrmods, file = "./outputs/secrmodels..rds")
save(secrmods, file = "./outputs/secrmodels.RData")

print(secr2)

#GOF test of top model
secrM <- secr.test (secr2, fit=TRUE, nsim = 99)
secrM

##Explore whether the outlier (order of magnitude difference in detection) is masking the effects of covariates

#Read in captures without HYE_M001 (outlier)
captures_X <- read.table("Capturefile2021_hyaenas_Final_X.txt",
                       header = F, 
                       col.names = c("Session", "ID", "Occasion", "Detector"))

### MAKE 2ND CAPTURE HISTORY OBJECT
#combine captures and traps into capthist
hyaenas21_X <- make.capthist(captures_X, 
                            traps,noccasions = 140)

summary(hyaenas21_X)

hyaenas21_X

#rerun competing models

secr10=secr.fit(hyaenas21_X, model=list(D~1, g0~h2, sigma~h2),
               mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, 
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr11=secr.fit(hyaenas21_X, model=list(D~1, g0~(h2+HumanIndex), sigma~h2),
               mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr12=secr.fit(hyaenas21_X, model=list(D~1, g0~(h2+PreyIndex), sigma~h2),
               mask=mask25.1,detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr13=secr.fit(hyaenas21_X, model=list(D~1, g0~(h2+HabitatFineScale), sigma~h2),
               mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, start = secr10, 
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr14=secr.fit(hyaenas21_X, model=list(D~1, g0~(h2+HabitatFineScale+HumanIndex), sigma~h2),
               mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr15=secr.fit(hyaenas21_X, model=list(D~Comm_log, g0~h2, sigma~h2),
               mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr16=secr.fit(hyaenas21_X, model=list(D~Trees, g0~h2, sigma~h2),
               mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr17=secr.fit(hyaenas21_X, model=list(D~Trees + Comm_log, g0~h2, sigma~h2),
               mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

secr18=secr.fit(hyaenas21_X, model=list(D~Trees*Comm_log, g0~h2, sigma~h2),
               mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr10,
               method = "Nelder-Mead",
               details = list(fastproximity = TRUE), 
               control = list(maxit = 19999), hcov=)

### EVALUATE RESULTS
secrmodsX <- secrlist(secr10, secr11, secr12, secr13, secr14, secr15, secr16, secr17, secr18)

#check convergence
print(secrmodsX$fit$convergence)

#check compatibility
AICcompatible(secrmodsX)

#check support
AIC(secrmodsX, criterion = "AICc", dmax = 7)

#save AIC table
AIC_tblLX <- as_tibble(AIC(secrmodsX, criterion = "AICc", dmax = 7), rownames = "Hyp") %>% 
  mutate(model = str_sub(model, 1, -46),
         dAICc = round(dAICc, 2), Weight = round(AICcwt, 2)) %>% 
  dplyr::select(-detectfn, -logLik, -AICc) %>%
  write_csv("./outputs/AICc_tableX.csv")

saveRDS(secrmods, file = "./outputs/secrmodelsX..rds")
save(secrmods, file = "./outputs/secrmodelsX.RData")

print(secr11)

#GOF test of top model
secrMX <- secr.test (secr11, fit=TRUE, nsim = 99)
secrMX

##Explore 3-class mixture models with outlier included

secr19=secr.fit(hyaenas21, model=list(D~1, g0~h3, sigma~h3),
                mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr20=secr.fit(hyaenas21, model=list(D~1, g0~(h3+HumanIndex), sigma~h3),
                mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr21=secr.fit(hyaenas21, model=list(D~1, g0~(h3+PreyIndex), sigma~h3),
                mask=mask25.1,detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr22=secr.fit(hyaenas21, model=list(D~1, g0~(h3+HabitatFineScale), sigma~h3),
                mask=mask25.1, detectfn = "HN",binomN = 1, CL=FALSE, trace=FALSE, start = secr19, 
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr23=secr.fit(hyaenas21, model=list(D~1, g0~(h3+HabitatFineScale+HumanIndex), sigma~h3),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr24=secr.fit(hyaenas21, model=list(D~Comm_log, g0~h3, sigma~h3),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr25=secr.fit(hyaenas21, model=list(D~Trees, g0~h3, sigma~h3),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr26=secr.fit(hyaenas21, model=list(D~Trees + Comm_log, g0~h3, sigma~h3),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

secr27=secr.fit(hyaenas21, model=list(D~Trees*Comm_log, g0~h3, sigma~h3),
                mask=mask25.1, detectfn = "HN", binomN = 1, CL=FALSE, trace=FALSE, start = secr19,
                method = "Nelder-Mead",
                details = list(fastproximity = TRUE), 
                control = list(maxit = 19999), hcov=)

### EVALUATE RESULTS
secrmods3 <- secrlist(secr19, secr20, secr21, secr22, secr23, secr24, secr25, secr26, secr27)

#check convergence
print(secrmods3$fit$convergence)

#check compatibility
AICcompatible(secrmods3)

#check support
AIC(secrmods3, criterion = "AICc", dmax = 7)

#save AIC table
AIC_tblL3 <- as_tibble(AIC(secrmods3, criterion = "AICc", dmax = 7), rownames = "Hyp") %>% 
  mutate(model = str_sub(model, 1, -46),
         dAICc = round(dAICc, 2), Weight = round(AICcwt, 2)) %>% 
  dplyr::select(-detectfn, -logLik, -AICc) %>%
  write_csv("./outputs/AICc_table3.csv")

saveRDS(secrmods, file = "./outputs/secrmodels3..rds")
save(secrmods, file = "./outputs/secrmodels3.RData")

print(secr21)