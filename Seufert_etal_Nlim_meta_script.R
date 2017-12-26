################################################################################
# Meta-analysis: effects of N limitation on crop photosynthesis and N allocation
# R code to reproduce results in Seufert V, Granath G, Muller C.
#
# Contact (R-code): gustaf.granath@gmail.com
################################################################################

# Statistical analyses
library(metafor) # tested with version 2.0-0

# Plotting and data management
library(ggplot2)
library(gridExtra)
library(grid)
library(readxl)
library(reshape2)


# Load functions and data ####

# function to calculate effect sizes (log response ratio) and its variance 
# within study and covariances
calc_resp <- function (obj=NULL) {
  var.rr <- ((obj$SDe^2/(obj$Ne*(obj$Xe^2))) +
               (obj$SDc^2/(obj$Nc*(obj$Xc^2))))
  var.control <- (obj$SDc^2/(obj$Nc*(obj$Xc^2)))
  rr <- log(obj$Xe/obj$Xc)
  return(cbind(obj, data.frame(var.rr, var.control, rr)))
}

# function to get variance-covariance matrix of within-study sampling dependence
# blockList=TRUE if you want within-study varcovar blocks as a list.
v_func <- function (dat, ind.study="ind.study", 
                    var.control="var.control", 
                    var.rr="var.rr",
                    blockList=FALSE) {
  dat.agg <- data.frame(st = unique(dat[[ind.study]]), 
                        len = rle(as.numeric(dat[[ind.study]]))$lengths, 
                        var = unique(dat[[var.control]]))
  
  require(Matrix)  
  ll = list()
  ll.diag = numeric()
  for (i in 1:NROW(dat.agg)) {  
    ll.diag <- dat[dat[[ind.study]] == dat.agg$st[i], var.rr]
    block <- matrix(dat.agg$var[i], ncol = dat.agg$len[i], nrow = dat.agg$len[i]) 
    diag(block) <- ll.diag 
    ll[[i]] <- block 
  }
  mat <- as.matrix(bdiag(lapply(ll, function (x) x) ))
  #diag(mat) <- dat[[var.rr]]
  if(blockList) {return(ll)} else 
    return(mat)
}

# Function to perform likelihood ratio test on many varibles and predictors
ma_lrt <- function (dataList, test.var) {
  res <- list()
  
  for (i in 1:length(dataList)) {
    dat <- dataList[[i]]
    
    # quadratic term is not in the original data and has to be added here
    dat$N_conc_exp2_quad <- dat$N_conc_exp2 * dat$N_conc_exp2 
    
    print(names(dataList)[[i]]) # print out variable name
    #dat <- dat.photo
    out <- list()
    for (j in 1:length(test.var)) {
      varName <- test.var[j] 
      
      if(varName == "N_conc_exp2") {  # testing the degree of N limitation
        form.full <- as.formula(paste("~ ", varName, sep = ""))
        form.red <- as.formula("~ 1")
      } else { # All other variables are tested while controlling for degree of N limitation
        # Only for leaf N on a mass basis we detect a significant qudratic term. This term should therefore be included
        # when testing the other predictors (but when the quadratic term is tested, hence the "& !()" part)
        if(names(dataList)[[i]] == "leafN_mass" & !(varName == "N_conc_exp2_quad")) {
          form.full <- as.formula(paste("~ N_conc_exp2 + N_conc_exp2_quad +", varName, sep = ""))
          form.red <- as.formula("~ N_conc_exp2 + N_conc_exp2_quad") 
        } else {
          form.full <- as.formula(paste("~ N_conc_exp2 + ", varName, sep = ""))
          form.red <- as.formula("~ N_conc_exp2 ")
        }
      }
      
      dat.sub <- dat[!(is.na(dat[[varName]])),] # remove NAs of the specific variable tested
      
      # Check if the predictor, if a factor, has more than 1 level. If not, go to next variable.
      if(is.factor(dat.sub[[varName]]) & length(levels(dat.sub[[varName]]))<2) {next}
      # If predictor is a factor, then extract level names and sample size for each group
      if(is.factor(dat.sub[[varName]])) {
        level.names <- paste(names(summary(dat.sub[[varName]])), collapse = " : ") # levels of the factor
        level.samplesize <-paste(summary(dat.sub[[varName]]), collapse = ", ") # sample size per level
      } else {
        level.names <- NA
        level.samplesize <- NA
      }
      
      V <- v_func(dat = dat.sub) # create var-covar matrix
      mod.null.reml <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, # to get variance components w/o predictors 
                              random = list(~ 1 | studyNr,~ 1 | obsNr),
                              data = dat.sub, method="REML") # full model
      mod.full.reml <- rma.mv(yi = rr, mods = form.full, V = V, intercept=TRUE, # to get correct variance components 
                         random = list(~ 1 | studyNr,~ 1 | obsNr),
                         data = dat.sub, method="REML") # full model
      mod.full <- rma.mv(yi = rr, mods = form.full, V = V, intercept=TRUE, # ML for testing fixed effects
                         random = list(~ 1 | studyNr,~ 1 | obsNr),
                         data = dat.sub, method="ML") # full model
      mod.red <- rma.mv(yi = rr, mods = form.red, V = V, intercept=TRUE,  # ML for testing fixed effects
                        random = list(~ 1 | studyNr,~ 1 | obsNr),
                        data = dat.sub, method="ML") # reduced mode
      out[[j]] <- data.frame(LRT = round(anova(mod.full, mod.red)$LRT,3), # LRT value
                             df = anova(mod.full, mod.red)$p.f-anova(mod.full, mod.red)$p.r, # degrees of freedom
                             Pvalue = round(anova(mod.full, mod.red)$pval,3), # p-value
                             sigma.study = round(mod.full.reml$sigma2[1],4), # variance components, between study
                             sigma.exp = round(mod.full.reml$sigma2[2],4), # variance components, beteween experiment
                             samples.study = mod.full$s.nlevels[1], # sample size, number of studies
                             samples.exp = mod.full$s.nlevels[2], # sample size, total number of experiments
                             factors = level.names, # levels of the factor
                             group.size = level.samplesize, # sample size per level
                             expl.var = ( (mod.null.reml$sigma2[1] +mod.null.reml$sigma2[2]) -
                                                    (mod.full.reml$sigma2[1] +mod.full.reml$sigma2[2]))/
                                          (mod.null.reml$sigma2[1] +mod.null.reml$sigma2[2])
      ) 
      names(out)[[j]] <-varName # store variable name
    }
    
    res[[i]] <- do.call(rbind, out) # make summary data frame for the variable and store it in a list
    names(res)[[i]] <-names(dataList)[[i]] # name list index after response variable
  }
  return(res)
}  

# Calculate effect of N limitation for models with only
# proportion N limitation as covariate
Nlim_eff <- function (mod) { 
  b.mat <- cbind(mod$b,mod$ci.lb, mod$ci.ub)
  maxNlim <- (((exp(b.mat))-1)*100)[1,]
  slope <- c(((exp(b.mat)-1)*100)[1,1]-
               (exp(sum(b.mat[1:2,1]))-1)*100,
             ((exp(b.mat)-1)*100)[1,1]-
               (exp(sum(b.mat[1,1], b.mat[2,2]))-1)*100,
             ((exp(b.mat)-1)*100)[1,1]-
               (exp(sum(b.mat[1,1], b.mat[2,3]))-1)*100)
  out <- data.frame(rbind(maxNlim,slope))  
  colnames(out) <- c("Effect", "lower", "upper")
  return(out)
}


# Load data
# save as list of data frames
vars <- c("photosyn", "leaf area", "leafN_area", "leafN_mass","chl cont",
          "RubCont", "SLA", "leaf starch", "leaf sug")
          #,"protein leaf")
rep.var <- list()

  for (i in 1:length(vars)) {
               nn <- vars[i]
               resp <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                               sheet = nn)
               # remove all NA text
               resp[resp=="NA"] <- NA
               
               # change to factors
               # Weird with NAs in "leg" column
               cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
                        "croptype", "species", "CO2", "monodicot", "C3C4")    
               resp[,cols] = lapply(resp[ ,cols],  function(x) as.factor(x))
               
               # Independent studies are defined as unique controls
               resp$ind.study <- factor(paste(resp$author, resp$Xc, sep = "_"))
               
               # Run function to get effect sizes and (co)variances
               resp <- calc_resp(obj=resp)
               
               # order data according to ind studies
               resp <- resp[order(resp$ind.study),]
               
               rep.var[[i]] <- resp
  }
names(rep.var) <- c("photosynthesis", "leaf_area", "leafN_area", 
                    "leafN_mass", "leafChl", "leafRubisco", 
                    "sla", "leafStarch", "leafSugar")

# add response names as a new column
rep.names <- as.list(setNames(names(rep.var), names(rep.var)))
rep.var <- Map(cbind, rep.var, response = rep.names)

# Assign plotting colors 
# need to load photodata first!!
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col.species = levels(factor(rep.var[["photosynthesis"]]$species))
cols = gg_color_hue(length(col.species)) # ten species for photosynthesis
col.df <- data.frame(cols, col.species, stringsAsFactors = FALSE)

# END FUNCTIONS AND LOADING DATA

# Number of species
stud <- lapply(rep.var, FUN = function (x) x$species )
levels(factor(unlist(stud)))

# Number of studies
stud <- lapply(rep.var, FUN = function (x) x$studyNr )
factor(unlist(stud))

# Test differences between area and mass for Chl, Rub, carbs ####

# Chl
dat.chl <- rep.var[["leafChl"]]
summary(factor(dat.chl$unit2)) # mostly area based
summary(factor(dat.chl$unit2))/nrow(dat.chl) # 13% mass based

# Test for differences
V <- v_func(dat = dat.chl)
mod.full <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + unit2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="ML")
mod.red <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
                   random = list(~ 1 | studyNr,~ 1 | obsNr),
                   data = dat.chl, method="ML")
anova(mod.full, mod.red) # weak evidence for differences so we merge data
summary(mod.full)

# Rubisco
dat.rub <- rep.var[["leafRubisco"]]
summary(factor(dat.rub$unit2)) # only ONE mass based so no need to test

# Starch
dat.star <- rep.var[["leafStarch"]]
summary(factor(dat.star$unit2)) # mostly area based
summary(factor(dat.star$unit2))/nrow(dat.star) # 12% mass based

# Test for differences
V <- v_func(dat = dat.star)
mod.full <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + unit2, V = V, intercept=TRUE, 
                   random = list(~ 1 | studyNr,~ 1 | obsNr),
                   data = dat.star, method="ML")
mod.red <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
                  random = list(~ 1 | studyNr,~ 1 | obsNr),
                  data = dat.star, method="ML")
anova(mod.full, mod.red) # weak evidence for differences so we merge data
summary(mod.full)

# Sugar
dat.sug <- rep.var[["leafSugar"]]
summary(factor(dat.sug$unit2)) # mostly area based
summary(factor(dat.sug$unit2))/nrow(dat.sug) # 26% mass based

# Test for differences
V <- v_func(dat = dat.sug)
mod.full <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + unit2, V = V, intercept=TRUE, 
                   random = list(~ 1 | studyNr,~ 1 | obsNr),
                   data = dat.sug, method="ML")
mod.red <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
                  random = list(~ 1 | studyNr,~ 1 | obsNr),
                  data = dat.sug, method="ML")
anova(mod.full, mod.red) # weak evidence for differences so we merge data
summary(mod.full)

# Leaf area
dat.la <- rep.var[["leaf_area"]]
summary(factor(dat.la$organ)) # mostly whole plant based
summary(factor(dat.la$organ))/nrow(dat.la) # 14% leaf based

# Test for differences
V <- v_func(dat = dat.la)
mod.full <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + organ, V = V, intercept=TRUE, 
                   random = list(~ 1 | studyNr,~ 1 | obsNr),
                   data = dat.la, method="ML")
mod.red <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
                  random = list(~ 1 | studyNr,~ 1 | obsNr),
                  data = dat.la, method="ML")
anova(mod.full, mod.red) # weak evidence for differences so we merge data
summary(mod.full)


# Test experiment setup factors and other predictors ####
test.var <- c("N_conc_exp2", "N_conc_exp2_quad", "CO2", # Treatment effects 
  "lengthlim","frequencyN", "Nsource", "pHcontr",  "facility", "medium", "potsize", # Possible experimental setup effects
  "leg", "croptype", "monodicot", "C3C4", "species") # Taxonomy groups 

# run lrt function
new <- ma_lrt(dataList = rep.var, test.var = test.var)  # runt test

# build large table with all response variables
large_tab <- do.call(rbind, new)
large_tab <- cbind(do.call(rbind, strsplit(rownames(large_tab), "[.]")), large_tab)
colnames(large_tab)[1:2] <- c("response", "predictor")
rownames(large_tab) <- NULL
write.csv(file = "large_table.csv", large_tab) #write table to file


# View data: plot responses in relation to predictors ####
# PLots to get an overview of the response ratio (RR)
# in relation to all predictors

# Define predictors to plot (+ the response column, "rr")
test.var <- c("rr", # response
              "N_conc_exp2", "CO2", # Treatment effects 
              "lengthlim","frequencyN", "Nsource", "pHcontr",  "facility", "medium", "potsize", # Possible experimental setup effects
              "leg", "croptype", "monodicot", "C3C4", "species") # Taxonomy groups 

# Plot all responses as separate files
lapply(rep.var, function (x) {
  fig.title <- x$response[1] # name of the figure
  file.name <- paste(fig.title,"_RR_all_var.png", sep="") # file name
  png(file.name, width = 18, height = 30, units ="cm", res=600) # initiate plot
  par(mfrow = c(5,3))
  X.sub <- x[, test.var] # store subset of data
    
    # run through all predictors  
    for (i in 1:(ncol(X.sub)-1) ) {
          tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub)[1+i]) # give plot title
          plot(X.sub[,"rr"] ~X.sub[,1+i], ylab = "log RR", xlab="", 
                main = tit, las=2)
          
          # extract sample size for each group, including NAs
          n <- as.numeric(table(X.sub[, 1+i]))
          nas <- sum(is.na(X.sub[, 1+i]))
          text(x = seq_along(n),y = -2, label = n)
          text(x = median(seq_along(n)),y = -2.5, label = paste("NAs=",nas), font=2)
          }
  dev.off()
  })


###########
# Compare results with full Bayes: MCMCglmm package ####
library(MCMCglmm)

# test with photosynthesis
dat.photo <- rep.var[[1]]

# setup var-covar matrix
Rsvd <- svd(V)
Rsvd <- Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.photo$exp <- as.factor(1:NROW(dat.photo))
Z <- model.matrix(~exp-1, dat.photo)%*%Rsvd
dat.photo$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

mod.leg <- MCMCglmm(rr ~ N_conc_exp2 + leg,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(mod.leg)
plot(mod.leg$VCV)

mod.Nlim <- MCMCglmm(rr ~ N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(mod.Nlim)
plot(mod.Nlim$VCV)

# fit same model with REML in metafor
V <- v_func(dat = dat.photo) # create var-covar matrix
mod.reml <- rma.mv(yi = rr, mods = N_conc_exp2, V = V, intercept=TRUE, 
                   random = list(~ 1 | studyNr,~ 1 | obsNr),
                   data = dat.photo, method="REML")

res.mc <- dat.photo$rr-predict(mod.Nlim) # raw residuals in MCMCglmm

# predicted minus predicted not including the within-study part
# so we get the residuals for "study dependence"
res.st.mc2 <- predict(mod.Nlim)-predict(mod.Nlim, marginal = ~studyNr ) 

# add study residuals to raw data to estiamte adjusted study values
res.st.mc2 <- dat.photo$rr + res.st.mc2 

# plot raw rr, adjusted rr (for study effect), and predicted slope (green)
plot(residuals.rma(mod.reml), res.mc)

plot(rr ~ N_conc_exp2, data=dat.photo, ylim=c(-1.5, 0.5))
points(res.st.mc2~ dat.photo$N_conc_exp2, col="red")
points(fitted.rma(mod.reml) ~dat.photo$N_conc_exp2, col="green")

# check that fitted and predicted give same value
plot(fitted.rma(mod.reml), predict(mod.Nlim)) 


# fix this below...for the leg model
m3 <- MCMCglmm(rr ~ newExpLength + facility + leg,  
               random=~idv(Z) + studyNr, 
               data=dat.photo[!(is.na(dat.photo$leg)),], 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)

summary(m3)

library(reshape2)
colnames(m3$Sol)
ref <- data.frame(m3$Sol[,185:233])
#colnames(ref) <- letters[1:4]
rdf <- melt(ref)
ggplot(data=rdf,aes(x=value, color=variable))+geom_density()
########### END MCMCglmm testing











# Photosynthesis ####
# save data
dat.photo <- rep.var[["photosynthesis"]]

# load data 
#library(readxl)
#dat.photo <- read_xls("Seufert_etal_Nlim_meta_data.xls",
#                       sheet = "photosyn")

# remove all NA text
#dat.photo[dat.photo=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
#cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
#         "croptype", "species", "CO2", "monodicot", "C3C4")    
#dat.photo[,cols] = lapply(dat.photo[ ,cols],  function(x) as.factor(x))


# Analyses show only a difference between shorter and long ("entire") experiments. 
# So we can use two categories instead of three to account for experiment length.
# There are, however, many NAs (studies that dont indicate the length of the experiment) 
# and we assign them to be, on average, in the middle of short and long experiments 
# (i.e. a value of 0.5)
# Short studies are coded as 0, long (entire) as 1, and the NAs as 0.5.
dat.photo$newExpLength <- ifelse(dat.photo$lengthlim == "<1/2" |
                                   dat.photo$lengthlim == ">1/2", 0, 1)
dat.photo$newExpLength[is.na(dat.photo$newExpLength)] <- 0.5

# fix data ####

# Independent studies are defined as unique controls
#dat.photo$ind.study <- factor(paste(dat.photo$author, dat.photo$Xc, sep = "_"))

# Run function to get effect sizes and (co)variances
#dat.photo <- calc_resp(obj=dat.photo)

# order data according to ind studies
#dat.photo <- dat.photo[order(dat.photo$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.photo, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")

# Explore data
plot(rr ~ Ne, dat.photo)
plot(rr ~ var.rr, dat.photo)
summary(lm(rr ~ var.rr, dat.photo))
dat.photo[dat.photo$var.rr>5,] # check observation 61

ggplot(data=dat.photo, aes(x=N_conc_exp2, y=rr)) +
  geom_point() +
  facet_wrap(~Nsource)

# covariation among explanatory variables
X.sub <- dat.photo[, c("rr", # RESPONSE
                       "N_conc_exp2", "lengthlim","frequencyN", "Nsource", # Possible Treatment effects 
                       "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                       "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("photo_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
#rot <- ifelse(is.factor(X.sub[,1+i]), 2,1) 
for (i in 1:(ncol(X.sub)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub)[1+i])
  plot(X.sub[,"rr"] ~X.sub[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2) #xaxt="n"
  #
  #labs <- levels(X.sub[,1+i])
  #text(cex=1, x=x-.25, y=-1.25, labs, xpd=TRUE, srt=90, pos=2) }
  #axis(1, at= seq(1, 10, by=1), labels = FALSE)
  #text(seq(1, 10, by=1), par("usr")[3] - 0.2, labels = lablist, srt = 45, pos = 1, xpd = TRUE)
  n <- as.numeric(table(X.sub[, 1+i]))
  nas <- sum(is.na(X.sub[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -2, label = n)
  text(x = median(seq_along(n)),y = -2.5, label = paste("NAs=",nas), font=2)
  
}
dev.off()
#X.sub$newExpLength <- dat.photo[, "newExpLength"] # add the "new" exp length variable

# Potsize, pH, frequency of application, general stress and medium, did not show any indications
# of being important. Lets look at other experiment set up variables.
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ facility + newExpLength + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ facility + newExpLength, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod3 <- rma.mv(yi = rr, mods = ~ newExpLength + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod4 <- rma.mv(yi = rr, mods = ~ facility + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")

anova(mod1, mod2) # N lim
anova(mod1, mod3) # facility
anova(mod1, mod4) # length
# when including study as random factor, there is no strong
# evidence of an effect of length and facility.
# Few studies are associated with a particular facility or study length

# No data on N source for two studies.
# Lets investigate N source effect.
dat.Nsour <- dat.photo[!(is.na(dat.photo$Nsource)),]
#dat.Nsour <- dat.photo[!(is.na(dat.photo$Nsource)) & !(dat.photo$leg == "yes, nod"),] # test w/o leg with nod. Now no N source effect
V2 <- v_func(dat = dat.Nsour)
mod5 <- rma.mv(yi = rr, mods = ~ Nsource + N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nsour, method="ML")
mod6 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nsour, method="ML")
anova(mod5, mod6) # N source effect
# Run REML for estimation
mod.Nsource <- rma.mv(yi = rr, mods = ~ Nsource + N_conc_exp2, V = V2, intercept=TRUE, 
                      random = list(~ 1 | studyNr,~ 1 | obsNr),
                      data = dat.Nsour, method="REML")
summary(mod.Nsource)
# strongest effect when nitrate and ammonium is combined. Significantly weaker effect if only
# ammonium is applied. However, we can see that the effect of N lim when ignoring N source, is close 
# to the effect of nitrate and nitrate + ammonium. Hence, ignoring N source is not changing the overall 
# estimate. Looking at the numbers, if nitrate+ammonium is used as the "correct" estiamte at maximum
# N limitation, then the difference is 5%, -34% vs -39%.


# REML models for estimation. 
# Compare estimation of the only N lim model, with 'many covariates' model
# and 'N source covariate' model.
mod0 <- rma.mv(yi = rr, mods = ~ Nsource + facility + newExpLength + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
mod1 <- rma.mv(yi = rr, mods = ~ Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]])
summary(mod0)
summary(mod1)
summary(mod2)
# Overall, maybe some (weak) evidence that longer studies weakens the effect
# and that chamber and outside gives a weaker response.But these results are linked 
# to specific studies.
# If the model weigh studies against sample size, the effect at maximum N limitation do not change,
# but the N limitation regression is less strong and P=0.08.


# % change under complete N limitation and slope of N lim effect
(exp(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub))-1)*100

ggplot(dat.Nsour, aes(x=N_conc_exp2, y=rr, color=Nsource)) +
  geom_point() +
  facet_wrap(~leg)


# Now test the other predictors of interest

# crop type
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
anova(mod1, mod2)
# no effect of crop type

# C3-C4
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
anova(mod1, mod2)
# no effect of C4

# CO2
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
anova(mod1, mod2)
# P=0.05. Barely significant effect but some evidence.
# However, if controlling for N source, then P=0.02 but estimate almost identical. 
# When uncluding N source we lose 2 studies though.
# lets check the model
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
summary(mod1)
c(exp(coef(mod1)[1]), exp(coef(mod1)[1]+coef(mod1)[3])) # effects at ambient and elevated CO2
# slightly stronger effect of N limitation for elevated CO2 but
# effect small (-33% for ambient CO2, -40% for elevated CO2)


# mono-dicot
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
anova(mod1, mod2)
# P=0.07. Not significant effect but some evidence. However, if controlling for N source, then 
# P = 0.25.
# lets check the model
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
summary(mod1)
# slightly stronger effect of N limitation for monocots but
# effect small (-30% for dicot, -36% for monocot). However, mono-dicot is confounded
# with N source and almost no effect if N source is accounted for. 

# legume
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
anova(mod1, mod2)
# Clear effect of legumes

# pH
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + pHcontr, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
anova(mod1, mod2)
# Nothing

# stress
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + stress, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="ML")
anova(mod1, mod2)
# Nothing

# frequencyN
dat.Nfreq <- dat.photo[!(is.na(dat.photo$frequencyN)),]
V2 <- v_func(dat = dat.Nfreq)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + frequencyN, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nfreq, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nfreq, method="ML")
anova(mod1, mod2)
# Nothing

# pot size
dat.pot <- dat.photo[!(is.na(dat.photo$potsize)),]
V2 <- v_func(dat = dat.pot)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + potsize, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.pot, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.pot, method="ML")
anova(mod1, mod2)
# Nothing

# medium
dat.medium <- dat.photo[!(is.na(dat.photo$medium)),]
V2 <- v_func(dat = dat.medium)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + medium, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.medium, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.medium, method="ML")
anova(mod1, mod2)
# Nothing


# Estimate coefs of N limitation, including legume status
# N source not included but could be. Dont change the results though
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")
summary(mod1)
c(exp(coef(mod1)[1]), exp(coef(mod1)[1]+coef(mod1)[3]),
  exp(coef(mod1)[1]+coef(mod1)[4])) # effects at: not legum, legum with nods, legume w/o nods 
# Legumes with nods removes the N limitation effect completely 

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.photo$rr)-1)*100,
           "PropNlim" =  1-dat.photo$N_conc_exp2)
dat.photo$N_conc_exp2.rev <- 1-dat.photo$N_conc_exp2
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.photo, method="REML")

pred.frame <-  predict(mod1, newmods = c(seq(0.25,1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.photo <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x=X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  scale_y_continuous( name = "Change in photosynthetic rate (%)", 
                      breaks = seq(-100,100,10)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimPhoto.png", fig.Nlim.photo, dpi = 300)




# Possible confounding factors
# time
mf.conf <- rma.mv(yi = rr, mods = ~ lengthlim , V = V, intercept=FALSE, random = ~ 1 | obsNr,
                     data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$lengthlim)
# maybe important but too many NAs

# N source
mf.conf <- rma.mv(yi = rr, mods = ~ Nsource, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$Nsource)
# seems important to control for

# N freq
mf.conf <- rma.mv(yi = rr, mods = ~ frequencyN, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$frequencyN)
# NOT important

# Facility
mf.conf <- rma.mv(yi = rr, mods = ~ facility, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$facility)
# Maybe but not vert consistent....weird result. And many NAs

# pH
mf.conf <- rma.mv(yi = rr, mods = ~ pHcontr, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$pHcontr)
# Doesnt seem to be important

# stress
mf.conf <- rma.mv(yi = rr, mods = ~ stress, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$stress)
# NOTHING

# potsize
mf.conf <- rma.mv(yi = rr, mods = ~ potsize, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$potsize)
# Smaller = smaller effect maybe...not that clear though?

# medium
mf.conf <- rma.mv(yi = rr, mods = ~ medium, V = V, intercept=FALSE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mf.conf)
table(dat.photo$medium)
# lower effect with soil...make sense I guess (some N in soil). Many NAs

# Explore and check covaraition Photosynthesis####
X.sub <- dat.photo[, c("rr", # RESPONSE
                     "N_conc_exp2", "lengthlim","frequencyN", "Nsource", # Possible Treatment effects 
                     "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                     "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("photo_RR_all.png", width = 1000, height = 1100)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub)-1) ) {
  plot(X.sub[,"rr"] ~X.sub[,1+i], ylab = "log RR", xlab = colnames(X.sub)[1+i] )
  n <- as.numeric(table(X.sub[, 1+i]))
  nas <- sum(is.na(X.sub[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -2, label = n)
  text(x = median(seq_along(n)),y = -2.5, label = paste("NAs=",nas), font=2)
}
dev.off()
X.sub$newExpLength <- dat.photo[, "newExpLength"] # add the "new" exp length variable

# Test experimental non-N factors that may effect the outcome
# we dont care about them but should control for it if needed

# Function to get Cramers V
get.V<-function(y){
  col.y<-ncol(y)
  V<-matrix(ncol=col.y,nrow=col.y)
  for(i in 1:col.y){
    for(j in 1:col.y){
      V[i,j]<-assocstats(table(y[,i],y[,j]))$cramer
    }
  }
  return(V)
}
library(vcd)
# Cramsers V for these two variables
get.V(X.sub[,c("facility", "lengthlim")])
table(X.sub[,c("facility", "lengthlim")])

# Cramsers V for N variables
get.V(X.sub[,c("frequencyN", "Nsource")])
table(X.sub[,c("frequencyN", "Nsource")])

# two very correlated
get.V(X.sub[,c("facility", "Nsource")])
table(X.sub[,c("facility", "Nsource")])
#### missing values need to be FOUND

get.V(X.sub[,c("species", "facility")])
table(X.sub[,c("species", "facility")])
# O. sativa in all facilities

get.V(X.sub[,c("species", "potsize")])
table(X.sub[,c("species", "potsize")])
# G. max dominates small

get.V(X.sub[,c("facility", "potsize")])
table(X.sub[,c("species", "facility", "potsize")])
# pot size correlates with facility, mostly
# small pots outside

get.V(X.sub[,c("species", "medium")])
table(X.sub[,c("species", "medium", "facility")])
# 

get.V(X.sub[,c("facility", "medium")])
table(X.sub[,c("facility", "medium")])
# little correlation

get.V(X.sub[,c("C3C4", "Nsource")])
table(X.sub[,c("C3C4", "Nsource")])
# C3C4 doesnt correlate with anything basiclly

get.V(X.sub[,c("species", "Nsource")])
table(X.sub[,c("species", "facility")])
# G. max in all N sources

#### Model testing ####

# O. sativa in all facilities.
table(X.sub[,c("species", "facility")])
# test facility effect
test.fac <- dat.photo[dat.photo$species == "O. sativa" | 
                        dat.photo$species == "Z. mays"| 
                        dat.photo$species == "G. max",]
V <- v_func(dat = test.fac)
mod1 <- rma.mv(yi = rr, mods = ~ facility, V = V, intercept=TRUE, random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.fac, method="REML")
summary(mod1)
# Outside dampen the N lim effect. Outside is perfectly correlated with only amonium addition. 
# I.e. we cant separate outside and amoniumm effect. 

# G. max has all N sources in growth chamber facility.
# test facility effect
table(X.sub[,c("species", "Nsource", "facility")])
test.Nsource <- dat.photo[dat.photo$species == "G. max" | dat.photo$species == "T. aestivum" & 
  dat.photo$facility == "growth chamber",]
V <- v_func(dat = test.Nsource)
mod2 <- rma.mv(yi = rr, mods = ~ Nsource, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.Nsource, method="REML")
summary(mod2)
# Small sample size but not much evidence for a strong ammonium effect (when added alone).
# CONLUCION: outdoor likely more the mechanisms than ammonium treatment.
# N source can be "ignored"

# test potsize
# when testing on all data, small pots dampens the negative N lim effect.
# here I test only within growth chambers because all sizes are represented
table(X.sub[,c("species", "potsize", "facility")])
# lets test within growth chambers
test.potsize <- dat.photo[dat.photo$facility == "growth chamber",]
V <- v_func(dat = test.potsize)
mod3 <- rma.mv(yi = rr, mods = ~ potsize, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.potsize, method="REML")
summary(mod3)
# Now small has a stronger negative effect, medium size as well, comapred to big. So the opposite.
# it is questionable if we can say anything general about pot size.

# test medium
# when testing on all data, small pots dampens the negative N lim effect.
# here I test only within growth chambers because all sizes are represented
table(X.sub[,c("species", "potsize", "facility")])
# lets test within growth chambers
test.medium <- dat.photo[dat.photo$facility == "growth chamber" & 
                            dat.photo$species == "G. max" |
                            dat.photo$species == "T. aestivum"|
                            dat.photo$species == "Z. mays",]
V <- v_func(dat = test.medium)
mod3 <- rma.mv(yi = rr, mods = ~ medium, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.medium, method="REML")
summary(mod3)
# No evidence of a medium effect


# test N frequency
table(X.sub[,c("species", "frequencyN", "Nsource")])
# lets test within growth chambers
V <- v_func(dat = dat.photo)
mod3 <- rma.mv(yi = rr, mods = ~ frequencyN, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
summary(mod3)
# No evidence of a N frequency effect

# test pH
table(X.sub[,c("species", "pHcontr", "facility")])
# lets test without outside exp
test.ph <- dat.photo[!(dat.photo$facility == "pots outside"),]
V <- v_func(dat = test.ph)
mod3 <- rma.mv(yi = rr, mods = ~ pHcontr, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.ph, method="REML")
summary(mod3)
# now all data
V <- v_func(dat = dat.photo)
mod3 <- rma.mv(yi = rr, mods = ~ pHcontr, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
summary(mod3)
# Small pH effect when all data included but this is linked to the many outdoor exp
# that did not control pH. When restricted to greenhouse/chambers no effect.


# test study length
table(X.sub[,c("species", "lengthlim", "facility")])
# lets test within growth chambers
V <- v_func(dat = dat.photo)
mod3 <- rma.mv(yi = rr, mods = ~ lengthlim, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
summary(mod3)
# longer studies = samaller effects
# but longer studies often outside
# lets test without outside
test.length <- dat.photo[!(dat.photo$facility == "pots outside"),]
V <- v_func(dat = test.length)
mod3 <- rma.mv(yi = rr, mods = ~ lengthlim, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.length, method="REML")
summary(mod3)
# still an effect
# do we see the same trend in a specific species that have all exp lengths, within a facility?
test.length <- dat.photo[dat.photo$species == "O. sativa" & dat.photo$facility == "greenhouse",]
V <- v_func(dat = test.length)
mod3 <- rma.mv(yi = rr, mods = ~ lengthlim, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.length, method="REML")
summary(mod3)
# same trend, although not significant but very low sample size
# Problem with many NAs
get.V(X.sub[,c("facility", "lengthlim")])
# facility and length not very correlated

# Only a difference between shorter and "entire" experiments. So we can make two categories.
# Many NAs and we put them in between
dat.photo$newExpLength <- ifelse(dat.photo$lengthlim == "<1/2" |
                                   dat.photo$lengthlim == ">1/2", 0, 1)
dat.photo$newExpLength[is.na(dat.photo$newExpLength)] <- 0.5

V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ newExpLength, V = V, intercept=TRUE, random = ~ 1 | obsNr,
                  data = dat.photo, method="REML")
summary(mod1)

# Full models
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ facility + newExpLength + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
plot(fitted.rma(mod1), rstandard.rma.mv(mod1)[["resid"]])
abline(h=0)
qqnorm(rstandard.rma.mv(mod1)[["resid"]])
abline(a=0,b=1)

### calculate influence diagnostics
inf <- influence(mod1)
plot(cooks.distance.rma.mv(mod1)) 
### plot the influence diagnostics
plot(inf, layout=c(8,1))summary(mod1)


# Check correlation between facility and newExpLength
X.sub$newExpLength <-  dat.photo$newExpLength
get.V(X.sub[,c("facility", "newExpLength", "medium")])
table(X.sub[,c("facility", "medium")])

# Crop type
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
summary(mod1)
# No differences. Slightly higher for oil seed but marginal significant

# What if we control for exp set ups
mod1 <- rma.mv(yi = rr, mods = ~ newExpLength + facility + croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
summary(mod1)
# No difference

# CO2
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ CO2, V = V, intercept=TRUE, random = ~ 1 | obsNr,
               data = dat.photo, method="REML")
summary(mod1)
# Not significant. Slightly stronger N lim effect with higher CO2

# What if we control for exp set ups
mod1 <- rma.mv(yi = rr, mods = ~ newExpLength + facility + CO2, V = V, intercept=TRUE, random = ~ 1 | obsNr,
               data = dat.photo, method="REML")
summary(mod1)
# Marginal significant. Slightly stronger N lim effect with higher CO2

# C3C4
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ C3C4, V = V, intercept=TRUE, random = ~ 1 | obsNr,
               data = dat.photo, method="REML")
summary(mod1)
# No effect.

# What if we control for exp set ups
mod1 <- rma.mv(yi = rr, mods = ~ newExpLength + facility + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr),
               data = dat.photo, method="REML")
summary(mod1)
# No effect.

# mono/dicot
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ monodicot, V = V, intercept=TRUE, random = ~ 1 | obsNr,
               data = dat.photo, method="REML")
summary(mod1)
# weak effect, marginally significant
# What if we control for exp set ups
mod1 <- rma.mv(yi = rr, mods = ~ newExpLength + facility + monodicot, V = V, 
               intercept=TRUE, random = ~ 1 | obsNr,
               data = dat.photo, method="REML")
summary(mod1)
# No effect of mono-dicot

# leg
V <- v_func(dat = dat.photo)
mod1 <- rma.mv(yi = rr, mods = ~ leg, V = V, intercept=TRUE, random = ~ 1 | obsNr,
               data = dat.photo, method="REML")
summary(mod1)
# Effect of nod-forming
# What if we control for exp set ups
mod1 <- rma.mv(yi = rr, mods = ~ facility + newExpLength + leg, V = V, 
               intercept=TRUE, random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.photo, method="REML")
summary(mod1)
plot(fitted.rma(mod1), rstandard.rma.mv(mod1)[["resid"]])
# Effect of nod-forming


###########
# with MCMCglmm
# match metafor results
library(MCMCglmm)
Rsvd<-svd(V)
Rsvd<-Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.photo$exp <- as.factor(1:NROW(dat.photo))
Z <- model.matrix(~exp-1, dat.photo)%*%Rsvd
dat.photo$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

m2 <- MCMCglmm(rr ~ newExpLength + facility + N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(m2)
plot(m2$VCV)

m2 <- MCMCglmm(rr ~ N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, prior=prior,nitt=100000,burnin=15000, pr=TRUE)
#res.mc <- dat.photo$rr-predict(m2) # raw residuals in MCMCglmm
res.st.mc2 <- predict(m2)-predict(m2, marginal = ~studyNr ) 
# predicted minus predicted not including the within-study part
# so we get the residuals for "study dependence"

res.st.mc2 <- dat.photo$rr + res.st.mc2 
# add study residuals to raw data to estiamte adjusted study values

# plot raw rr, adjusted rr (for study effect), and predicted slope (green)
plot(residuals.rma(mod1), res.mc)
plot(rr ~ N_conc_exp2, data=dat.photo, ylim=c(-1.5, 0.5))
points(res.st.mc2~ dat.photo$N_conc_exp2, col="red")
points(fitted.rma(mod1) ~dat.photo$N_conc_exp2, col="green")

plot(fitted.rma(mod1), predict(m2)) # check that fitted and predicted give same value


m3 <- MCMCglmm(rr ~ newExpLength + facility + leg,  
               random=~idv(Z) + studyNr, 
               data=dat.photo[!(is.na(dat.photo$leg)),], 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)
m4 <- MCMCglmm(rr ~ species,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)

m4 <- MCMCglmm(rr ~ CO2,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)

summary(m3)
plot(m4$VCV)
plot(m4$Sol)

library(reshape2)
colnames(m3$Sol)
ref <- data.frame(m3$Sol[,185:233])
#colnames(ref) <- letters[1:4]
rdf <- melt(ref)
ggplot(data=rdf,aes(x=value, color=variable))+geom_density()
###########


# Leaf area ####
#library(xlsx)
library(readxl)
#vars <- c("photosyn","chl cont","starAct", "starCont", "leaf area", "leafN_area", "leafN_mass","protein leaf", "leaf sug", "leaf starch", "SLA")

#dat.leafA <- read.xlsx("Seufert_etal_Nlim_meta_data.xls",
#                       sheetIndex = "leaf area", stringsAsFactors=FALSE)
dat.leafA <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                       sheet = "leaf area")

#str(dat.photo)
#merge(dat.photo, dat.leafA, by="author_species_exp")
#colnames(dat.photo)
#colnames(dat.leafA)

# remove all NA text
dat.leafA[dat.leafA=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.leafA[,cols] = lapply(dat.leafA[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.leafA$ind.study <- factor(paste(dat.leafA$author, dat.leafA$Xc, sep = "_"))

# Calculate effect sizes
dat.leafA <- calc_resp(obj=dat.leafA)

# order data according to independent studies
dat.leafA <- dat.leafA[order(dat.leafA$ind.study),]

# Get known variance-covariance matrix
V <- v_func(dat = dat.leafA, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")

#plot rr against N and variance
plot(rr ~ Ne, dat.leafA)
plot(rr ~ var.rr, dat.leafA)
summary(lm(rr ~ var.rr, dat.leafA))

# Explore and check covariation LEAF AREA ####
X.sub.leaf <- dat.leafA[, c("rr", # RESPONSE
                       "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                       "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                       "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("leafA_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.leaf)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.leaf)[1+i])
  plot(X.sub.leaf[,"rr"] ~X.sub.leaf[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2)
  n <- as.numeric(table(X.sub.leaf[, 1+i]))
  nas <- sum(is.na(X.sub.leaf[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -2, label = n)
  text(x = median(seq_along(n)),y = -2.5, label = paste("NAs=",nas), font=2)
}
dev.off()
#dat.leafA[is.na(dat.leafA$Nsource),] # NAs of N source

# Analyze effects with metafor

# Final N lim effect regression model
V <- v_func(dat = dat.leafA)
mod1 <- rma.mv(yi = rr, mods = ~  Nsource + leg + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ leg + Nsource, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod3 <- rma.mv(yi = rr, mods = ~ Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
V2 <- v_func(dat = dat.leafA[!(is.na(dat.leafA$Nsource)),]) # remove NAs
mod4 <- rma.mv(yi = rr, mods = ~ leg + N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA[!(is.na(dat.leafA$Nsource)),], method="ML")

anova(mod1, mod2) # N lim
anova(mod1, mod3) # legume
anova(mod1, mod4) # Nsource
# when including study as random factor, there is no strong

# No effects of facility, frequency, medium, pHcontr, experimental length
V <- v_func(dat = dat.leafA)
mod0 <- rma.mv(yi = rr, mods = ~  N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod1 <- rma.mv(yi = rr, mods = ~  facility + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~  pHcontr + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
#mod3 <- rma.mv(yi = rr, mods = ~  frequencyN + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafA, method="ML") # many NAs
mod4 <- rma.mv(yi = rr, mods = ~  medium + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
#mod5 <- rma.mv(yi = rr, mods = ~  lengthlim + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafA, method="REML") # some NAs but no effect

#mod6 <- rma.mv(yi = rr, mods = ~  potsize + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafA, method="ML") 
# some NAs. Some evidence and more pot size is further explored below

anova(mod0, mod1) # facility
anova(mod0, mod2) # pH
#anova(mod0, mod3) # N application frequency
anova(mod0, mod4) # medium

# pot size less clear. Some NAs but there may be an effect - weaker effect 
# in smaller pots. This seems to be driven by few studies though.
# Also, partly correlated with nitrate addition.
dat.Npot <- dat.leafA[!(is.na(dat.leafA$Nsource) | is.na(dat.leafA$potsize)),]
dat.Npot <- droplevels(dat.Npot)
V2 <- v_func(dat = dat.Npot)
mod5 <- rma.mv(yi = rr, mods = ~ potsize+Nsource + N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Npot, method="ML")
mod6 <- rma.mv(yi = rr, mods = ~  Nsource + N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Npot, method="ML")
anova(mod5, mod6) # potsize effect


# REML models for estimation. Compare estimation of the only N lim model, with 'many covariates' model
mod1 <- rma.mv(yi = rr, mods = ~ leg + Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="REML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="REML")
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]]) # looks good!
summary(mod1)
summary(mod2)

# Test for only N lim effect
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
anova(mod1, mod2)

# N source seems to be important but two studies are missing data.
# Lets investigate this effect.
table(dat.Nsour[,c("croptype", "Nsource")])

dat.Nsour <- dat.leafA[!(is.na(dat.leafA$Nsource)),]
#dat.Nsour <- dat.leafA[!(is.na(dat.leafA$Nsource)) & !(dat.leafA$leg == "yes, nod"),]
#dat.Nsour <- dat.leafA[!(is.na(dat.leafA$Nsource)) & !(dat.leafA$studyNr == "study10"),]
dat.Nsour <- droplevels(dat.Nsour)
V2 <- v_func(dat = dat.Nsour)
mod5 <- rma.mv(yi = rr, mods = ~ leg+Nsource + N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nsour, method="ML")
mod6 <- rma.mv(yi = rr, mods = ~ leg +N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nsour, method="ML")
anova(mod5, mod6) # N source effect
# Run REML for estimation
dat.Nsour$quad.Nconc <- dat.Nsour$N_conc_exp2*dat.Nsour$N_conc_exp2
mod.Nsource <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + Nsource, V = V2, intercept=TRUE, 
                      random = list(~ 1 | studyNr,~ 1 | obsNr),
                      data = dat.Nsour, method="REML")

# % change under complete N limitation and slope of N lim effect
(exp(cbind(mod.Nsource$b,mod.Nsource$ci.lb, mod.Nsource$ci.ub))-1)*100

ggplot(dat.Nsour, aes(x=N_conc_exp2, y=rr, color=Nsource)) +
  geom_point() +
  facet_wrap(~leg)
# Now test the other predictors of interest

# crop type
# controlling for N source
table(dat.Nsour[,c("croptype", "Nsource")])

dat.Nsour <- dat.leafA[!(is.na(dat.leafA$Nsource)),]
dat.Nsour <- droplevels(dat.Nsour)
V2 <- v_func(dat = dat.Nsour)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + croptype, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nsour, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V2, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.Nsour, method="ML")
anova(mod1, mod2)
# no effect of crop type

# C3-C4
mod1 <- rma.mv(yi = rr, mods = ~ Nsource+N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ Nsource+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
anova(mod1, mod2)
# no effect of C4

# CO2
mod1 <- rma.mv(yi = rr, mods = ~ Nsource+N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ Nsource+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
anova(mod1, mod2)
# Not significant

# mono-dicot
mod1 <- rma.mv(yi = rr, mods = ~ Nsource+N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ Nsource+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
anova(mod1, mod2)
# Not significant

# legume
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="ML")
anova(mod1, mod2)
# P=0.06 (0.08 if controlling for N source)
# so some indications, lets look at the model
# Estimate coefs
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="REML")
summary(mod1)
# Legumes with nods almost removes the N limitation effect completely 

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.leafA$rr)-1)*100,
                        "PropNlim" =  1-dat.leafA$N_conc_exp2)
dat.leafA$N_conc_exp2.rev <- 1-dat.leafA$N_conc_exp2
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafA, method="REML")
nrow(dat.leafA) # sample size

pred.frame <-  predict(mod1, newmods = c(seq(0.4,1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.leafA <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x=X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf area (%)", 
                      breaks = seq(-100,100,10)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimLeafA.png", fig.Nlim.leafA, dpi = 300)






##############################################33
# test N source
# strong overall effect of nitrate. No obvious correlation with other factors
# but mainly 2-3 studies that have this strong effect.
table(X.sub.leaf[,c("Nsource", "frequencyN")])
V <- v_func(dat = dat.leafA)
mod1 <- rma.mv(yi = rr, mods = ~ facility + Nsource, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafA, method="REML")
summary(mod1)
# lets test within a species that has all N sources
test.length <- dat.leafA[(dat.leafA$species == "G. max"),]
V <- v_func(dat = test.length)
mod3 <- rma.mv(yi = rr, mods = ~ Nsource, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = test.length, method="REML")
summary(mod3)
# still consistent effect


# with MCMCglmm
# match metafor results
library(MCMCglmm)
V <- v_func(dat = dat.leafA)
Rsvd<-svd(V)
Rsvd<-Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.leafA$exp <- as.factor(1:NROW(dat.leafA))
Z <- model.matrix(~exp-1, dat.leafA)%*%Rsvd
dat.leafA$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

m2 <- MCMCglmm(rr ~ Nsource + N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.leafA[!(is.na(dat.leafA$Nsource)),], 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(m2)
plot(m2$Sol)


# Leaf N area ####
#library(xlsx)
library(readxl)
#vars <- c("photosyn","chl cont","starAct", "starCont", "leaf area", "leafN_area", "leafN_mass","protein leaf", "leaf sug", "leaf starch", "SLA")

#dat.leafNA <- read.xlsx("Seufert_etal_Nlim_meta_data.xls",
#                       sheetIndex = "leafN_area", stringsAsFactors=FALSE)
dat.leafNA <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                        sheet = "leafN_area")

#str(dat.photo)
#merge(dat.photo, dat.leafA, by="author_species_exp")
#colnames(dat.photo)
#colnames(dat.leafA)

# remove all NA text
dat.leafNA[dat.leafNA=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.leafNA[,cols] = lapply(dat.leafNA[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.leafNA$ind.study <- factor(paste(dat.leafNA$author, dat.leafNA$Xc, sep = "_"))

# Calc effect size
dat.leafNA <- calc_resp(obj=dat.leafNA)

# order data according to ind studies
dat.leafNA <- dat.leafNA[order(dat.leafNA$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.leafNA, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.leafNA)
plot(rr ~ var.rr, dat.leafNA)
summary(lm(rr ~ var.rr, dat.leafNA))

# Explore and check covariation LEAF N AREA ####
X.sub.leafNA <- dat.leafNA[, c("rr", # RESPONSE
                            "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                            "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                            "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("leafNA_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.leafNA)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.leafNA)[1+i])
  plot(X.sub.leafNA[,"rr"] ~X.sub.leafNA[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-1.75, 0.25))
  n <- as.numeric(table(X.sub.leafNA[, 1+i]))
  nas <- sum(is.na(X.sub.leafNA[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -1.4, label = n)
  text(x = median(seq_along(n)),y = -1.6, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for only N lim effect
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
anova(mod1, mod2)

# No effects of facility, frequency, medium, pH, N source
V <- v_func(dat = dat.leafNA)
mod0 <- rma.mv(yi = rr, mods = ~  N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
mod1 <- rma.mv(yi = rr, mods = ~  facility + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~  pHcontr + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
mod3 <- rma.mv(yi = rr, mods = ~  frequencyN + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML") # many NAs
mod4 <- rma.mv(yi = rr, mods = ~  medium + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
mod5 <- rma.mv(yi = rr, mods = ~  Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="ML")
#mod6 <- rma.mv(yi = rr, mods = ~  lengthlim + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafNA, method="ML") # some NAs but NO effect

#mod7 <- rma.mv(yi = rr, mods = ~  potsize + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafNA, method="ML") # some NAs but NO effect
# Many NAs. Some evidence and more pot size is further explored below

anova(mod0, mod1) # facility
anova(mod0, mod2) # pH
anova(mod0, mod3) # N application frequency
anova(mod0, mod4) # medium
anova(mod0, mod5) # Nsource

# More on testing N source
# indicated to be important for photosynthesis and leaf area.
# Tricky here because only ammonium treatment mostly in experiments with low
# N limitation, and only in outdoor experiments.
# Degree of N limitation should be of greatest importance so we 
# control for that when testing the effect of N source.
table(X.sub.leafNA[,c("Nsource", "N_conc_exp2")])
table(X.sub.leafNA[,c("Nsource", "facility")])
V <- v_func(dat = dat.leafNA)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+Nsource, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="REML")
summary(mod1)
# Nitrate may enhance the N limitation effect, but weak evidence. 
# % decrease at 100% N limitation
1-exp(c(summary(mod1)$b[1], sum(summary(mod1)$b[c(1,3)]), sum(summary(mod1)$b[c(1,4)]) )) 

# legume effect
V <- v_func(dat = dat.leafNA)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
anova(mod1, mod2)
summary(mod1)
# Strong effect. Weaker effect of N limitation for nod-forming plants,
# and stronger effect (greater reduction of leaf N under N limitation) for legumes
# without nods, compared to non-legumes

# CO2
V <- v_func(dat = dat.leafNA)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
anova(mod1, mod2)
summary(mod1)
# elevated CO2 gives stronger effect under N limitation, but weakly significant (p=0.05)
# % decrease at 100% N limitation
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="REML")
1-exp(c(summary(mod1)$b[1], sum(summary(mod1)$b[c(1,3)]) )) 


# Crop type
V <- v_func(dat = dat.leafNA)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
anova(mod1, mod2)
summary(mod1)
#seems to be no difference between crop types

# monodicot
V <- v_func(dat = dat.leafNA)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
anova(mod1, mod2)
summary(mod1)
# nothing

# C3-C4
V <- v_func(dat = dat.leafNA)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="ML")
anova(mod1, mod2)
summary(mod1)
# An effect of this phylogenetic division
# % decrease at 100% N limitation (C3 vs C4)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNA, method="REML")
1-exp(c(summary(mod1)$b[1], sum(summary(mod1)$b[c(1,3)]))) 
#which C4 plants are there
unique(dat.leafNA[dat.leafNA$C3C4 == "C4","species"])
# only corn!!

# REML models for estimation. 
# Compare estimation of the only N lim model, with 'many covariates' model
mod1 <- rma.mv(yi = rr, mods = ~ CO2 + leg + Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="REML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="REML")
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]]) # looks good!
summary(mod1)
summary(mod2)

# % change under complete N limitation and slope of N lim effect
(exp(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub))-1)*100
# At what N limitation is there no effect?
solve(mod2$b[2], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,3], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,2], -mod2$b[1])

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.leafNA$rr)-1)*100,
                        "PropNlim" =  1-dat.leafNA$N_conc_exp2)
dat.leafNA$N_conc_exp2.rev <- 1-dat.leafNA$N_conc_exp2
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNA, method="REML")
nrow(dat.leafNA) # sample size

pred.frame <-  predict(mod1, newmods = c(seq(min(nlim.plot$PropNlim),1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.leafNA <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x=X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf N per unit area (%)", 
                      breaks = seq(-100,100,10)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimLeafNA.png", fig.Nlim.leafNA, dpi = 300)

# with MCMCglmm
# match metafor results
library(MCMCglmm)
V <- v_func(dat = dat.leafNA)
Rsvd<-svd(V)
Rsvd<-Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.leafNA$exp <- as.factor(1:NROW(dat.leafNA))
Z <- model.matrix(~exp-1, dat.leafNA)%*%Rsvd
dat.leafNA$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

m2 <- MCMCglmm(rr ~ Nsource + N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.leafNA, 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(m2)
plot(m2$Sol)


# Leaf N mass ####
#library(xlsx)
library(readxl)
#vars <- c("photosyn","chl cont","starAct", "starCont", "leaf area", "leafN_area", "leafN_mass","protein leaf", "leaf sug", "leaf starch", "SLA")

#dat.leafNM <- read.xlsx("Seufert_etal_Nlim_meta_data.xls",
#                        sheetIndex = "leafN_mass", stringsAsFactors=FALSE)
dat.leafNM <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                        sheet = "leafN_mass")
#str(dat.photo)
#merge(dat.photo, dat.leafA, by="author_species_exp")
#colnames(dat.photo)
#colnames(dat.leafA)

# remove all NA text
dat.leafNM[dat.leafNM=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.leafNM[,cols] = lapply(dat.leafNM[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.leafNM$ind.study <- factor(paste(dat.leafNM$author, dat.leafNM$Xc, sep = "_"))

# Calc effect size
dat.leafNM <- calc_resp(obj=dat.leafNM)

# order data according to ind studies
dat.leafNM <- dat.leafNM[order(dat.leafNM$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.leafNM, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.leafNM)
plot(rr ~ var.rr, dat.leafNM)
summary(lm(rr ~ var.rr, dat.leafNM))

# Explore and check covariation LEAF N MASS ####
X.sub.leafNM <- dat.leafNM[, c("rr", # RESPONSE
                               "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                               "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                               "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("leafNM_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.leafNM)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.leafNM)[1+i])
  plot(X.sub.leafNM[,"rr"] ~X.sub.leafNM[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-1.75, 0.25))
  n <- as.numeric(table(X.sub.leafNM[, 1+i]))
  nas <- sum(is.na(X.sub.leafNM[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -1.4, label = n)
  text(x = median(seq_along(n)),y = -1.6, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for only N lim effect
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
mod3 <- rma.mv(rr ~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
anova(mod1, mod2) # N lim effect
anova(mod2, mod3) # N lim quadratic effect, Evidence of non-linear response.

# No effects of facility, frequency, medium, pH, N source
V <- v_func(dat = dat.leafNM)
mod0 <- rma.mv(yi = rr, mods = ~  N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
mod1 <- rma.mv(yi = rr, mods = ~  facility + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~  pHcontr + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
#mod3 <- rma.mv(yi = rr, mods = ~  frequencyN + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafNM, method="ML") # some NAs. No effect
#mod4 <- rma.mv(yi = rr, mods = ~  medium + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafNM, method="REML") # effect!
mod5 <- rma.mv(yi = rr, mods = ~  Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML")
mod6 <- rma.mv(yi = rr, mods = ~  lengthlim + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="ML") # some NAs but NO effect
#mod7 <- rma.mv(yi = rr, mods = ~  potsize + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafNM, method="ML") # effect!!!
#Some NAs. Some evidence of an effect

anova(mod0, mod1) # facility
anova(mod0, mod2) # pH
anova(mod0, mod5) # Nsource
anova(mod0, mod6) # lengthlim. Odd effect - only the intermediate group different

# potsize and medium is correlated.
table(X.sub.leafNM[,c("medium", "potsize")])
dat.leafNM[dat.leafNM$medium == "hydroponic", "studyNr"]
# effect dominated by one study and only 3 studies with  hydroponic
# More on testing medium
table(X.sub.leafNM[,c("medium", "potsize")])
table(X.sub.leafNM[,c("medium", "facility")])
V <- v_func(dat = dat.leafNM)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+medium, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="REML")
summary(mod1)
# hydroponic medium reduces the effect of N limitation on leaf N mass.
# medium inert only n=2 so too few data points to say something
# % decrease at 100% N limitation (hydroponic, sand, soil)
1-exp(c(summary(mod1)$b[1], sum(summary(mod1)$b[c(1,4)]), sum(summary(mod1)$b[c(1,5)]) )) 

# legume effect
V <- v_func(dat = dat.leafNM)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
anova(mod1, mod2)
summary(mod1)
# Strong effect. Weaker effect of N limitation for nod-forming plants,
# and stronger effect (greater reduction of leaf N under N limitation) for legumes
# without nods, compared to non-legumes

# CO2
dat.leafNMco <- dat.leafNM[!(is.na(dat.leafNM$CO2)),]
V <- v_func(dat = dat.leafNMco)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNMco, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNMco, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect

# Crop type
V <- v_func(dat = dat.leafNM)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
anova(mod1, mod2)
summary(mod1)
#seems to be no difference between crop types

# monodicot
V <- v_func(dat = dat.leafNM)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
anova(mod1, mod2)
summary(mod1)
# nothing

# C3-C4
V <- v_func(dat = dat.leafNM)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.leafNM, method="ML")
anova(mod1, mod2)
summary(mod1)
# NO effect

# REML models for estimation. 
# Compare estimation of the only N lim model, with 'many covariates' model
mod1 <- rma.mv(yi = rr, mods = ~ medium + leg + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="REML")
mod2 <- rma.mv(rr~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="REML")
mod3 <- rma.mv(rr~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="REML")
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]]) # Not good!
plot(fitted.rma(mod3), rstandard.rma.mv(mod3)[["resid"]]) # looks good! Quadratic term better.
summary(mod1)
summary(mod2)
summary(mod3)

# % change under complete N limitation and slope of N lim effect
# need to refit with two terms to get estiamte as full N lim
dat.leafNM$quad.Nconc <- dat.leafNM$N_conc_exp2*dat.leafNM$N_conc_exp2
mod3 <- rma.mv(rr~ N_conc_exp2 + quad.Nconc, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="REML")
(exp(cbind(mod3$b,mod3$ci.lb, mod3$ci.ub))-1)*100
# At what N limitation is there no effect?
solve(mod2$b[2], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,3], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,2], -mod2$b[1])

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.leafNM$rr)-1)*100,
                        "PropNlim" =  1-dat.leafNM$N_conc_exp2)
dat.leafNM$N_conc_exp2.rev <- 1-dat.leafNM$N_conc_exp2
dat.leafNM$N_conc_exp2.rev.q <- dat.leafNM$N_conc_exp2.rev*dat.leafNM$N_conc_exp2.rev
mod1 <- rma.mv(rr ~ N_conc_exp2.rev + N_conc_exp2.rev.q , V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.leafNM, method="REML")
pred.frame <-  predict(mod1, newmods = 
                         cbind(seq(min(nlim.plot$PropNlim),1,0.001),
                               seq(min(nlim.plot$PropNlim),1,0.001)^2), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100
nrow(dat.leafNM) # sample size

#mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.leafNM, method="REML")

#pred.frame <-  predict(mod1, newmods = c(seq(min(nlim.plot$PropNlim),1,0.001)), addx=TRUE)
#pred.frame <- do.call(cbind.data.frame, pred.frame)
#pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.leafNM <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x= X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf N per unit mass (%)", 
                      breaks = seq(-100,100,10)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimleafNM.png", fig.Nlim.leafNM, dpi = 300)

# END LEAF N MASS

# Chl content ####
#library(xlsx)
library(readxl)
#vars <- c("photosyn","chl cont","starAct", "starCont", "leaf area", "leafN_area", "leafN_mass","protein leaf", "leaf sug", "leaf starch", "SLA")

dat.chl <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                        sheet = "chl cont")

# remove all NA text
dat.chl[dat.chl=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.chl[,cols] = lapply(dat.chl[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.chl$ind.study <- factor(paste(dat.chl$author, dat.chl$Xc, sep = "_"))

# Calc effect size
dat.chl <- calc_resp(obj=dat.chl)

# order data according to ind studies
dat.chl <- dat.chl[order(dat.chl$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.chl, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.chl)
plot(rr ~ var.rr, dat.chl)
summary(lm(rr ~ var.rr, dat.chl))

# Explore and check covariation LEAF CHL####
X.sub.chl <- dat.chl[, c("rr", # RESPONSE
                               "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                               "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                               "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("chl_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.chl)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.chl)[1+i])
  plot(X.sub.chl[,"rr"] ~X.sub.chl[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-1.75, 0.25))
  n <- as.numeric(table(X.sub.chl[, 1+i]))
  nas <- sum(is.na(X.sub.chl[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -1.4, label = n)
  text(x = median(seq_along(n)),y = -1.6, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for only N lim effect
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="ML")
mod3 <- rma.mv(rr ~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="ML")
anova(mod1, mod2) # N lim effect. Yes.
anova(mod2, mod3) # N lim quadratic effect. No effect.

# No effects of facility, frequency, medium, pH, N source
V <- v_func(dat = dat.chl)
mod0 <- rma.mv(yi = rr, mods = ~  N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="ML")
#mod1 <- rma.mv(yi = rr, mods = ~  facility + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.chl, method="ML")
#mod2 <- rma.mv(yi = rr, mods = ~  pHcontr + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.chl, method="ML")
#mod3 <- rma.mv(yi = rr, mods = ~  frequencyN + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.chl, method="ML") # some NAs. No effect
mod4 <- rma.mv(yi = rr, mods = ~  medium + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="ML") # effect!
#mod5 <- rma.mv(yi = rr, mods = ~  Nsource + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.chl, method="ML")
#mod6 <- rma.mv(yi = rr, mods = ~  lengthlim + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.chl, method="ML") # some NAs but NO effect
#mod7 <- rma.mv(yi = rr, mods = ~  potsize + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.chl, method="ML") # some NAs but NO effect

anova(mod0, mod4) # medium. no effect

# potsize and medium is correlated.
table(X.sub.chl[,c("medium", "potsize")])

# More on testing length and potsize
table(X.sub.chl[,c("lengthlim", "potsize")])
table(X.sub.chl[,c("leg", "potsize")])
dat.chl.len <- dat.chl[!(is.na(dat.chl$lengthlim)),]
V <- v_func(dat = dat.chl.len)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+lengthlim, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.len, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.len, method="ML")
anova(mod1,mod2) 
#not significant and odd effect (largest at intermidiate length)
dat.chl.pot <- dat.chl[!(is.na(dat.chl$potsize)),]
V <- v_func(dat = dat.chl.pot)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+potsize, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.pot, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.pot, method="ML")
anova(mod1,mod2) #significant. But many NAs (16!) AND may be an effect of legme status
summary(mod1)
table(X.sub.chl[,c("leg", "potsize")])
# Nod-forming plants have only small pots!

# legume effect
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence (P=0.07). Weaker effect of N limitation for nod-forming plants.
# lets look at leg within croptype "oilseed" separately
dat.chl.leg <- dat.chl[dat.chl$croptype=="oilseed",]
V <- v_func(dat = dat.chl.leg)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.leg, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.leg, method="ML")
anova(mod1,mod2) # Effect of legume status


# Crop type
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2 +leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence for a stronger response in cereals compared to fibre and oilseed crop.
# P=0.05
# lets look at leg within croptype "oilseed" separately
dat.chl.ct <- dat.chl[dat.chl$leg=="no",]
V <- v_func(dat = dat.chl.ct)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.ct, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.ct, method="ML")
anova(mod1,mod2) # Effect of croptype
summary(mod1)
# lets test crop type and legume together
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2 +leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod3 <- rma.mv(yi = rr, mods = ~ leg+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
anova(mod1, mod2) # legume. P=0.05
anova(mod1, mod3) # crop type. P=0.03

# CO2
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect


# monodicot
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
anova(mod1, mod2)
summary(mod1)
# Evidence for stronger response (larger reduction in chlorophyll) in monocots.
# P=0.01. BUT correlated with C4 plants (see below)

# C3-C4
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence for a stronger response for C4 plants.
# P=0.06. Two C4 species, Z. mays and S. bicolor
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl, method="REML")
summary(mod1)
# % change under complete N limitation and slope of N lim effect
(exp(mod1$b[1,1])-1)*100
(exp(sum(mod1$b[c(1,3),1]))-1)*100
# However, mono and dicot are correlated: Ther are no dicot C4 plants!
table(X.sub.chl[,c("monodicot", "leg")])
# lets look at mono-dicot separately
dat.chl.modi <- dat.chl[dat.chl$C3C4=="C3",]
V <- v_func(dat = dat.chl.modi)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.modi, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.modi, method="ML")
anova(mod1,mod2) #NO evidence of an mono-dicot effect
# and now C3-C4
dat.chl.c <- dat.chl[dat.chl$monodicot=="mono",]
V <- v_func(dat = dat.chl.c)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.c, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.chl.c, method="ML")
anova(mod1,mod2) #NO evidence of an C3-C4 effect

# REML models for estimation. 
# Compare estimation of the only N lim model, with 'many covariates' model
V <- v_func(dat = dat.chl)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg +  croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="REML")
# No effect of legume if in a model together with croptype
mod2 <- rma.mv(rr~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="REML")
plot(fitted.rma(mod1), rstandard.rma.mv(mod1)[["resid"]]) 
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]])
summary(mod1)
summary(mod2)

# % change under complete N limitation and slope of N lim effect
(exp(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub))-1)*100
# At what N limitation is there no effect?
solve(mod2$b[2], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,3], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,2], -mod2$b[1])

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.chl$rr)-1)*100,
                        "PropNlim" =  1-dat.chl$N_conc_exp2)
dat.chl$N_conc_exp2.rev <- 1-dat.chl$N_conc_exp2
nrow(dat.chl) # sample size

mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.chl, method="REML")

pred.frame <-  predict(mod1, newmods = c(seq(min(nlim.plot$PropNlim),1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.chl <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x= X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf Chl per unit area (%)", 
                      breaks = seq(-100,100,10)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimchl.png", fig.Nlim.chl, dpi = 300)

# END Chl content

# Rubisco content ####
#library(xlsx)
library(readxl)
#vars <- c("photosyn","chl cont","starAct", "starCont", "leaf area", "leafN_area", "leafN_mass","protein leaf", "leaf sug", "leaf starch", "SLA")

#dat.rub <- read.xlsx("Seufert_etal_Nlim_meta_data.xls",
#                     sheetIndex = "RubCont", stringsAsFactors=FALSE)
dat.rub <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                     sheet = "RubCont")

# remove all NA text
dat.rub[dat.rub=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.rub[,cols] = lapply(dat.rub[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.rub$ind.study <- factor(paste(dat.rub$author, dat.rub$Xc, sep = "_"))

# Calc effect size
dat.rub <- calc_resp(obj=dat.rub)

# order data according to ind studies
dat.rub <- dat.rub[order(dat.rub$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.rub, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.rub)
plot(rr ~ var.rr, dat.rub)
summary(lm(rr ~ var.rr, dat.rub))

# Explore and check covariation LEAF RUBISCO####
X.sub.rub <- dat.rub[, c("rr", # RESPONSE
                         "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                         "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                         "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("rub_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.rub)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.rub)[1+i])
  plot(X.sub.rub[,"rr"] ~X.sub.rub[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-2.5, 0.5))
  n <- as.numeric(table(X.sub.rub[, 1+i]))
  nas <- sum(is.na(X.sub.rub[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -1.7, label = n)
  text(x = median(seq_along(n)),y = -1.9, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for only N lim effect
V <- v_func(dat = dat.rub)
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
mod3 <- rma.mv(rr ~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
anova(mod1, mod2) # N lim effect. Yes....kind of. P=0.04
anova(mod2, mod3) # N lim quadratic effect. No effect.

# No effects of facility, frequency, medium, pH, 
# Indications N source, length and potsize but not significant
# and also correlated (also with legume status) and rather low sample sizes in general
V <- v_func(dat = dat.rub)
mod0 <- rma.mv(yi = rr, mods = ~  N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
mod1 <- rma.mv(yi = rr, mods = ~  facility + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~  pHcontr + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
mod3 <- rma.mv(yi = rr, mods = ~  frequencyN + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
mod4 <- rma.mv(yi = rr, mods = ~  medium + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML")
#mod5 <- rma.mv(yi = rr, mods = ~  Nsource + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.rub, method="ML") # Only nitrate lower, BUT only 2 samples in this category
#mod6 <- rma.mv(yi = rr, mods = ~  lengthlim + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML") # some NAs, 
#mod7 <- rma.mv(yi = rr, mods = ~  potsize + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="ML") # some NAs but NO effect

anova(mod0, mod1) # facility. no effect
anova(mod0, mod2) # pH. no effect
anova(mod0, mod3) # frequnency N. no effect
anova(mod0, mod4) # medium. no effect

# potsize, Nsource and length is correlated.
table(X.sub.rub[,c("lengthlim", "potsize")])
table(X.sub.rub[,c("Nsource", "potsize")])
table(X.sub.rub[,c("leg", "potsize")])

# More on testing length and potsize
table(X.sub.rub[,c("lengthlim", "potsize")])
dat.rub.len <- dat.rub[!(is.na(dat.rub$lengthlim)),]
V <- v_func(dat = dat.rub.len)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+lengthlim, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub.len, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub.len, method="ML")
anova(mod1,mod2) 
#not significant and odd effect (largest at intermidiate length)
dat.rub.nsou <- dat.rub[!(is.na(dat.rub$Nsource)),]
V <- v_func(dat = dat.rub.nsou)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+Nsource, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub.nsou, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub.nsou, method="ML")
anova(mod1,mod2) 
#not significant
dat.rub.pot <- dat.rub[!(is.na(dat.rub$potsize)),]
V <- v_func(dat = dat.rub.pot)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2+potsize, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub.pot, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub.pot, method="ML")
anova(mod1,mod2) 
#not significant


# legume effect
V <- v_func(dat = dat.rub)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence (P=0.10). Weaker effect of N limitation for nod-forming plants.
# Few studies with nod-forming plants so very low power.
# Very significant effect (P<0.0001) if not including study level random effect


# Crop type
V <- v_func(dat = dat.rub)
mod1 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect bewtween oilseed and cereal.

# CO2 ONLY 3 samples with elevated. So cant say much
V <- v_func(dat = dat.rub)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect


# monodicot
V <- v_func(dat = dat.rub)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.rub, method="ML")
anova(mod1, mod2)
summary(mod1)
# Stronger response (larger reduction in rubisco) in monocots.
# But only significant if study effect is not included (P<0.0001)
# Correlated with legume status though

# C3-C4. Only C3 plants.

# REML models for estimation. 
# Compare estimation of the only N lim model, with 'many covariates' model
V <- v_func(dat = dat.rub)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="REML")
mod2 <- rma.mv(rr~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="REML")
plot(fitted.rma(mod1), rstandard.rma.mv(mod1)[["resid"]]) 
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]])
summary(mod1)
summary(mod2)

# % change under complete N limitation and slope of N lim effect
(exp(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub))-1)*100
# At what N limitation is there no effect?
solve(mod2$b[2], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,3], -mod2$b[1])
solve(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub)[2,2], -mod2$b[1])

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.rub$rr)-1)*100,
                        "PropNlim" =  1-dat.rub$N_conc_exp2)
dat.rub$N_conc_exp2.rev <- 1-dat.rub$N_conc_exp2
nrow(dat.rub) # sample size

mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.rub, method="REML")

pred.frame <-  predict(mod1, newmods = c(seq(min(nlim.plot$PropNlim),1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.rub <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x= X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf rub per unit area (%)", 
                      breaks = seq(-100,100,10)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimrub.png", fig.Nlim.rub, dpi = 300)

# END rub content

# Leaf sugar content ####
dat.sug <- rep.var[["leafSugar"]]

# library(readxl)
# dat.sug <- read_xls("Seufert_etal_Nlim_meta_data.xls",
#                      sheet = "leaf sug")
# 
# # remove all NA text
# dat.sug[dat.sug=="NA"] <- NA
# 
# # change to factors
# # Weird with NAs in "leg" column
# cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
#          "croptype", "species", "CO2", "monodicot", "C3C4")    
# dat.sug[,cols] = lapply(dat.sug[ ,cols],  function(x) as.factor(x))
# 
# # fix data ####
# 
# # Independent studies are defined as unique controls
# dat.sug$ind.study <- factor(paste(dat.sug$author, dat.sug$Xc, sep = "_"))
# 
# # Calc effect size
# dat.sug <- calc_resp(obj=dat.sug)
# 
# # order data according to ind studies
# dat.sug <- dat.sug[order(dat.sug$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.sug, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.sug)
plot(rr ~ var.rr, dat.sug)
summary(lm(rr ~ var.rr, dat.sug))

# Explore and check covariation LEAF SUGARS ####
X.sub.sug <- dat.sug[, c("rr", # RESPONSE
                         "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                         "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                         "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("sug_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.sug)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.sug)[1+i])
  plot(X.sub.sug[,"rr"] ~X.sub.sug[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-1.5, 2.3))
  n <- as.numeric(table(X.sub.sug[, 1+i]))
  nas <- sum(is.na(X.sub.sug[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -1.1, label = n)
  text(x = median(seq_along(n)),y = -1.3, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for differences between sugar types
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + sugar, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
anova(mod1, mod2) # No effect of different sugars.


# Test for only N lim effect
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
mod3 <- rma.mv(rr ~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
anova(mod1, mod2) # N lim effect. No.
anova(mod2, mod3) # N lim quadratic effect. No effect.

# No effects of facility, frequency, medium, pH, N source, length
# Some indications of a potsize effect, but only two studies are
# behind this effect
V <- v_func(dat = dat.sug)
mod0 <- rma.mv(yi = rr, mods = ~  N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
#mod1 <- rma.mv(yi = rr, mods = ~  facility + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sug, method="ML")
#mod2 <- rma.mv(yi = rr, mods = ~  pHcontr + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sug, method="ML")
#mod3 <- rma.mv(yi = rr, mods = ~  frequencyN + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sug, method="ML") # some NAs. No effect
mod4 <- rma.mv(yi = rr, mods = ~  medium + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML") # No effect!
mod5 <- rma.mv(yi = rr, mods = ~  Nsource + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML")
mod6 <- rma.mv(yi = rr, mods = ~  lengthlim + N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="ML") # some NAs but NO effect
#mod7 <- rma.mv(yi = rr, mods = ~  potsize + N_conc_exp2, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sug, method="ML") # some indication  of an effect

anova(mod0, mod4) # medium. no effect
anova(mod0, mod5) # N source. no effect
anova(mod0, mod6) # Length of experiment. no effect

# legume effect
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug.leg, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug.leg, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect

# N fix status and croptype correlated
table(X.sub.sug[,c("croptype", "leg")])

# Crop type
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence for a stronger response in legumes compared to cereal and oilseed crop.
# P=0.07. BUT only ONE study with legumes

# CO2
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect

# monodicot
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect

# C3-C4
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sug, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect


# REML models for estimation. 
#dat.sug.sub <- dat.sug[dat.sug$rr<1,]
#dat.sug.sub <- dat.sug[!(dat.sug$studyNr == "study50"),]
#V <- v_func(dat = dat.sug.sub)
V <- v_func(dat = dat.sug)
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="REML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="REML")
hist(rstandard.rma.mv(mod1)[["resid"]]) 
summary(mod1)

# % change under N limitation (without N lim as covariate)
(exp(cbind(mod1$b,mod1$ci.lb, mod1$ci.ub))-1)*100
# % change under max N limitation (N lim as covariate)
(exp(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub))-1)*100

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.sug$rr)-1)*100,
                        "PropNlim" =  1-dat.sug$N_conc_exp2)
dat.sug$N_conc_exp2.rev <- 1-dat.sug$N_conc_exp2
nrow(dat.sug) # sample size

mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sug, method="REML")

pred.frame <-  predict(mod1, newmods = c(seq(min(nlim.plot$PropNlim),1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.sug <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x= X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf sugar (%)", 
                      breaks = seq(-100,1000,50)) +
  xlab("Proportion N limitation") +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimsug.png", fig.Nlim.sug, dpi = 300)

# END leaf sugar content

# Starch content ####
#library(xlsx)
library(readxl)

#dat.star <- read.xlsx("Seufert_etal_Nlim_meta_data.xls",
#                     sheetIndex = "leaf starch", stringsAsFactors=FALSE)
dat.star <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                      sheet = "leaf starch")

# remove all NA text
dat.star[dat.star=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.star[,cols] = lapply(dat.star[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.star$ind.study <- factor(paste(dat.star$author, dat.star$Xc, sep = "_"))

# Calc effect size
dat.star <- calc_resp(obj=dat.star)

# order data according to ind studies
dat.star <- dat.star[order(dat.star$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.star, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.star)
plot(rr ~ var.rr, dat.star)
summary(lm(rr ~ var.rr, dat.star))

# Explore and check covariation LEAF STARCH ####
X.sub.star <- dat.star[, c("rr", # RESPONSE
                         "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                         "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                         "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("starch_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.star)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.star)[1+i])
  plot(X.sub.star[,"rr"] ~X.sub.star[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-0.5, 1.5))
  n <- as.numeric(table(X.sub.star[, 1+i]))
  nas <- sum(is.na(X.sub.star[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -0.2, label = n)
  text(x = median(seq_along(n)),y = -0.4, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for only N lim effect
# Observation-level variance is close to zero,
# and study-variance is included although data are sparse
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.star, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr, ~ 1 | obsNr),
               data = dat.star, method="ML")
mod3 <- rma.mv(rr ~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.star, method="ML")
anova(mod1, mod2) # N lim effect. Yes
anova(mod2, mod3) # N lim quadratic effect. No

# Sample size so small that there is no point testing effects of different covariates. 
# The overall figure do not indicate any obvious effects that cannot be 
# attributed to one single study.
# potsize, Nsource and length is correlated.
table(X.sub.star[,c("lengthlim", "potsize")])
table(X.sub.star[,c("Nsource", "potsize")])
table(X.sub.star[,c("leg", "potsize")])

# legume effect
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
anova(mod1, mod2)
summary(mod1)
# No. But no studies with nof-forming plants

# Crop type
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ croptype+N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
anova(mod1, mod2)
summary(mod1)
# Legume effect but only 2 data points so doesnt mean anything

# CO2
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect


# monodicot
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effects

# C3-C4
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.star, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect

# REML models for estimation. 
V <- v_func(dat = dat.star)
mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.star, method="REML")
mod2 <- rma.mv(rr~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.star, method="REML")
plot(fitted.rma(mod1), rstandard.rma.mv(mod1)[["resid"]]) 
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]])
summary(mod1)
summary(mod2)

# % change under complete N limitation and slope of N lim effect (absolute v% change)
Nlim_eff(mod1)
#realtive % change
exp(cbind(mod1$b,mod1$ci.lb, mod1$ci.ub)[2,])*100

# Negative effect of N lim so at what N limitation level
# is there no effect?
solve(mod1$b[2], -mod1$b[1]) # at N lim 1.11
solve(cbind(mod1$b,mod1$ci.lb, mod1$ci.ub)[2,3], -mod1$b[1])
solve(cbind(mod1$b,mod1$ci.lb, mod1$ci.ub)[2,2], -mod1$b[1])
# But CI includes 1 (lower at N lim 0.76)

# plot overall N lim effect
# make data frame for plotting first
nlim.plot <- data.frame("Neffect" = (exp(dat.star$rr)-1)*100,
                        "PropNlim" =  1-dat.star$N_conc_exp2)
dat.star$N_conc_exp2.rev <- 1-dat.star$N_conc_exp2
nrow(dat.star) # sample size

mod1 <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.star, method="REML")

pred.frame <-  predict(mod1, newmods = c(seq(min(nlim.plot$PropNlim),1,0.001)), addx=TRUE)
pred.frame <- do.call(cbind.data.frame, pred.frame)
pred.frame$pred <- (exp(pred.frame$pred)-1)*100

library(ggplot2)
fig.Nlim.star <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
  geom_point() +
  geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev)) +
  geom_ribbon(data=pred.frame,aes(x= X.N_conc_exp2.rev, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
              inherit.aes=FALSE, alpha=0.3) +
  geom_hline(yintercept = 0, lty=2) +
  scale_y_continuous( name = "Change in leaf starch (%)", 
                      breaks = seq(-25,600,50)) +
  xlab("Proportion N limitation") +
  xlim(c(0,1))+
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12))
ggsave("nlimstarch.png", fig.Nlim.star, dpi = 300)

# END starch content


# SLA ####
#library(xlsx)
library(readxl)
#vars <- c("photosyn","chl cont","starAct", "starCont", "leaf area", "leafN_area", "leafN_mass","protein leaf", "leaf sug", "leaf starch", "SLA")

#dat.sla <- read.xlsx("Seufert_etal_Nlim_meta_data.xls",
#                        sheetIndex = "SLA", stringsAsFactors=FALSE)
dat.sla <- read_xls("Seufert_etal_Nlim_meta_data.xls",
                     sheet = "SLA")

#str(dat.photo)
#merge(dat.photo, dat.leafA, by="author_species_exp")
#colnames(dat.photo)
#colnames(dat.leafA)

# remove all NA text
dat.sla[dat.sla=="NA"] <- NA

# change to factors
# Weird with NAs in "leg" column
cols = c("Nsource", "frequencyN", "potsize", "facility", "medium", "leg", "lengthlim", "pHcontr",
         "croptype", "species", "CO2", "monodicot", "C3C4")    
dat.sla[,cols] = lapply(dat.sla[ ,cols],  function(x) as.factor(x))

# fix data ####

# Independent studies are defined as unique controls
dat.sla$ind.study <- factor(paste(dat.sla$author, dat.sla$Xc, sep = "_"))

# Calc effect size
dat.sla <- calc_resp(obj=dat.sla)

# order data according to ind studies
dat.sla <- dat.sla[order(dat.sla$ind.study),]


# Fix known variance-covariance matrix
V <- v_func(dat = dat.sla, ind.study = "ind.study", 
            var.control ="var.control", var.rr = "var.rr")
#plot rr against N and variance
plot(rr ~ Ne, dat.sla)
plot(rr ~ var.rr, dat.sla)
summary(lm(rr ~ var.rr, dat.sla))

# Explore and check covariation SLA ####
X.sub.sla <- dat.sla[, c("rr", # RESPONSE
                               "N_conc_exp2", "lengthlim", "frequencyN", "Nsource", # Possible Treatment effects 
                               "species", "leg", "monodicot", "C3C4","croptype", # taxonomy aggregation 
                               "CO2", "pHcontr",  "facility", "medium", "potsize" # Possible experimental set up effects
)]
# PLot matrix of the predictors
png("sla_RR_all.png", width = 18, height = 30, units ="cm", res=600)
par(mfrow = c(5,3))
for (i in 1:(ncol(X.sub.sla)-1) ) {
  tit <- ifelse(i==1, "N lim (0=max limitation)",colnames(X.sub.sla)[1+i])
  plot(X.sub.sla[,"rr"] ~X.sub.sla[,1+i], ylab = "log RR", xlab="", 
       main = tit, las=2, ylim = c(-1.75, 0.25))
  n <- as.numeric(table(X.sub.sla[, 1+i]))
  nas <- sum(is.na(X.sub.sla[, 1+i]))
  #nas = sum(X.sub[, 1+i] == "NA")
  text(x = seq_along(n),y = -1.4, label = n)
  text(x = median(seq_along(n)),y = -1.6, label = paste("NAs=",nas), font=2)
}
dev.off()

# Test for only N lim effect
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML")
mod3 <- rma.mv(rr ~ poly(N_conc_exp2,2), V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML")
anova(mod1, mod2) # N lim effect
anova(mod2, mod3) # N lim quadratic effect
# NO effect of N limitation

# Is there an overall effect?
mod1 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="REML")
summary(mod1)
# No

# No effects of facility, frequency, medium, pH, N source
V <- v_func(dat = dat.sla)
mod0 <- rma.mv(yi = rr, mods = ~  1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML")
#mod1 <- rma.mv(yi = rr, mods = ~  facility, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sla, method="ML") # some NAs. No effect.
mod2 <- rma.mv(yi = rr, mods = ~  pHcontr, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML")
#mod3 <- rma.mv(yi = rr, mods = ~  frequencyN, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sla, method="ML") # some NAs. No effect
mod4 <- rma.mv(yi = rr, mods = ~  medium, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML")
#mod5 <- rma.mv(yi = rr, mods = ~  Nsource, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sla, method="ML")  # some NAs. No effect
mod6 <- rma.mv(yi = rr, mods = ~  lengthlim, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="ML") # some NAs but NO effect
#mod7 <- rma.mv(yi = rr, mods = ~  potsize, V = V, intercept=TRUE, 
#               random = list(~ 1 | studyNr,~ 1 | obsNr),
#               data = dat.sla, method="ML") # some NAs but NO effect


anova(mod0, mod2) # pH. Only 3 data points for pH control so not 
anova(mod0, mod4) # medium. Only one data point for inert. If removed no effect.
anova(mod0, mod6) # lengthlim. No effect.
# No obvious eveidence of large effects.

# legume effect
V <- v_func(dat = dat.sla)
mod1 <- rma.mv(yi = rr, mods = ~ leg, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence. Weaker effect of N limitation for nod-forming plants,
# and stronger effect (greater reduction of SLA under N limitation) for legumes
# without nods, compared to non-legumes

# CO2
dat.slaco <- dat.sla[!(is.na(dat.sla$CO2)),]
V <- v_func(dat = dat.slaco)
mod1 <- rma.mv(yi = rr, mods = ~ CO2, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.slaco, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.slaco, method="ML")
anova(mod1, mod2)
summary(mod1)
# No effect

# Crop type
V <- v_func(dat = dat.sla)
mod1 <- rma.mv(yi = rr, mods = ~ croptype, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
anova(mod1, mod2)
summary(mod1)
#seems to be no difference between crop types

# monodicot
V <- v_func(dat = dat.sla)
mod1 <- rma.mv(yi = rr, mods = ~ monodicot, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
anova(mod1, mod2)
summary(mod1)
# nothing

# C3-C4
V <- v_func(dat = dat.sla)
mod1 <- rma.mv(yi = rr, mods = ~ C3C4, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
mod2 <- rma.mv(yi = rr, mods = ~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | obsNr, ~ 1 | studyNr),
               data = dat.sla, method="ML")
anova(mod1, mod2)
summary(mod1)
# Some evidence that C4 plants increase SLA under N limitation, compared to
# C3 plants. But only 4 data points for C4 plants not reliable.

# REML models for estimation. 
# Compare estimation of the only N lim model, with 'many covariates' model
mod1 <- rma.mv(yi = rr, mods = ~ leg, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="REML")
mod2 <- rma.mv(rr~ 1, V = V, intercept=TRUE, 
               random = list(~ 1 | studyNr,~ 1 | obsNr),
               data = dat.sla, method="REML")
plot(fitted.rma(mod1), rstandard.rma.mv(mod1)[["resid"]])
plot(fitted.rma(mod2), rstandard.rma.mv(mod2)[["resid"]]) 
summary(mod1)
summary(mod2)

# % change under N limitation for legumes
(exp(mod1$b[1,1])-1)*100 # Not legumes
(exp(sum(mod1$b[1:2,1]))-1)*100 # Legumes with nods
(exp(sum(mod1$b[c(1,3),1]))-1)*100 # Legumes without nods

# % change under N limitation
(exp(cbind(mod2$b,mod2$ci.lb, mod2$ci.ub))-1)*100

# END SLA

# Figure 6 in ms - photoVSleaf y NleafVSleaf ####
# First get photo and leaf data for figure (a)
leaf.sub <- rep.var["leaf_area"]$leaf_area[,c("author_species_exp", "rr", "var.rr", "var.control")]
#leaf.sub <-dat.leafA[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leaf", "var.rr.leaf", "var.control.leaf")
comb.dat.a <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
#comb.dat.a <- merge(dat.photo, leaf.sub, by="author_species_exp")
comb.dat.a <- droplevels(comb.dat.a)

# Second get N leaf per area and leaf data for figure (b)
comb.dat.b <- merge(rep.var["leafN_area"]$leafN_area, leaf.sub, by="author_species_exp")
#comb.dat.a <- merge(dat.photo, leaf.sub, by="author_species_exp")
comb.dat.b <- droplevels(comb.dat.b)

# plot figure 6
p.col.a <- col.df[col.df$col.species %in% comb.dat.a$species,1] 
fig6a <- ggplot(comb.dat.a, aes(x=rr.leaf, y=rr, group = species)) +
            geom_point(aes(shape=leg, color = species), size=4) +
            geom_abline(slope=1) +
            geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
            geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
            xlim(c(-3.0, 0.5)) +
            ylim(c(-3.0, 0.5)) +
            ylab(expression(paste("Change in photosynthesis (log"[e],italic(rr), ")"))) +
            xlab(expression(paste("Change in leaf area (log"[e],italic(rr), ")"))) +
            scale_color_manual(name = "Species", values = p.col.a, 
                     labels = c(expression(italic("G. hirsutum")), 
                                expression(italic("G. max")),
                                expression(italic("O. sativa")), 
                                expression(italic("S. bicolor")), 
                                expression(italic("T. aestivum")), 
                                expression(italic("Z. mays")) )) +
            scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
            theme(
              axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
              axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
              axis.title.x = element_text(size=14, vjust=5),
              axis.title.y = element_text(size=14, vjust=-5),
              axis.line.x = element_line(color="black", size = 1),
              axis.line.y = element_line(color="black", size = 1),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              legend.text.align = 0,
              legend.text = element_text(size=14),
              legend.title = element_text(size=14),
              axis.ticks.length=unit(-0.20, "cm")  ) +
              annotation_custom(
                grob = textGrob(label = "(a)", gp = gpar(fontsize = 23)),
                ymin = 0.3,      # Vertical position of the textGrob
                #  ymax = 5,
                xmin = -6.5)


p.col.b <- col.df[col.df$col.species %in% comb.dat.b$species,1] 
fig6b <- ggplot(comb.dat.b, aes(x=rr.leaf, y=rr, group = species)) +
                      geom_point(aes(shape=leg, color = species), size=4) +
                      geom_abline(slope=1) +
                      geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
                      geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
                      xlim(c(-3.0, 0.5)) +
                      ylim(c(-3.0, 0.5)) +
                      ylab(expression(paste("Change in ",N[L]," per unit area (log"[e],italic(rr), ")"))) +
                      xlab(expression(paste("Change in leaf area (log"[e],italic(rr), ")"))) +
                      scale_color_manual(name = "Species", values = p.col.b, 
                                         labels = c(expression(italic("G. hirsutum")), 
                                                    expression(italic("G. max")),
                                                    expression(italic("H. vulgare")),                                                     
                                                    expression(italic("O. sativa")), 
                                                    expression(italic("T. aestivum")))) +
                      scale_shape_discrete(name = "Legume status",
                          breaks=c("no","yes, non nod","yes, nod"),
                          labels = c("not legume", "legume, no nods", "legume with nods")) +
                      theme(
                        axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                        axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                        axis.title.x = element_text(size=14, vjust=5),
                        axis.title.y = element_text(size=14, vjust=-5),
                        axis.line.x = element_line(color="black", size = 1),
                        axis.line.y = element_line(color="black", size = 1),
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black", fill=NA, size=1),
                        panel.background = element_blank(),
                        legend.text.align = 0,
                        legend.text = element_text(size=14),
                        legend.title = element_text(size=14),
                        axis.ticks.length=unit(-0.20, "cm")  ) +
                        guides(color = guide_legend(order=2),
                               shape = guide_legend(order=1)) +
                        annotation_custom(
                          grob = textGrob(label = "(b)", gp = gpar(fontsize = 23)),
                          ymin = 0.3,      # Vertical position of the textGrob
                          #  ymax = 5,
                          xmin = -6.5)
              
png("figure6ab.png", width=18, height=24, units="cm", res=300)
grid.arrange(fig6a,fig6b, ncol=1, nrow =2)
dev.off()

# Figure 7 in ms - photoVSleafNarea y photoVSleafNmass ####
# First get photo and leaf N area  data for figure (a)
leaf.sub <- rep.var["leafN_area"]$leafN_area[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat.a <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat.a <- droplevels(comb.dat.a)

# Second get N leaf per area and leaf data for figure (b)
leaf.sub <- rep.var["leafN_mass"]$leafN_mass[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat.b <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat.b <- droplevels(comb.dat.b)

# plot figure 7
p.col.a <- col.df[col.df$col.species %in% comb.dat.a$species,1] 
fig7a <- ggplot(comb.dat.a, aes(x=rr.leafN, y=rr, group = species)) +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  xlim(c(-1.75, 0.5)) +
  ylim(c(-2.75, 0.5))+
  ylab(expression(paste("Change in photosynthesis (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in leaf N per unit area (log"[e],italic(rr), ")"))) +
  scale_color_manual(name = "Species", values = p.col.a, 
                     labels = c(expression(italic("B. napus")), 
                                expression(italic("G. hirsutum")),
                                expression(italic("G. max")), 
                                expression(italic("H. vulgare")), 
                                expression(italic("O. sativa")), 
                                expression(italic("P. vulgaris")),
                                expression(italic("T. aestivum")),
                                expression(italic("Z. mays")) )) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  theme(
    axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x = element_text(size=14, vjust=5),
    axis.title.y = element_text(size=14, vjust=-5),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    legend.text.align = 0,
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.ticks.length=unit(-0.20, "cm")  )  +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1)) +
  annotation_custom(
    grob = textGrob(label = "(a)", gp = gpar(fontsize = 23)),
    ymin = 0.3,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -4)
   

p.col.b <- col.df[col.df$col.species %in% comb.dat.b$species,1] 
fig7b <- ggplot(comb.dat.b, aes(x=rr.leafN, y=rr, group = species)) +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  xlim(c(-1.75, 0.5)) +
  ylim(c(-2.75, 0.5))+
  ylab(expression(paste("Change in photosynthesis (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in leaf N per unit mass (log"[e],italic(rr), ")"))) +
  scale_color_manual(name = "Species", values = p.col.b, 
                     labels = c(expression(italic("G. hirsutum")), 
                                expression(italic("G. max")),
                                expression(italic("M. esculenta")),                                                     
                                expression(italic("O. sativa")), 
                                expression(italic("T. aestivum")),
                                expression(italic("Z. mays")) )) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  theme(
    axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x = element_text(size=14, vjust=5),
    axis.title.y = element_text(size=14, vjust=-5),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    legend.text.align = 0,
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.ticks.length=unit(-0.20, "cm")  ) +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1)) +
  annotation_custom(
    grob = textGrob(label = "(b)", gp = gpar(fontsize = 23)),
    ymin = 0.3,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -4)

png("figure7ab.png", width=36, height=13, units="cm", res=300)
grid.arrange(fig7a,fig7b, ncol=2, nrow =1)
dev.off()


# Figure 8 in ms - photo, Chl, Rubisco ####
# First get leaf N area  and photo data for figure (a)
leaf.sub <- rep.var["leafChl"]$leafChl[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat.a <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat.a <- droplevels(comb.dat.a)

# Second get leafRubisco and photo for figure (b)
leaf.sub <- rep.var["leafRubisco"]$leafRubisco[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat.b <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat.b <- droplevels(comb.dat.b)

# Second get Chl and Rubisco for figure (b)
leaf.sub <- rep.var["leafRubisco"]$leafRubisco[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat.c <- merge(rep.var["leafChl"]$leafChl, leaf.sub, by="author_species_exp")
comb.dat.c <- droplevels(comb.dat.c)

# plot figure 8
p.col.a <- col.df[col.df$col.species %in% comb.dat.a$species,1] 
fig8a <- ggplot(comb.dat.a, aes(x=rr.leafN, y=rr, group = species)) +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  ylab(expression(paste("Change in photosynthesis (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in leaf chl per unit area (log"[e],italic(rr), ")"))) +
  xlim(c(-2.5, 0.5)) +
  ylim(c(-3.0, 0.5)) +  
  scale_color_manual(name = "Species", values = p.col.a, 
                     labels = c(expression(italic("B. napus")), 
                                expression(italic("G. hirsutum")),
                                expression(italic("G. max")), 
                                expression(italic("O. sativa")), 
                                expression(italic("T. aestivum")),
                                expression(italic("Z. mays")) )) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  theme(
    axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x = element_text(size=14, vjust=5),
    axis.title.y = element_text(size=14, vjust=-5),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    legend.text.align = 0,
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.ticks.length=unit(-0.20, "cm")  )  +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1)) +
  annotation_custom(
    grob = textGrob(label = "(a)", gp = gpar(fontsize = 23)),
    ymin = 0.3,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -5.5)


p.col.b <- col.df[col.df$col.species %in% comb.dat.b$species,1] 
fig8b <- ggplot(comb.dat.b, aes(x=rr.leafN, y=rr, group = species)) +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  ylab(expression(paste("Change in photosynthesis (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in leaf Rubisco per unit area (log"[e],italic(rr), ")"))) +
  xlim(c(-2.5, 0.5)) +
  ylim(c(-3.0, 0.5)) +  
  scale_color_manual(name = "Species", values = p.col.b, 
                     labels = c(expression(italic("G. max")),
                                expression(italic("O. sativa")), 
                                expression(italic("T. aestivum")) )) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  theme(
    axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x = element_text(size=14, vjust=5),
    axis.title.y = element_text(size=14, vjust=-5),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    legend.text.align = 0,
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.ticks.length=unit(-0.20, "cm")  ) +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1)) +
  annotation_custom(
    grob = textGrob(label = "(b)", gp = gpar(fontsize = 23)),
    ymin = 0.3,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -5.5)

p.col.c <- col.df[col.df$col.species %in% comb.dat.c$species,1] 
fig8c <- ggplot(comb.dat.c, aes(x=rr.leafN, y=rr, group = species)) +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, color = "gray") +
  ylab(expression(paste("Change in leaf chl per unit area (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in Rubisco per unit area (log"[e],italic(rr), ")"))) +
  ylim(c(-3.0, 0.5)) +
  xlim(c(-2.5, 0.5)) +
  scale_color_manual(name = "Species", values = p.col.b, 
                     labels = c(expression(italic("G. max")),
                                expression(italic("O. sativa")), 
                                expression(italic("T. aestivum")) )) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  theme(
    axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x = element_text(size=14, vjust=5),
    axis.title.y = element_text(size=14, vjust=-5),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    legend.text.align = 0,
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
    axis.ticks.length=unit(-0.20, "cm")  ) +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=1)) +
  annotation_custom(
    grob = textGrob(label = "(c)", gp = gpar(fontsize = 23)),
    ymin = 0.3,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -5.5)

png("figure8abc.png", width=18, height=39, units="cm", res=300)
grid.arrange(fig8a, fig8b, fig8c, ncol=1, nrow =3)
dev.off()



# photo vs leaf area ####
leaf.sub <- rep.var["leaf_area"]$leaf_area[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leaf", "var.rr.leaf", "var.control.leaf")
comb.dat <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)

#leaf.sub <-dat.leafA[,c("author_species_exp", "rr", "var.rr", "var.control")]
#colnames(leaf.sub) <- c("author_species_exp", "rr.leaf", "var.rr.leaf", "var.control.leaf")
#comb.dat <- merge(dat.photo, leaf.sub, by="author_species_exp")
#comb.dat <- droplevels(comb.dat)
#head(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                  # ID variables - all the variables to keep but not split apart on
                  #id.vars=c("rr", "rr.leaf"),
                  # The source columns
                  measure.vars=c("rr", "rr.leaf"),
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="resp",
                  value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leaf"),
     variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leaf"),
                                  variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb","Nsource")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 
#comb.dat.long <-  comb.dat.long[!(comb.dat.long$Nsource == "nit"),]

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")

# test with correlation between photosynthesis and leaf area
# within studies
cor.r = -0.4

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1))

vcov(mod1, type = "obs")


#extra leaf A
V <- v_func(dat = dat.leafA)
Rsvd<-svd(V)
Rsvd<-Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.leafA$exp <- as.factor(1:NROW(dat.leafA))
Z <- model.matrix(~exp-1, dat.leafA)%*%Rsvd
dat.leafA$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

m2 <- MCMCglmm(rr ~ Nsource + N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.leafA[!(is.na(dat.leafA$Nsource)),], 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(m2)
plot(m2$Sol)


###
plot(rr ~ rr.leaf, comb.dat)
abline(a=0,b=1)
points(rr ~ rr.leaf, comb.dat, col= as.numeric(comb.dat$species))

# Leaf N area vs leaf area ####
# First get photo and leaf data for figure (a)
leaf.sub <- rep.var["leaf_area"]$leaf_area[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leaf", "var.rr.leaf", "var.control.leaf")
comb.dat <- merge(rep.var["leafN_area"]$leafN_area, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
nrow(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leaf"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leaf"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leaf"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")

# test with correlation between N leaf per area and leaf area
# within studies
cor.r = -0.4

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1))

vcov(mod1, type = "obs")


#extra with MCMCglmm
V <- v_func(dat = dat.leafA)
Rsvd<-svd(V)
Rsvd<-Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.leafA$exp <- as.factor(1:NROW(dat.leafA))
Z <- model.matrix(~exp-1, dat.leafA)%*%Rsvd
dat.leafA$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

m2 <- MCMCglmm(rr ~ Nsource + N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.leafA[!(is.na(dat.leafA$Nsource)),], 
               prior=prior,nitt=100000,burnin=15000, pr=TRUE)
summary(m2)
plot(m2$Sol)


# photo vs leaf N area ####
leaf.sub <- rep.var["leafN_area"]$leafN_area[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
nrow(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")

# test with correlation between photosynthesis and leaf area
# within studies
cor.r = -0.35

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1))

varMat <- vcov(mod1, type = "obs")
mod1$X %*% coef(mod1) + varMat

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")


library(MCMCglmm)
V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control.comb", var.rr="var.rr.comb",
            blockList = FALSE)

C<-list(S1=matrix(c(1,0.5, 0.5,1),2,2), S2=matrix(c(2,0, 0,2),2,2))
Cinv_sparse<-as(list2bdiag(lapply(C,solve)), "dgCMatrix")
#Assuming the data are sorted response type with group (i.e. in the same 
#                                                       order as C) then you can fit
comb.dat.long$observation<-1:nrow(comb.dat.long)
rownames(Cinv_sparse)<-1:nrow(comb.dat.long)
m.biv <- MCMCglmm::MCMCglmm(rr.comb ~ resp-1, random=~observation, 
                            rcov=~us(resp):mv.ind.stud, data=comb.dat.long, 
                            ginverse=list(observation=Cinv_sparse), 
                            prior=list(R=list(V=diag(2), nu=0), G=list(G1=list(V=1, fix=1))))

Rsvd<-svd(V)
Rsvd<-Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
comb.dat.long$exp <- as.factor(1:NROW(comb.dat.long))
Z <- model.matrix(~exp-1, comb.dat.long)%*%Rsvd
comb.dat.long$Z <- Z
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

m2 <- MCMCglmm(rr.comb ~ resp-1,  
               random=~idv(Z), rcov=~us(resp):units, 
               data=comb.dat.long, 
               nitt=100000,burnin=15000, pr=TRUE)
summary(m2)
plot(m2$Sol)

# END PHOTO vs LEAF A

# photo vs leaf N mass ####
leaf.sub <- rep.var["leafN_mass"]$leafN_mass[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
nrow(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")

# test with correlation between photosynthesis and leaf area
# within studies
cor.r = -0.4

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)


mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1))
# There are two outliers. Very low photosynthesis.
# However, results almost the identlical if they are excluded.


varMat <- vcov(mod1, type = "obs")
mod1$X %*% coef(mod1) + varMat

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")


# photo vs leaf chl area ####
leaf.sub <- rep.var["leafChl"]$leafChl[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
nrow(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")

# test with correlation between photosynthesis and leaf area
# within studies
cor.r = 0.35

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)


mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# No difference in mean response (P=0.10)
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1)) # strong correlation
# There are one outlier. Very low photosynthesis.
# However, results almost the identlical if they are excluded.


varMat <- vcov(mod1, type = "obs")
mod1$X %*% coef(mod1) + varMat

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")

# End photo vs leaf chl area


# photo vs starch ####
leaf.sub <-dat.star[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(dat.photo, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
head(comb.dat)
nrow(comb.dat)

PvsLstararea <- ggplot(comb.dat, aes(x=rr.leafN, y=rr, group = species)) +
  #geom_smooth() +
  geom_point(aes(shape=leg, color = species), size=4) +
  #geom_abline(slope=1) +
  ylab("RR photosynthesis") +
  xlab("RR leaf starch") +
  xlim(c(-0.3, 1.5)) +
  ylim(c(-1.5, 0.3)) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods"))
ggsave("PvsLstarch_area.png", PvsLstararea, dpi = 300, width = 15, height=12, units = "cm")  

# multivariate corr test
library(reshape2)
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
#colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control.comb", var.rr="var.rr.comb")

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# Difference in mean response (P=0.0001)
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1)) # strong correlation
# There is one outlier. Very strong photosynthesis effect but not as much for starisco.

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")

# End photo vs starch

# photo vs rubisco ####
leaf.sub <- rep.var["leafRubisco"]$leafRubisco[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(rep.var["photosynthesis"]$photosynthesis, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
nrow(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")

# test with correlation between photosynthesis and leaf area
# within studies
cor.r = 0.4

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# Difference in mean response (P=0.0001)
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1)) # strong correlation
# There is one outlier. Very strong photosynthesis effect but not as much for rubisco.

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")

# End photo vs rubisco

# Chl leaf area vs leaf rubisco ####
leaf.sub <- rep.var["leafRubisco"]$leafRubisco[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(rep.var["leafChl"]$leafChl, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
nrow(comb.dat)

# multivariate corr test
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control", var.rr="var.rr")
# test with correlation between photosynthesis and leaf area
# within studies
cor.r = 0.4

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# No difference in mean response (P=0.10)
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1)) # strong correlation

varMat <- vcov(mod1, type = "obs")
mod1$X %*% coef(mod1) + varMat

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")

# End chl vs rubisco


# Chl leaf area vs leaf starch ####
leaf.sub <-dat.star[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(dat.chl, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
head(comb.dat)
nrow(comb.dat)

ChlvsLstar <- ggplot(comb.dat, aes(x=rr.leafN, y=rr, group = species)) +
  #geom_smooth() +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  ylab("RR leaf Chl area") +
  xlab("RR leaf starch") +
  xlim(c(-2.5, 0.75)) +
  ylim(c(-2, 0.75)) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods"))
ggsave("ChlvsLstar_area.png", ChlvsLstar, dpi = 300, width = 15, height=12, units = "cm")  

# multivariate corr test
library(reshape2)
comb.dat.long <- melt(comb.dat,
                      # ID variables - all the variables to keep but not split apart on
                      #id.vars=c("rr", "rr.leaf"),
                      # The source columns
                      measure.vars=c("rr", "rr.leafN"),
                      # Name of the destination column that will identify the original
                      # column that the measurement came from
                      variable.name="resp",
                      value.name="rr.comb")

comb.dat.long$var.rr.comb <- melt(comb.dat, measure.vars=c("var.rr", "var.rr.leafN"),
                                  variable.name="resp", value.name="var.rr")$var.rr
comb.dat.long$var.control.comb <- melt(comb.dat, measure.vars=c("var.control", "var.control.leafN"),
                                       variable.name="resp", value.name="var.control.comb")$var.control.comb
comb.dat.long <- comb.dat.long[order(comb.dat.long$ind.study),]

comb.dat.long <- comb.dat.long[, c("ind.study", "studyNr", "resp", "rr.comb", "var.rr.comb", "var.control.comb")]

comb.dat.long$mv.ind.stud <- factor(paste(comb.dat.long$ind.study, comb.dat.long$resp, sep = "_") )
#comb.dat.long$obsNr <- factor(paste(comb.dat.long$obsNr, comb.dat.long$resp, sep = "_") )
#colnames(comb.dat.long)[5:6] <-c("var.rr", "var.control") 

V <- v_func(dat = comb.dat.long, ind.study="mv.ind.stud", 
            var.control="var.control.comb", var.rr="var.rr.comb")

mod1 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="UN",
               data = comb.dat.long, method="REML")
mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study, ~ 1 | studyNr), struct="DIAG",
               data = comb.dat.long, method="REML")

summary(mod1) # full model with correlation between the responses
# No difference in mean response (P=0.10)
# test covariance
chisq <- as.numeric(2*(logLik(mod1)-logLik(mod2)))
0.5*(1-pchisq(chisq, 1)) # strong correlation

varMat <- vcov(mod1, type = "obs")
mod1$X %*% coef(mod1) + varMat

# test with correlation
cor.r = -0.2

covaris = list()
for (i in 1:nlevels(comb.dat.long$ind.study)) {
  studid <- comb.dat.long$ind.study == levels(comb.dat.long$ind.study)[i]
  vari <- comb.dat.long[ studid, "var.rr.comb"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

mod2 <- rma.mv(yi = rr.comb, mods = ~ resp, V = V, intercept=TRUE, 
               #random = list(~ 1 | obsNr, ~ 1 | studyNr),
               random = list(~ resp | ind.study), struct="UN",
               data = comb.dat.long, method="REML")
summary(mod1)
vcov(mod1, type = "obs")

# End chl vs star

# N leaf mass vs leaf rubisco ####
leaf.sub <-dat.rub[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(dat.leafNM, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
head(comb.dat)
nrow(comb.dat)

all(comb.dat$species %in% col.df$col.species) # check if all species are in the color list
p.col <- col.df[col.df$col.species %in% comb.dat$species,1]  
NMleafvsLrub <- ggplot(comb.dat, aes(x=rr.leafN, y=rr, group = species)) +
  #geom_smooth() +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  ylab(expression(paste("Change in leaf N per unit mass (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in rubisco per unit area (log"[e],italic(rr), ")"))) +
  ylim(c(-2.5, 0.5)) +
  xlim(c(-2.5, 0.5)) +
  scale_color_manual(name = "Species", values = p.col) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  annotation_custom(
    grob = textGrob(label = "a)", gp = gpar(fontsize = 20)),
    ymin = 0.25,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -5.5)
NMleafvsLrub.gt <- ggplot_gtable(ggplot_build(NMleafvsLrub))
NMleafvsLrub.gt$layout$clip[NMleafvsLrub.gt$layout$name == "panel"] <- "off"
grid.draw(NMleafvsLrub.gt)
# end 

# N leaf mass vs leaf chl ####
leaf.sub <-dat.chl[,c("author_species_exp", "rr", "var.rr", "var.control")]
colnames(leaf.sub) <- c("author_species_exp", "rr.leafN", "var.rr.leafN", "var.control.leafN")
comb.dat <- merge(dat.leafNM, leaf.sub, by="author_species_exp")
comb.dat <- droplevels(comb.dat)
head(comb.dat)
nrow(comb.dat)

all(comb.dat$species %in% col.df$col.species) # check if all species are in the color list
p.col <- col.df[col.df$col.species %in% comb.dat$species,1]  
NMleafvsLchl <- ggplot(comb.dat, aes(x=rr.leafN, y=rr, group = species)) +
  #geom_smooth() +
  geom_point(aes(shape=leg, color = species), size=4) +
  geom_abline(slope=1) +
  ylab(expression(paste("Change in leaf N per unit mass (log"[e],italic(rr), ")"))) +
  xlab(expression(paste("Change in leaf chl per unit area (log"[e],italic(rr), ")"))) +
  ylim(c(-2.5, 0.5)) +
  xlim(c(-2.5, 0.5)) +
  scale_color_manual(name = "Species", values = p.col) +
  scale_shape_discrete(name = "Legume status",
                       breaks=c("no","yes, non nod","yes, nod"),
                       labels = c("not legume", "legume, no nods", "legume with nods")) +
  annotation_custom(
    grob = textGrob(label = "b)", gp = gpar(fontsize = 20)),
    ymin = 0.25,      # Vertical position of the textGrob
    #  ymax = 5,
    xmin = -5.5)
NMleafvsLchl.gt <- ggplot_gtable(ggplot_build(NMleafvsLchl))
NMleafvsLchl.gt$layout$clip[NMleafvsLchl.gt$layout$name == "panel"] <- "off"
grid.draw(NMleafvsLchl.gt)

# xtra figure
png("figure_NvdRub_Chl.png", width=32, height=12, units="cm", res=300)
grid.arrange(NMleafvsLrub.gt, NMleafvsLchl.gt, ncol=2, nrow =1)
dev.off()

# end 

# Overview: full N lim. Fig 3 in MS ####

# list of data frames. 'rep.var' already loaded in the begining.
  # rep.var <- list(photosynthesis=dat.photo, leaf_area=dat.leafA, leafN_area = dat.leafNA, 
  #                leafN_mass = dat.leafNM, leafChl = dat.chl, leafRubisco = dat.rub, 
  #                sla = dat.sla, leafStarch = dat.star, leafSugar = dat.sug)

# add response names as a new column
rep.names <- as.list(setNames(names(rep.var), names(rep.var)))
rep.var <- Map(cbind, rep.var, response = rep.names)

# run model for all responses and save output of the effect
res <- lapply(rep.var, function (x) {
              # var-covar matrix
              V <- v_func(dat = x, ind.study = "ind.study", 
                   var.control ="var.control", var.rr = "var.rr")
              # model. Leaf N on mass needs a quadratic term
              if(x$response[1] == "leafN_mass") {
                  x$quad.Nconc <- x$N_conc_exp2*x$N_conc_exp2
                  mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + quad.Nconc, V = V, intercept=TRUE, 
                           random = list(~ 1 | studyNr,~ 1 | obsNr),
                           data = x, method="REML")
                                } else  {
                  mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2, V = V, intercept=TRUE, 
                           random = list(~ 1 | studyNr,~ 1 | obsNr),
                           data = x, method="REML")
              }

              # % change under complete N limitation and slope of N lim effect
              return( c(((exp(cbind(mod.p$b,mod.p$ci.lb, mod.p$ci.ub))-1)*100)[1,], summary(mod.p)$s.nlevels))
            })
res <- as.data.frame(do.call(rbind, res))
colnames(res) <- c("Effect", "lo", "up", "n.studies", "n.total")
res$response <- rownames(res)
res$response <- factor(res$response, levels = rev(res$response))

# plot
fig3 <- ggplot(res, aes(y=Effect, x=response)) +
  geom_errorbar(lwd=0.7, width=0, aes(ymin=lo, ymax=up)) +
  geom_point(shape=21, size=3, fill="white") +
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  xlab("") +
  ylab(expression("Change under complete N limitation ("*"%)")) +
  theme(axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.title.x = element_text(size=14, vjust=5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        axis.ticks.length=unit(-0.20, "cm") 
        
      ) +
  scale_y_continuous(breaks = seq(-80, 120, by=10)) +
                     #,minor_breaks = seq(-75, 5, by=10), limits=c(-80, 120))+
  scale_x_discrete(labels=c("leaf sugar", "leaf starch", "SLA", "Rubisco", 
                            "chlorophyll", expression(N[L]*" per unit mass"),
                            expression(N[L]*" per unit area"), "leaf area", "photosynthesis")) +
  coord_flip()
fig3
ggsave("fig3_ms.png", fig3, width = 25, height = 10, units = "cm", dpi = 300)

# Overview: N lim trends. Figure 4 in ms####

  # list of data frames
    #rep.var <- list(photosynthesis=dat.photo, leaf_area=dat.leafA, leafN_area = dat.leafNA, 
    #           leafN_mass = dat.leafNM, leafChl = dat.chl, leafRubisco = dat.rub)

# Select variables for figure4 and add response names as a new column
fig4.vars <- names(rep.var)[names(rep.var) %in% c("photosynthesis", "leaf_area", "leafN_area", "leafN_mass", "leafChl", 
                                                  "leafRubisco" ,  "leafStarch"
                      )] 
rep.var.fig4 <- rep.var[fig4.vars]
rep.names <- as.list(setNames(fig4.vars, fig4.vars))
rep.var.fig4 <- Map(cbind, rep.var.fig4, response = rep.names)

# run model for all responses and save output of the effect  
res <- lapply(rep.var.fig4, function (x) {
                      x$N_conc_exp2.rev <- 1-x$N_conc_exp2
                      # var-covar matrix
                      V <- v_func(dat = x, ind.study = "ind.study", 
                                    var.control ="var.control", var.rr = "var.rr")
                      if(x$response[1] == "leafN_mass") {
                        x$N_conc_exp2.rev.q <- x$N_conc_exp2.rev*x$N_conc_exp2.rev # qudratic variable
                        mod <- rma.mv(rr ~ N_conc_exp2.rev + N_conc_exp2.rev.q , V = V, intercept=TRUE, 
                                       random = list(~ 1 | studyNr,~ 1 | obsNr),
                                       data = x, method="REML")
                        pred.frame <-  predict(mod, newmods = 
                                                 cbind(seq(min(1-x$N_conc_exp2),1,0.001),
                                                       seq(min(1-x$N_conc_exp2),1,0.001)^2), addx=TRUE)
                        
                                      } else {mod <- rma.mv(yi = rr, mods = ~ N_conc_exp2.rev, V = V, intercept=TRUE, 
                                                            random = list(~ 1 | studyNr,~ 1 | obsNr),
                                                            data = x, method="REML")
                                                               
                                              pred.frame <-  predict(mod, newmods = c(seq(min(1-x$N_conc_exp2),1,0.001)), 
                                                                                      addx=TRUE)}
                      
                      pred.frame <- do.call(cbind.data.frame, pred.frame)
                      pred.frame$pred <- (exp(pred.frame$pred)-1)*100
                      return(pred.frame)
                  })

# plot overall N lim effect
fig4.list <- list(1)
resp.ynames <- list("Change in photosynthesis (%)", "Change in leaf area (%)", expression("Change in "*N[L]*" per unit area (%)"),
                    expression("Change in "*N[L]*" per unit mass (%)"), expression("Change in chlorophyll per unit area (%)"),  
                    expression("Change in Rubisco per unit area (%)"), expression("Change in leaf starch (%)"))
fig.pan <- list("(a)","(b)","(c)","(d)","(e)","(f)", "(g)")
for ( i in 1:length(rep.var.fig4)) {
  # store the data frames for the response
  nlim.plot <- data.frame("Neffect" = (exp(rep.var.fig4[[i]]$rr)-1)*100,
                          "PropNlim" =  (1-rep.var.fig4[[i]]$N_conc_exp2) * 100)
  pred.frame <- res[[i]] # model predictions
  if(i < 7) {
             fig4af <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
                    geom_point(shape=21, size=3, fill="white") +
                    geom_hline(yintercept = 0, lty=2) +
                    geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev*100)) +
                    geom_ribbon(data=pred.frame,aes(x=X.N_conc_exp2.rev*100, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
                                inherit.aes=FALSE, alpha=0.3) +
                    scale_y_continuous( name = resp.ynames[[i]], 
                                      breaks = seq(-100, 100,20), limits = c(-100, 60)) +
                    xlim(c(0,100)) +
                    xlab(expression("N-limitation ("*"%)")) +
                    theme(
                      axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                      axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
                      axis.title.x = element_text(size=14, vjust=5),
                      axis.title.y = element_text(size=14, vjust=-5),
                      axis.line.x = element_line(color="black", size = 1),
                      axis.line.y = element_line(color="black", size = 1),
                      axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.background = element_blank(),
                      axis.ticks.length=unit(-0.20, "cm")  ) +
                    annotation_custom(
                                  grob = textGrob(label = fig.pan[[i]], gp = gpar(fontsize = 23)),
                                  ymin = 55,      # Vertical position of the textGrob
                                  #  ymax = 5,
                                  xmin = -90)
                    fig4.list[[i]] <- ggplot_gtable(ggplot_build(fig4af))
                    fig4.list[[i]]$layout$clip[fig4.list[[i]]$layout$name == "panel"] <- "off"
  } else {
    fig4g <- ggplot(data=nlim.plot, aes(y=Neffect, x=PropNlim))+
      geom_point(shape=21, size=3, fill="white") +
      geom_hline(yintercept = 0, lty=2) +
      geom_line(data=pred.frame, aes(y=pred, x=X.N_conc_exp2.rev*100)) +
      geom_ribbon(data=pred.frame,aes(x=X.N_conc_exp2.rev*100, ymin=(exp(ci.lb)-1)*100,ymax=(exp(ci.ub)-1)*100),
                  inherit.aes=FALSE, alpha=0.3) +
      scale_y_continuous( name = resp.ynames[[i]], 
                          breaks = seq(0,350,50), limits = c(0, 325)) +
      xlim(c(0,100)) +
      xlab(expression("N-limitation ("*"%)")) +
      theme(
        axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.title.x = element_text(size=14, vjust=5),
        axis.title.y = element_text(size=14, vjust=-5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        axis.ticks.length=unit(-0.20, "cm")  ) +
      annotation_custom(
        grob = textGrob(label = fig.pan[[i]], gp = gpar(fontsize = 23)),
        ymin = 320,      # Vertical position of the textGrob
        #  ymax = 5,
        xmin = -90)
    fig4.list[[i]] <- ggplot_gtable(ggplot_build(fig4g))
    fig4.list[[i]]$layout$clip[fig4.list[[i]]$layout$name == "panel"] <- "off"
               }
}

png("figure4a_f.png", width=28, height=45, units="cm", res=300)
grid.arrange(fig4.list[[1]],fig4.list[[2]], fig4.list[[3]], fig4.list[[4]], 
             fig4.list[[5]], fig4.list[[6]], fig4.list[[7]], ncol=2, nrow =4)
dev.off()

# Legume/CO2 modify N lim effect: Figure 5 in ms ####
# rep.var <- list(photosynthesis=dat.photo, leaf_area=dat.leafA, leafN_area = dat.leafNA, 
#                 leafN_mass = dat.leafNM, leafChl = dat.chl, leafRubisco = dat.rub, 
#                 sla = dat.sla, leafStarch = dat.star, leafSugar = dat.sug)
# # add response names as a new column
# rep.names <- as.list(setNames(names(rep.var), names(rep.var)))
# rep.var <- Map(cbind, rep.var, response = rep.names)

# Legume status

# construct model matrix for predictions
pred.leafNM <- rbind(c(0,0,0,0), c(0,0,1,0), c(0,0,0,1)) # qudratic term
pred.other <- rbind(c(0,0,0), c(0,1,0), c(0,0,1))

# run model for all responses and save output of the effect
rep.var.sub <- rep.var[-8] # remove leaf starch as there are no studies with legumes with nodes
res <- lapply(rep.var.sub, function (x) {
  # var-covar matrix
  V <- v_func(dat = x, ind.study = "ind.study", 
              var.control ="var.control", var.rr = "var.rr")
  # model. Leaf N on mass needs a quadratic term
  if(x$response[1] == "leafN_mass") {
    x$quad.Nconc <- x$N_conc_exp2*x$N_conc_exp2
    mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + quad.Nconc + leg, V = V, intercept=TRUE, 
                    random = list(~ 1 | studyNr,~ 1 | obsNr),
                    data = x, method="REML")
    pred.df <- predict(mod.p, newmods = pred.leafNM) # make predictions
    n.lev <- c(sum(mod.p$X[,1]) - sum(mod.p$X[,4:5]), sum(mod.p$X[,4]), sum(mod.p$X[,5])) # sample size per level 
  } else  {
    mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + leg, V = V, intercept=TRUE, 
                    random = list(~ 1 | studyNr,~ 1 | obsNr),
                    data = x, method="REML")
    pred.df <- predict(mod.p, newmods = pred.other) # make predictions
    n.lev <- c(sum(mod.p$X[,1]) - sum(mod.p$X[,3:4]), sum(mod.p$X[,3]), sum(mod.p$X[,4])) # sample size per level 
    }
    pred.frame <- do.call(cbind.data.frame, pred.df)
    pred <- (exp(pred.frame[,c(1,3:4)])-1)*100
    
  # % change under complete N limitation and slope of N lim effect
    # rbind(rep(summary(mod.p)$s.nlevels,3))
  return(cbind(pred, n.lev, "study" = summary(mod.p)$s.nlevels[1], 
         "experiments" = summary(mod.p)$s.nlevels[2], "response" = x$response[1]))
})

res <- as.data.frame(do.call(rbind, res), stringsAsFactors = TRUE)
colnames(res)[1:3] <- c("Effect", "lo", "up")
res$treat <- factor(rep(c("not legum", "legume with nods", "legume, no nods"), 8), 
                        levels =  c( "legume with nods",  "legume, no nods", "not legum"  ))
res$response <- factor(res$response, levels = rev(res$response)) # get same order as legend

# plot
fig5a <- ggplot(res, aes(y=Effect, x=response,group=treat)) +
  geom_errorbar(lwd=0.7, width=0, aes(ymin=lo, ymax=up),position = position_dodge(width = 0.5)) +
  geom_point(aes(shape=treat),size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  xlab("") +
  ylab(expression("Change under complete N-limitation ("*"%)")) +
  scale_shape_manual(name = "Legume status",
                       breaks=c("not legum","legume, no nods","legume with nods"),
                       labels = c("not legume", "legume, no nods", "legume with nods"),
                       values = c(17,15,16)) +
  theme(
        axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.title.x = element_text(size=14, vjust=5),
        axis.title.y = element_text(size=14, vjust=-5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        axis.ticks.length=unit(-0.20, "cm"),
        legend.position=c(.75, .5),
        legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", 
                                         fill = "white", size = 0.7, linetype = "solid")) +
  scale_y_continuous(breaks = seq(-80, 200, by=20), 
                    # minor_breaks = seq(-90, 290, by=10), 
                     limits=c(-100, 200)) +
  scale_x_discrete(labels=c("leaf sugar", "SLA", "Rubisco", 
                            "chlorophyll", expression(N[L]*" per unit mass"),
                            expression(N[L]*" per unit area"), "leaf area", "photosynthesis")) +
  #geom_text(aes(y=200,label=n.lev),hjust=0, vjust=0.35, position = position_dodge(width = 0.5), size=5) +
  geom_text(aes(y=200,label=n.lev),hjust=0.3, vjust=rep(c(-0.2, 0.1, 0.4),8), position = position_dodge(width = 0.5), size=5) +
  geom_text(aes(y=-100, x=8.4, label= "(a)"), size = 12) +
  coord_flip()
fig5a

# CO2 treat

# construct model matrix for predictions
pred.leafNM <- rbind(c(0,0,0), c(0,0,1)) # qudratic term
pred.other <- rbind(c(0,0), c(0,1))

# run model for all responses and save output of the effect
rep.var.sub <- rep.var #[-8]
res <- lapply(rep.var.sub, function (x) {
  # var-covar matrix
  V <- v_func(dat = x, ind.study = "ind.study", 
              var.control ="var.control", var.rr = "var.rr")
  # model. Leaf N on mass needs a quadratic term
  if(x$response[1] == "leafN_mass") {
    x$quad.Nconc <- x$N_conc_exp2*x$N_conc_exp2
    mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + quad.Nconc + CO2, V = V, intercept=TRUE, 
                    random = list(~ 1 | studyNr,~ 1 | obsNr),
                    data = x, method="REML")
    pred.df <- predict(mod.p, newmods = pred.leafNM) # make predictions
    n.lev <- c(sum(mod.p$X[,1]) - sum(mod.p$X[,4]), sum(mod.p$X[,4])) # sample size per level 
  } else  {
    mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + CO2, V = V, intercept=TRUE, 
                    random = list(~ 1 | studyNr,~ 1 | obsNr),
                    data = x, method="REML")
    pred.df <- predict(mod.p, newmods = pred.other) # make predictions
    n.lev <- c(sum(mod.p$X[,1]) - sum(mod.p$X[,3]), sum(mod.p$X[,3])) # sample size per level 
  }
  pred.frame <- do.call(cbind.data.frame, pred.df)
  pred <- (exp(pred.frame[,c(1,3:4)])-1)*100
  
  # % change under complete N limitation and slope of N lim effect
  # rbind(rep(summary(mod.p)$s.nlevels,3))
  return(cbind(pred, n.lev, "study" = summary(mod.p)$s.nlevels[1], 
               "experiments" = summary(mod.p)$s.nlevels[2], "response" = x$response[1]))
})

res <- as.data.frame(do.call(rbind, res), stringsAsFactors = TRUE)
colnames(res)[1:3] <- c("Effect", "lo", "up")
res$treat <- factor(rep(c("ambient", "elevated"), 9), 
                    levels =  c("elevated", "ambient"))
res$response <- factor(res$response, levels = rev(res$response)) # get same order as legend

# plot
fig5b <- ggplot(res, aes(y=Effect, x=response,group=treat)) +
  geom_errorbar(lwd=0.7, width=0, aes(ymin=lo, ymax=up),position = position_dodge(width = 0.5)) +
  geom_point(aes(shape=treat),size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  xlab("") +
  ylab(expression("Change under complete N-limitation ("*"%)")) + 
  scale_shape_manual(name = "CO2", 
                     breaks = c("ambient", "elevated"),
                     labels = c(expression(CO[2]*" ambient "), expression(CO[2]*" elevated")),
                     values = c(15, 16)) +
  theme(
        axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.title.x = element_text(size=14, vjust=5),
        axis.title.y = element_text(size=14, vjust=-5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        axis.ticks.length=unit(-0.20, "cm"),
        legend.position=c(.75, .5),
        legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", 
                                          fill = "white", size = 0.7, linetype = "solid")) +
  scale_y_continuous(breaks = seq(-80, 200, by=20), 
                     # minor_breaks = seq(-90, 290, by=10), 
                     limits=c(-100, 200)) +
  scale_x_discrete(labels=c("leaf sugar", "leaf starch", "SLA", "Rubisco", 
                            "chlorophyll", expression(N[L]*" per unit mass"),
                            expression(N[L]*" per unit area"), "leaf area", "photosynthesis")) +
  geom_text(aes(y=200,label=n.lev),hjust=0.3, vjust=rep(c(-0.1, 0.4),9), position = position_dodge(width = 0.5), size=5) +
  geom_text(aes(y=-100, x=9.3, label= "(b)"), size = 12) +
  coord_flip()
fig5b

png("figure5ab.png", width=26, height=32, units="cm", res=300)
grid.arrange(fig5a, fig5b, ncol=1, nrow =2)
dev.off()


# Fig 1S N source####
# N source

# construct model matrix for predictions
pred.leafNM <- rbind(c(0,0,0,0,0), c(0,0,1,0,0), c(0,0,0,1,0), c(0,0,0,0,1)) # qudratic term
pred.other <- rbind(c(0,0,0), c(0,1,0), c(0,0,1))

# run model for all responses and save output of the effect
rep.var.sub <- rep.var[-c(8:9)] # remove leaf starch as there are no studies with legumes with nodes
res <- lapply(rep.var.sub, function (x) {
  # var-covar matrix
  V <- v_func(dat = x, ind.study = "ind.study", 
              var.control ="var.control", var.rr = "var.rr")
  # model. Leaf N on mass needs a quadratic term
  if(x$response[1] == "leafN_mass") {
    x$quad.Nconc <- x$N_conc_exp2*x$N_conc_exp2
    mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + quad.Nconc + Nsource, V = V, intercept=TRUE, 
                    random = list(~ 1 | studyNr,~ 1 | obsNr),
                    data = x, method="REML")

    pred.df <- predict(mod.p, newmods = pred.leafNM) # make predictions

    n.lev <- c(sum(mod.p$X[,1]) - sum(mod.p$X[,4:6]), sum(mod.p$X[,4]), sum(mod.p$X[,5]), sum(mod.p$X[,6])) # sample size per level 
  } else  {
    mod.p <- rma.mv(yi = rr, mods = ~ N_conc_exp2 + Nsource, V = V, intercept=TRUE, 
                    random = list(~ 1 | studyNr,~ 1 | obsNr),
                    data = x, method="REML")
    pred.df <- predict(mod.p, newmods = pred.other) # make predictions
    n.lev <- c(sum(mod.p$X[,1]) - sum(mod.p$X[,3:4]), sum(mod.p$X[,3]), sum(mod.p$X[,4])) # sample size per level 
  }
  pred.frame <- do.call(cbind.data.frame, pred.df)
  pred <- (exp(pred.frame[,c(1,3:4)])-1)*100

    # % change under complete N limitation and slope of N lim effect
  # rbind(rep(summary(mod.p)$s.nlevels,3))
  return(cbind(pred, n.lev, "study" = summary(mod.p)$s.nlevels[1], 
               "experiments" = summary(mod.p)$s.nlevels[2], "response" = x$response[1]))
})

res <- as.data.frame(do.call(rbind, res), stringsAsFactors = TRUE)
colnames(res)[1:3] <- c("Effect", "lo", "up")
res$treat <- factor(c(rep(c("NH4+", "NO3-", "NH4+ - NO3-"), 3), c("NH4+", "NO3-", "NH4+ - NO3-", "urea"),
                      rep(c("NH4+", "NO3-", "NH4+ - NO3-"), 3)),
                    levels =  c("NH4+", "NO3-", "NH4+ - NO3-", "urea"))
res$response <- factor(res$response, levels = rev(res$response)) # get same order as legend

# plot
res$treat <- relevel(res$treat, ref = "urea") 
figS1 <- ggplot(res, aes(y=Effect, x=response,group=treat)) +
  geom_errorbar(lwd=0.7, width=0, aes(ymin=lo, ymax=up),position = position_dodge(width = 0.5)) +
  geom_point(aes(shape=treat),size = 3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  xlab("") +
  ylab(expression("Change under complete N-limitation ("*"%)")) +
  scale_shape_manual(name = "N source",
                     breaks=c(  "NH4+ - NO3-", "NO3-","NH4+","urea"),
                     #labels = c("NH4+", "NO3-", "NH4+ - NO3-", "urea"),
                     labels = c("NH4+ - NO3-", "NO3-", "NH4+", "urea"),
                     values = c(18,16,15,17)) +
  theme(
    axis.text.x  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y  = element_text(size=14, color="black", margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.title.x = element_text(size=14, vjust=5),
    axis.title.y = element_text(size=14, vjust=-5),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    axis.ticks.length=unit(-0.20, "cm"),
    legend.position=c(.75, .5),
    legend.text = element_text(size=18),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.background = element_rect(color = "black", 
                                     fill = "white", size = 0.7, linetype = "solid")) +
  scale_y_continuous(breaks = seq(-90, 50, by=20), 
                     # minor_breaks = seq(-90, 290, by=10), 
                     limits=c(-100, 50)) +
  scale_x_discrete(labels=c("SLA", "Rubisco", 
                            "chlorophyll", expression(N[L]*" per unit mass"),
                            expression(N[L]*" per unit area"), "leaf area", "photosynthesis")) +
  geom_text(aes(y=50,label=n.lev),hjust=0, vjust=0, 
            position = position_dodge(width = 0.75), size=5) +
  coord_flip()
figS1

png("figureS1.png", width=26, height=16, units="cm", res=300)
grid.arrange(figS1, ncol=1, nrow =1)
dev.off()



# Test analyses ####
## Test how to impute between response dependence within studies
require(metafor)

# make data
set.seed(1)
dat <- data.frame(exp = gl(8,4), resp = gl(2,2, length=32), yi = rnorm(32,10,4),
      vi = rnorm(32,4,0.5), vi.cov = rnorm(32,0.5,0.25) )
dat$vi.cov <- rep(c(0.1,0.15),8, each=2) 
dat$st.depend <- interaction(dat$exp, dat$resp)

# known var-cov matrix (within-study dependence for each response)  
V = v_func(dat = dat, study.dep = "st.depend", var.y = "vi" , var.y.cov = "vi.cov")

# run model
res.no.covar <- rma.mv(yi, V, mods = ~ resp - 1, 
              random = ~ resp | exp, 
              struct="UN", data=dat, method="REML")
print(res.no.covar, digits=3)
# correlation = -0.613

# test with correlation
cor.r = 0.5

covaris = list()
for (i in 1:nlevels(dat$exp)) {
  studid <- dat$exp == levels(dat$exp)[i]
  vari <- dat[ studid, "vi"]
  
  # make cor matrix
  a <- matrix(nrow=length(vari), ncol=length(vari))
  diag(a) <- 1
  a[lower.tri(a)] <-  cor.r
  a[upper.tri(a)] <- cor.r
  
  # calculate covariance for the given cor coef 
  b <- sqrt(vari) %*% t(sqrt(vari))
  a_cov <- b * a  # covariance matrix
  # cov2cor(a_cov) #test
  
  # store var-cov block for each study
  covaris[[i]] <- a_cov 
}

# var-cov matrix for the between response correlation
# within each study
V2 <- as.matrix(bdiag(lapply(covaris, function (x) x) ))
# cov2cor(V2) #test if cor is correct

# add correlation between responses to 
# within-response matrix
V = ifelse(V==0, V2, V)

# run mode with complete var-covar matrix
res.with.covar <- rma.mv(yi, V, mods = ~ resp - 1, 
                       random = ~ resp | exp, 
                       struct="UN", data=dat, method="REML")
print(res.no.covar, digits=3)
# correlation still -0.613....



v_func <- function (dat, study.dep, var.y , var.y.cov) {
 
  dat.agg <- data.frame(st = unique(dat[[study.dep]]), len = rle(as.numeric(dat[[study.dep]]))$lengths, 
                        var = unique(dat[[var.y.cov]]))
  require(Matrix)  
  ll = list()
  for (i in 1:NROW(dat.agg)) {  
    ll[[i]] <- matrix(dat.agg$var[i], ncol = dat.agg$len[i], nrow = dat.agg$len[i]) 
  }
  mat <- as.matrix(bdiag(lapply(ll, function (x) x) ))
  diag(mat) <- dat[[var.y]]
  return(mat)
}
