################################################################################
# Nitrogen and crop photosynthesis: A meta-analysis of crop response patterns 
# to nitrogen limitation for improved model representation
#
# R code to reproduce results in Seufert V, Granath G, Muller C.
#
# Contact (R-code): gustaf.granath@gmail.com
################################################################################

# Load packages

# Statistical analyses
library(metafor) # tested with version 2.0-0

# Plotting and data management
library(ggplot2)
library(gridExtra)
library(grid)
library(readxl)
library(reshape2)

# Code overview
# First run section 'Load functions and data' to import data and load functions. This 
# section ends with 'END LOAD DATA AND FUNCTIONS'
# After this, all other sections can be run.

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
# END LOAD DATA AND FUNCTIONS


# Test differences between area and mass for Chl, Rub, sugar, leaf area ####
# Mostly area based but a few mass based experiments.
# For leaf area, mostly at the plant level but a few experiments with only leaf level.
# Below we check if this can be ignored in our analyses and that we can include
# all experiments to increase sample size.

# Chlorophyll
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

# sugar
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
# Complete test of all predictors defined in the object test.var
# Produces the file large_table.csv

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
# PLots to get an overview of the response ratio (RR) for all responses
# in relation to all predictors
# Produces png files (one per response) showing the results

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


# Compare results with full Bayes: MCMCglmm package ####
# Code to compare a the results with a Bayesian approach using the 
# MCMCglmm package
library(MCMCglmm)

# Here an example test with photosynthesis
# but other responses can be tested as well
dat.photo <- rep.var[["photosynthesis"]]

# setup fixed var-covar matrix for MCMCglmm
V <- v_func(dat = dat.photo) # create var-covar matrix
Rsvd <- svd(V)
Rsvd <- Rsvd$v%*%(t(Rsvd$u)*sqrt(Rsvd$d))
dat.photo$exp <- as.factor(1:NROW(dat.photo))
Z <- model.matrix(~exp-1, dat.photo)%*%Rsvd
dat.photo$Z <- Z
# define priors
prior = list(R = list(V = diag(1)*1e-6, nu = 2), 
             G = list(G1 = list(V = 1, fix = 1),
                      G2 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 25^2)))

mod.leg <- MCMCglmm(rr ~ N_conc_exp2 + leg,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, prior=prior,
               nitt=100000,burnin=15000, # adjust for sufficient sampling 
               pr=TRUE)
summary(mod.leg)
plot(mod.leg$VCV)

# check study effect
colnames(mod.leg$Sol)
ref <- data.frame(mod.leg$Sol[,187:236])
#colnames(ref) <- letters[1:4]
rdf <- melt(ref)
ggplot(data=rdf,aes(x=value, color=variable))+geom_density()

# run simpler model
mod.Nlim <- MCMCglmm(rr ~ N_conc_exp2,  
               random=~idv(Z) + studyNr, 
               data=dat.photo, prior=prior,
               nitt=100000,burnin=15000, 
               pr=TRUE)
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

### END MCMCglmm testing

# Figure 3 in ms - Overview: full N lim. ####

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

# Figure 4 in ms - Overview: N lim trends.####

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

# Figure 5 in ms - Legume/CO2 modify N lim effect####
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