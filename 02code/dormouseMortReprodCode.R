# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2020-01-28
# LAST MODIFIED: 2024-08-15
# DESCRIPTION: Reproducibility code for Berg et al. (2024). 
# COMMENTS: Temporal Bayesian Survival Trajectory Analyses (BaSTA) 
#           on the dormouse data.
# ================================ CODE START ================================ #
# ================ #
# ==== SETUP: ==== 
# ================ #
# Character vector of packages installed (used to determine whether to install
# the required packages)
instPacks <- installed.packages()[, 1]

# Install devtools to load packages from GitHub:
if (!"devtools" %in% instPacks) {
  install.packages("devtools")
  library(devtools)
}

# BaSTA for survival analysis:
if (!"BaSTA" %in% instPacks) {
  install_git("https://github.com/fercol/basta2.0", subdir = "pkg/")
}

# paramDemo for demographic functions used after BaSTA:
if (!"paramDemo" %in% instPacks) {
  install_git("https://github.com/fercol/paramDemo", subdir = "pkg/")
}

# snowfall for parallel computing:
if (!"snowfall" %in% instPacks) {
  install.packages("snowfall")
}

# RColorBrewer for plotting colors:
if (!"RColorBrewer" %in% instPacks) {
  install.packages("RColorBrewer")
}

# Load libraries:
library(BaSTA)
library(snowfall)
library(RColorBrewer)
library(paramDemo)

# Set working directory (change accordingly):
setwd("Path to.../DormouseMortReproducibility/")
setwd("~/FERNANDO/PROJECTS/4.PACKAGES/BergEtal2025DorMouseReprodCode/")

# Logical to save plots:
saveResults <- FALSE

# ===================== #
# ==== DATA PREP.: ====
# ===================== #
# Load BaSTA capture recapture matrix:
bastaDat <- read.csv("03data/tables/DormouseBaSTAdata.csv", header = TRUE,
                     stringsAsFactors = FALSE)

# Sexes:
sexes <- c(f = "Female", m = "Male")

# Sex label:
sxlab <- c("f", "m")

# Years of start and end:
ystart <- c(1999, 2007, 2015)
yend <- c(2006, 2014, 2022)

# Number of intervals:
nyints <- length(ystart)

# Interval labels:
intslab <- sapply(1:nyints, function(yy) {
  return(sprintf("%s-%s", ystart[yy], yend[yy]))
}) 

# Year vector:
Ycols <- grep("X", colnames(bastaDat))

# ==================== #
# ==== FUNCTIONS: ====
# ==================== #
# Kullback-Leibler discrepancies:
CalcKL <- function(x1, x2, length = 1000) {
  # 99 Quantiles:
  qx1 <- quantile(x1, c(0.001, 0.999))
  qx2 <- quantile(x2, c(0.001, 0.999))
  
  # Find range:
  xRange <- range(c(qx1, qx2))
  
  # Sequence of values:
  exseq <- seq(xRange[1], xRange[2], length = length)
  dx <- diff(exseq)[1]
  nseq <- length(exseq)
  
  # Densities:
  dx1 <- density(x1, n = nseq, from = xRange[1], to = xRange[2])
  dx2 <- density(x2, n = nseq, from = xRange[1], to = xRange[2])
  
  # Standardize densities:
  dx1$y <- dx1$y / sum(dx1$y * dx)
  dx2$y <- dx2$y / sum(dx2$y * dx)
  
  # Find zeros:
  ide <- which(dx1$y > 0 & dx2$y > 0)
  
  # KLd:
  kld1 <- sum(dx1$y[ide] * log(dx1$y[ide] / dx2$y[ide]) * dx)
  kld2 <- sum(dx2$y[ide] * log(dx2$y[ide] / dx1$y[ide]) * dx)
  
  # Standardized KL:
  qKlc1 <- (1 + (1 - exp(-2 * kld1)^(1 / 2))) / 2
  qKlc2 <- (1 + (1 - exp(-2 * kld2)^(1 / 2))) / 2
  
  # Mean standardized KL:
  mqKl <- (qKlc1 + qKlc2) / 2
  
  # Output:
  outList <- c(klfm = kld1, klmf = kld2, qklfm = qKlc1, 
               qklmf = qKlc2, mqKl = mqKl)
  return(outList)
}

# Truncated normal pdf:
dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

# Truncated normal quantile function:
qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# =========================================== #
# ==== BaSTA ANALYSIS BY YEAR INTERVALS: ====
# =========================================== #
# Adult models to run:
models <- c("GO", "WE", "LO")
nmods <- length(models)

# Output for model selection matrix:
DICmat <- list(Female = data.frame(Period = intslab, GO = rep(NA, nyints), 
                                   WE = rep(NA, nyints), LO = rep(NA, nyints)),
               Male = data.frame(Period = intslab, GO = rep(NA, nyints), 
                                 WE = rep(NA, nyints), LO = rep(NA, nyints)))

# List outputs:
coefList <- list()
dataList <- list()
outList <- list()

# Matrix of number of individuals per sex and period:
Nmat <- data.frame(Period = intslab, Female = rep(NA, nyints),
                   Male = rep(NA, nyints))

# Run analysis:
for (yy in 1:nyints) {
  # Vector of years within period:
  yv <- ystart[yy]:yend[yy]
  
  # Number of years:
  ny <- length(yv)
  
  # sublists for results:
  outList[[intslab[yy]]] <- list()
  coefList[[intslab[yy]]] <- list()
  
  # Model selection list:
  modList <- list()
  
  # Run analyses per sex:
  for (sx in sxlab) {
    # Index of rows to subset data within period:
    idkeep <- which((bastaDat$Birth < yend[yy] | 
                       bastaDat$Max.Birth < yend[yy]) & 
                      (bastaDat$Death > ystart[yy] | 
                         bastaDat$Min.Death > ystart[yy]) & bastaDat$Sex == sx)
    
    # Index of columns to include:
    idcols <- c(colnames(bastaDat)[1:3], sprintf("X%s", yv), 
                colnames(bastaDat)[c(Ycols[length(Ycols)] + 1):ncol(bastaDat)])
    
    # Subset data:
    ybastadat <- droplevels(bastaDat[idkeep, idcols])
    
    # Run data check to make sure there are no remaining issues:
    ydcheck <- DataCheck(ybastadat, studyStart = ystart[yy], 
                         studyEnd = yend[yy])
    
    # Fix data based on data check (if necessary):
    yfixdat <- FixCMRdata(ybastadat, studyStart = ystart[yy], 
                          studyEnd = yend[yy],
                          autofix = rep(1, 6), silent = TRUE)
    
    # Final analysis data:
    ydat <- yfixdat$newData
    
    # Number of records used for BaSTA:
    Nmat[[sexes[sx]]][yy] <- nrow(ydat)
    
    # Make sure that the year columns are numeric:
    for (bd in c("Min.Birth", "Max.Birth", "Min.Death", "Max.Death")) {
      ydat[[bd]] <- as.numeric(ydat[[bd]])
    }
    
    # Empty DIC vector for sex and period:
    DICs <- rep(NA, nmods)
    names(DICs) <- models
    
    # Prepare list of coefficients:
    coefList[[intslab[yy]]][[sexes[sx]]] <- list()
    
    # Run model selection:
    for (imod in 1:nmods) {
      outi <- basta(object = ydat, studyStart = ystart[yy], 
                    studyEnd = yend[yy], model = models[imod],
                    shape = "bathtub", nsim = 8, 
                    parallel = TRUE, ncpus = 8,
                    niter = 45000, burnin = 5001, thinning = 50)
      
      # Fill up coefficient list:
      coefList[[intslab[yy]]][[sexes[sx]]][[models[imod]]] <- outi$coeff
      
      cat("\n=======================================\n")
      cat(sprintf("\nPeriod %s\nSex %s\nModel %s\n", yy, sexes[sx], 
                  models[imod]))
      summary(outi)
      cat("\n=======================================\n")
      plot(outi, plot.type = "gof")

      # Fill up model list:
      modList[[sexes[sx]]][[models[imod]]] <- outi
      
      # Fill up DIC matrix:
      if (!is.na(outi$DIC[1])) DICs[models[imod]] <- outi$DIC["DIC"]
      DICmat[[sexes[sx]]][[models[imod]]][yy] <- DICs[models[imod]]
    }
    
    # Find model with lowest DIC:
    idmod <- which(DICs == min(DICs, na.rm = TRUE))
    
    # Save selected model in output list:
    out <- modList[[sexes[sx]]][[models[idmod]]]
    outList[[intslab[yy]]][[sexes[sx]]] <- out
  }
}

# Save final results:
if (saveResults) {
  save(list = c("outList", "coefList", "Nmat", "DICmat", "sexes", "sxlab", 
                "nyints", "intslab"), 
       file = "04results/rdata/DormouseMortModelSelection.RData")
}

# ========================== #
# ==== EXTRACT RESULTS: ====
# ========================== #
# ----------------------------- #
# ---- Summary data table: ----
# ----------------------------- #
# save the DIC matrix:
DICfin <- data.frame(Sex = rep(sexes, each = 3), N = c(Nmat$Female, Nmat$Male),
                     rbind(DICmat$Female, DICmat$Male))

if (saveResults) {
  write.csv(DICfin, file = "04results/tables/Table1_DICs.csv", 
            row.names = FALSE)
}

# -------------------------- #
# ---- Life expectancy: ----
# -------------------------- #
npost <- nrow(outList$`1999-2006`$Female$params)
LEpost <- matrix(NA, npost, nyints * 2, 
                 dimnames = list(NULL, 
                                 paste("LE", rep(sexes, each = nyints),
                                       rep(sprintf("Per%s", 1:nyints), 2), 
                                       sep = ".")))

for (sx in 1:2) {
  sxName <- sexes[sx]
  for (iy in 1:nyints) {
    # parameter names:
    cnames <- colnames(outList[[iy]][[sxName]]$params)
    
    # Extract theta matrix:
    theMat <- outList[[iy]][[sxName]]$params[, which(cnames != "pi")]
    
    # Extract model:
    model <- outList[[iy]][[sxName]]$modelSpecs["model"]
    
    # Calculate life expectancy posterior:
    exv <- apply(theMat, 1, function(thei) {
      CalcRemainLifeExp(theta = thei, x = 0, dx = 0.001, xmax = 20, 
                        model = model, shape = "bathtub")["RemLExp"]
    })
    
    # life expectancy summary statistics:
    exstats <- c(Mean = mean(exv), SD = sd(exv), 
                 Lower = c(quantile(exv, 0.025, names = FALSE)),
                 Upper = c(quantile(exv, 0.975, names = FALSE)))
    
    # Store summary statistics:
    temp <- data.frame(Sex = sxName, Period = intslab[iy], t(exstats))
    if (sx == 1 & iy == 1) {
      lifeExp <- temp
    } else {
      lifeExp <- rbind(lifeExp, temp)
    }
    
    # Store life expectancy posteriors:
    LEpost[, sprintf("LE.%s.Per%s", sxName, iy)] <- exv
  }
}

# Calculate KL:
klList <- list()
klMat <- data.frame(Periods = paste(intslab[c(1, 1, 2)], 
                                    intslab[c(2, 3, 3)], sep = "/"),
                    LifeExp = rep(NA, 3))
for (sx in sexes) {
  klList[[sx]] <- list()
  klMat0 <- klMat
  ii <- 0
  for (yi in 1:(nyints - 1)) {
    for (yj in (yi + 1):nyints) {
      exvij <- LEpost[, c(sprintf("LE.%s.Per%s", sx, yi), 
                          sprintf("LE.%s.Per%s", sx, yj))]
      ii <- ii + 1
      kli <- CalcKL(x1 = exvij[, 1], x2 = exvij[, 2], length = 10000)
      klMat0$LifeExp[ii] <- (kli["mqKl"] - 0.5) * 2
    }
  }
  klList[[sx]] <- klMat0
}

KLperiodLE <- rbind(data.frame(sex = rep("Female", 3), klList$Female),
                    data.frame(sex = rep("Male", 3), klList$Male))


# KL between sexes within periods:
KLsexLE <- data.frame(Periods = intslab, lifeExp = rep(NA, 3))
for (iy in 1:nyints) {
  exvij <- LEpost[, c(sprintf("LE.Female.Per%s", iy), 
                      sprintf("LE.Male.Per%s", iy))]
  ii <- ii + 1
  kli <- CalcKL(x1 = exvij[, 1], x2 = exvij[, 2], length = 10000)
  KLsexLE$lifeExp[iy] <- (kli["mqKl"] - 0.5) * 2
}

# ---------------------- #
# ---- Aging rates: ----
# ---------------------- #
ARages <- c(2, 4, 6)
nARage <- length(ARages)
ARpostNames <- paste(rep(sprintf("AR.Age%s", ARages), each = 2 * nyints), 
                     rep(paste(rep(sexes, each = nyints), 
                               rep(sprintf("Per%s", 1:nyints), 2), sep = "."), nARage),
                     sep = ".")
nARpost <- length(ARpostNames)
ARpost <- matrix(NA, npost, nARpost, dimnames = list(NULL, ARpostNames))

for (sx in 1:2) {
  sxName <- sexes[sx]
  for (iy in 1:nyints) {
    # Parameter names:
    cnames <- colnames(outList[[iy]][[sxName]]$params)
    
    # Theta matrix:
    theMat <- outList[[iy]][[sxName]]$params[, which(cnames != "pi")]
    
    # Model:
    model <- outList[[iy]][[sxName]]$modelSpecs["model"]
    
    # Aging rates at AR ages:
    arv <- apply(theMat, 1, function(thei) {
      CalcAgeingRateMort(theta = thei, x = ARages, 
                         model = model, shape = "bathtub")[, "AR"]
    })
    
    # Summary statistics:
    arstats <- cbind(Mean = apply(arv, 1, mean), SD = apply(arv, 1, sd), 
                     Lower = apply(arv, 1, quantile, probs = 0.025, 
                                   names = FALSE),
                     Upper = apply(arv, 1, quantile, probs = 0.975, 
                                   names = FALSE))
    
    # Store summary statistics:
    temp <- data.frame(Sex = rep(sxName, nARage), 
                       Period = rep(intslab[iy], nARage), 
                       Age = ARages, arstats)
    if (sx == 1 & iy == 1) {
      ageingRate <- temp
    } else {
      ageingRate <- rbind(ageingRate, temp)
    }
    
    # Store posterior values:
    postName <- sprintf("AR.Age%s.%s.Per%s", ARages, sxName, iy)
    ARpost[, postName] <- t(arv)
  }
}

# KL between periods:
klList <- list()
klMat <- data.frame(Periods = paste(intslab[c(1, 1, 2)], 
                                    intslab[c(2, 3, 3)], sep = "/"),
                    matrix(0, 3, nARage, 
                           dimnames = list(NULL, sprintf("Age%s", ARages))))
for (sx in sexes) {
  klList[[sx]] <- list()
  klMat0 <- klMat
  for (ia in 1:nARage) {
    ii <- 0
    for (yi in 1:(nyints - 1)) {
      for (yj in (yi + 1):nyints) {
        ii <- ii + 1
        iname <- sprintf("AR.Age%s.%s.Per%s", ARages[ia], sx, yi)
        jname <- sprintf("AR.Age%s.%s.Per%s", ARages[ia], sx, yj)
        kli <- CalcKL(x1 = ARpost[, iname], x2 = ARpost[, jname], 
                      length = 10000)
        klMat0[[sprintf("Age%s", ARages[ia])]][ii] <- (kli["mqKl"] - 0.5) * 2
      }
    }
    klList[[sx]] <- klMat0
  }  
}

# KL table to be stored:
KLperiodAR <- data.frame(Sex = rep(sexes, each = 3), 
                         rbind(klList$Female, klList$Male))


# ------------------------------------------------------- #
# ---- Juvenile parameters and age indep. parameter: ----
# ------------------------------------------------------- #
# Parameter names:
parNames <- c("a0", "a1", "c")
nPars <- length(parNames)

# Extract juvenile parameters.
sexPer <- paste(rep(sexes, each = nyints), rep(sprintf("Per%s", 1:nyints), 2),
                sep = ".")
paramsPostNames <- paste(rep(parNames, each = nyints * 2), 
                         rep(sexPer, nPars), sep = ".")
nParPost <- length(paramsPostNames)

# Matrix of posterior values:
paramPost <- matrix(NA, npost, nParPost, 
                    dimnames = list(NULL, paramsPostNames))

for (ip in 1:nPars) {
  for (sx in sexes) {
    for (iy in 1:nyints) {
      ipname <- sprintf("%s.%s.Per%s", parNames[ip], sx, iy)
      paramPost[, ipname] <- outList[[iy]][[sx]]$params[, parNames[ip]]
    }
  }
}


# Calculate KL:
klList <- list()
klMat <- data.frame(Periods = paste(intslab[c(1, 1, 2)], 
                                    intslab[c(2, 3, 3)], sep = "/"),
                    matrix(0, 3, nPars, dimnames = list(NULL, parNames)))
for (sx in sexes) {
  klList[[sx]] <- list()
  klMat0 <- klMat
  for (ip in 1:nPars) {
    parSx <- sprintf("%s.%s", parNames[ip], sx)
    ii <- 0
    for (yi in 1:(nyints - 1)) {
      for (yj in (yi + 1):nyints) {
        ii <- ii + 1
        ipname <- sprintf("%s.%s.Per%s", parNames[ip], sx, yi)
        jpname <- sprintf("%s.%s.Per%s", parNames[ip], sx, yj)
        kli <- CalcKL(x1 = paramPost[, ipname], x2 = paramPost[, jpname], 
                      length = 10000)
        klMat0[[parNames[ip]]][ii] <- (kli["mqKl"] - 0.5) * 2
      }
    }
    klList[[sx]] <- klMat0
  }  
}

KLperiodPARS <- rbind(data.frame(sex = rep("Female", 3), klList$Female),
                      data.frame(sex = rep("Male", 3), klList$Male))


# -------------------------------------------------- #
# ---- Estimated parameters per sex and period: ----
# -------------------------------------------------- #
for (sx in sexes) {
  for (iy in 1:nyints) {
    # Extract model:
    model <- outList[[iy]][[sx]]$modelSpecs["model"]
    
    # Extract coefficients:
    coefsi <- outList[[iy]][[sx]]$coefficients[, c(1:4, 7)]
    
    # number of coefficients:
    ncoefs <- nrow(coefsi)
    
    # Create data frame:
    temp <- data.frame(Sex = rep(sx, ncoefs), 
                       Period = rep(intslab[iy], ncoefs), 
                       Model = rep(model, ncoefs),
                       Parameter = rownames(coefsi), coefsi)
    if (sx == sexes[1] & iy == 1) {
      paramMat <- temp
    } else {
      paramMat <- rbind(paramMat, temp)
    }
  }
}

rownames(paramMat) <- NULL



# ----------------------- #
# ---- Save results: ----
# ----------------------- #
if (saveResults) {
  save(list = c("LEpost", "lifeExp", "ARages", "ARpost", "ageingRate", 
                "KLperiodLE", "KLsexLE", "KLperiodAR", "paramPost", 
                "KLperiodPARS", "paramMat"), 
       file = "04results/rdata/DormouseMortResultsSummaries.RData")
}

# =================================== #
# ==== RESULTS TABLES FOR PAPER: ====
# =================================== #
# Period values of Life expectancy and aging rates:
vars <- c("LifeExp", sprintf("AgingRateAge%s", ARages))
nvars <- length(vars)

for (sx in sexes) {
  for (iv in 1:nvars) {
    if (iv == 1) {
      temp <- lifeExp[which(lifeExp$Sex == sx), ]
    } else {
      temp <- ageingRate[which(ageingRate$Sex == sx & 
                                 ageingRate$Age == ARages[iv-1]), -3]
    }
    temp2 <- data.frame(Sex = temp[, 1], Variable = rep(vars[iv], 3),
                        temp[, -c(1)])
    if (sx == sexes[1] & iv == 1) {
      periodTab <- temp2
    } else {
      periodTab <- rbind(periodTab, temp2)
    }
  }
}

periodTab[, c("Mean", "SD", "Lower", "Upper")] <- 
  signif(periodTab[, c("Mean", "SD", "Lower", "Upper")], 3)

if (saveResults) {
  write.csv(periodTab, file = "04results/tables/Table2_varsByPeriod.csv",
            row.names = FALSE)
}

# KL values of Life expectancy and aging rates:
KLperiod <- cbind(KLperiodPARS, lifeExp = KLperiodLE[, -c(1:2)], KLperiodAR[, -c(1:2)])
KLperiod[, -c(1, 2)] <- signif(KLperiod[, -c(1, 2)], 2)

if (saveResults) {
  write.csv(KLperiod, file = "04results/tables/Table3_KLperiod.csv",
            row.names = FALSE)
  
}

# Parameter matrix:
paramMat[, -c(1:4)] <- signif(paramMat[, -c(1:4)], 3)
if (saveResults) {
  write.csv(paramMat, file = "04results/tables/TableS2S3_paramMat.csv",
            row.names = FALSE)
}

# ======================= #
# ==== PLOT RESULTS: ====
# ======================= #
# ----------------------------- #
# Mortality and survival plots:
# ----------------------------- #
lty <- 1:nyints
cols <- brewer.pal(4, "Set1")
xlim <- c(0, 5)
ylim <- list(surv = c(0, 1), mort = c(0, 3))
xv <- seq(0, 10, 0.001)
dy <- c(surv = 0.2, mort = 0.5)
sexes <- c(f = "Female", m = "Male")
dens <- seq(20, 30, length = nyints)
angle <- seq(45, 135, length = nyints)
yaxname <- c(surv = "Survival", mort = "Hazard rate")

laymat <- cbind(c(0, 2, 3, 0), c(4, 5, 6, 1), c(7, 8, 9, 1))
widths <- c(0.2, 1, 1)
heights <- c(0.1, 0.75, 0.75, 0.25)
mar <- c(2, 1, 1, 1)
plw <- 6.5
plh <- plw / (sum(widths) / sum(heights))
cex.lab <- 2
cex.ax <- 1.5
if (saveResults) {
  pdf(file = "04results/plots/paperPlots/Fig2.pdf", width = plw,
      height = plh)
}

layout(laymat, widths = widths, heights = heights)

# X-axis label:
par(mar = mar * c(0, 1, 0, 1))
plot(xlim, c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
text(mean(xlim), 0.95, "Age (years)", cex = cex.lab, xpd = NA)
xl <- xlim[1] + diff(xlim) * c(0.05, 0.35, 0.65)
xoff <- 0.1
yl <- 0.2
text(x = mean(xlim), yl + 0.3, "Periods", cex = 1.5, font = 2)
for (il in 1:3) {
  lines(xl[il] + diff(xlim) * c(0, xoff), rep(yl, 2), col = cols[il],
        lwd = 4, lty = lty[il])
  text(x = xl[il] + diff(xlim) * (xoff + 0.015), y = yl, intslab[il], adj = 0, 
       cex = 1.25)
}

# Y-axis label:
for (idem in c("surv", "mort")) {
  par(mar = mar * c(1, 0, 1, 0))
  plot(c(0, 1), ylim[[idem]], col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.25, mean(ylim[[idem]]), yaxname[idem], srt = 90, cex = cex.lab)
  axis(side = 4, at = seq(ylim[[idem]][1], ylim[[idem]][2], dy[idem]), 
       pos = 0.6, las = 2, lwd = NA, cex.axis = 1.2)
  
}

ilet <- 0
for (sx in c("f", "m")) {
  # Sex labels:
  par(mar = mar * c(0, 1, 0, 1))
  plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.5, 0.5, sexes[sx], cex = cex.lab)
  
  sxname <- sprintf("Sex%s", sx)
  for (idem in c("surv", "mort")) {
    ilet <- ilet + 1
    par(mar = mar)
    plot(xlim, ylim[[idem]], col = NA, xlab = "", ylab = "", axes = FALSE)
    text(xlim[1] + diff(xlim) * 0.05, ylim[[idem]][2] - 
           diff(ylim[[idem]]) * 0.1, 
         labels = sprintf("%s)", letters[ilet]), cex = 2, font = 2,
         adj = 0)
    for (ii in 1:nyints) {
      cuts <- outList[[ii]][[sexes[sx]]]$cuts$nocov
      x <- outList[[ii]][[sexes[sx]]]$x[cuts]
      yy <- outList[[ii]][[sexes[sx]]][[idem]]$nocov[, cuts]
      
      polygon(c(x, rev(x)), c(yy[2, ], rev(yy[3, ])), 
              col = adjustcolor(cols[ii], alpha.f = 0.25), 
              border = NA)
      lines(x, yy[1, ], lty = lty[ii], col = cols[ii], lwd = 2)
    }
    axis(side = 1, at = seq(xlim[1], xlim[2], 1), 
         pos = ylim[[idem]][1], cex.axis = 1.2)
    axis(side = 2, at = seq(ylim[[idem]][1], ylim[[idem]][2], dy[idem]), 
         pos = xlim[1], labels = NA)
    for (jj in 1:2) {
      lines(rep(xlim[jj], 2), ylim[[idem]])
      lines(xlim, rep(ylim[[idem]][jj], 2))
    }
  }
  
}

if (saveResults) dev.off()


# ------------------------------ #
# Parameter posterior densities:
# ------------------------------ #
postMat <- cbind(paramPost, ARpost, LEpost)
nPostCol <- ncol(postMat)

# General categories:
genVars <- c(parNames, sprintf("AR.Age%s", ARages), "LE")
nvars <- length(genVars)
lowVars <- c(-Inf, 0, 0, -Inf, -Inf, -Inf, 0)
names(lowVars) <- genVars

# Layout:
laymat <- rbind(cbind(c(0, rep(1, nvars)), c(0, 1:nvars + 1), 
                      c(nvars + 2, 1:nvars + nvars + 2), 
                      c(2 * nvars + 3, 1:nvars + 2 * nvars + 3)),
                c(0, 0, rep(nvars * 3 + 4, 2)))
widths <- c(0.1, 0.1, 1, 1)
heights <- c(0.1, rep(0.35, nvars), 0.15)
whratio <- sum(widths) / sum(heights)

# Parameter expression:
parExpr <- expression(italic(a)[0], italic(a)[1], italic(c), AR[2], AR[4], 
                      AR[6], italic(e[0]))

# Margins:
mar <- c(3, 1, 1, 1)

# Pdf width:
pdfw <- 6

# Plot:
if (saveResults) {
  pdf(file = "04results/plots/paperPlots/Fig3.pdf", width = pdfw,
      height = pdfw / whratio)
}
layout(laymat, widths = widths, heights = heights)

# y axis label:
par(mar = mar * c(1, 0, 1, 0))
plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
text(0.5, 0.5, "Posterior densities", srt = 90, cex = 1.5)

# Parameter names:
for (ip in 1:nvars) {
  plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.75, 0.85, parExpr[ip], cex = 1.5, xpd = NA, family = "serif")
  
}

isx <- 0
for (sx in sexes) {
  isx <- isx + 1
  par(mar = mar * c(0, 1, 0, 1))
  plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.5, 0.5, sexes[isx], cex = 1.5)
  
  par(mar = mar)
  for (ip in 1:nvars) {
    xlim <- c(NA, NA)
    ylim <- c(0, NA)
    dens <- list()
    for (yi in 1:nyints) {
      varName <- sprintf("%s.%s.Per%s", genVars[ip], sx, yi)
      yy <- postMat[, varName]
      
      qyy <- quantile(yy, probs = c(0.01, 0.99))
      idq <- which(yy >= qyy[1] & yy <= qyy[2])
      muy <- mean(yy[idq])
      sdy <- sd(yy[idq])
      yv <- seq(qyy[1], qyy[2], length = 200)
      q95 <- qtnorm(p = c(0.025, 0.975), mean = muy, sd = sdy, 
                    lower = lowVars[genVars[ip]])
      id95 <- which(yv >= q95[1] & yv <= q95[2])
      dyy <- dtnorm(yv, mean = muy, sd = sdy, lower = lowVars[genVars[ip]])
      dyy <- dyy / max(dyy)
      dd <- list(x = yv, y = dyy, id95 = id95)
      xlim[1] <- min(xlim[1], dd$x, na.rm = TRUE)
      xlim[2] <- max(xlim[2], dd$x, na.rm = TRUE)
      ylim[2] <- max(ylim[2], dd$y, na.rm = TRUE)
      dens[[intslab[yi]]] <- dd
    }
    plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
    for (yi in 1:nyints) {
      xx <- dens[[intslab[yi]]]$x
      yy <- dens[[intslab[yi]]]$y
      id95 <- dens[[intslab[yi]]]$id95
      polygon(c(xx[id95], rev(xx[id95])), c(yy[id95], rep(0, length(id95))), 
              col = adjustcolor(cols[yi], alpha.f = 0.25), border = NA)
      lines(xx, yy, col = cols[yi], lwd = 2)
      
    }
    Axis(x = xlim, side = 1, pos = 0)
  }
}

# Legend:
xseqleg <- seq(0, 0.7, length = nyints)
dleg <- diff(xseqleg)[1]
yleg <- 0.35
par(mar = mar * c(0, 1, 0, 1))
plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
for (yi in 1:nyints) {
  lines(xseqleg[yi] + c(0, dleg * 0.45), rep(yleg, 2), lwd = 2, col = cols[yi])
  text(xseqleg[yi] + dleg * 0.5, yleg, intslab[yi], adj = 0)
}
text(0.5, 1.1, "Variable value", cex = 1.5, xpd = NA)
if (saveResults) dev.off()