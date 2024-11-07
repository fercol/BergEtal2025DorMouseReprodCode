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

# BaSTA for survival analysis:
if (!"BaSTA" %in% instPacks) {
  install.packages("BaSTA")
}

# paramDemo for demographic functions used after BaSTA:
if (!"paramDemo" %in% instPacks) {
  install.packages("devtools")
  library(devtools)
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

# Logical to save plots:
saveResults <- TRUE

# ===================== #
# ==== DATA PREP.: ====
# ===================== #
# Load BaSTA capture recapture matrix:
bastaDat <- read.csv("DormouseBaSTAdata.csv", header = TRUE,
                     stringsAsFactors = FALSE)

# Sexes:
sexes <- c(f = "Female", m = "Male")

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

# =========================================== #
# ==== BaSTA ANALYSIS BY YEAR INTERVALS: ====
# =========================================== #
# Prepare output list:
outList <- list()
datList <- list()

# Run analysis:
for (yy in 1:nyints) {
  # Find year ranges:
  yv <- ystart[yy]:yend[yy]
  
  # Number of years included in this step:
  ny <- length(yv)
  
  # Subset recapture matrix:
  Ysub <- as.matrix(bastaDat[, sprintf("X%s", yv)])
  
  # Number of recaptures per individual in subset recapture matrix:
  Ysum <- c(Ysub %*% rep(1, ny))
  
  # Find individuals that should be included in period:
  idkeep <- which((bastaDat$Birth < yend[yy] | 
                     bastaDat$Max.Birth < yend[yy]) & 
                    (bastaDat$Death > ystart[yy] | 
                       bastaDat$Min.Death > ystart[yy]) )
  
  # Columns from main matrix to be included in this period:
  idcols <- c(colnames(bastaDat)[1:3], sprintf("X%s", yv), 
              colnames(bastaDat)[c(Ycols[length(Ycols)] + 1):ncol(bastaDat)])
  
  # Period BaSTA data:
  ybastadat <- droplevels(bastaDat[idkeep, idcols])
  
  # Data check:
  ydcheck <- DataCheck(ybastadat, studyStart = ystart[yy], studyEnd = yend[yy])
  
  # Fix any issues with the resulting data:
  yfixdat <- FixCMRdata(ybastadat, studyStart = ystart[yy], studyEnd = yend[yy],
                        autofix = rep(1, 6), silent = TRUE)
  
  # Include resulting data in the data list:
  datList[[intslab[yy]]] <- yfixdat$newData
  
  # Extract period analysis data frame:
  ydat <- yfixdat$newData
  
  # Ensure that the min and max birth and death columns are numeric:
  for (bd in c("Min.Birth", "Max.Birth", "Min.Death", "Max.Death")) {
    ydat[[bd]] <- as.numeric(ydat[[bd]])
  }
  
  # Run BaSTA on subseted data:
  out <- basta(object = ydat, studyStart = ystart[yy], 
               studyEnd = yend[yy], formulaMort = ~ Sex - 1, 
               shape = "bathtub", nsim = 8, parallel = TRUE, ncpus = 8,
               niter = 60000, burnin = 10001, thinning = 50)
  
  # Goodness of fit plots:
  plot(out, plot.type = "gof")
  
  # Store BaSTA output:
  outList[[intslab[yy]]] <- out
}

# ============================= #
# ==== SUMMARY DATA TABLE: ====
# ============================= #
# Create summary data table:
for (iy in 1:nyints) {
  sxNum <- c(f = 0, m = 0)
  for (sx in c("f", "m")) {
    ssx <- c(f = "Sexf", m = "Sexm")[sx]
    sxName <- sexes[sx]
    idsx <- which(datList[[iy]]$Sex == sx)
    sxNum[sx] <- length(idsx)
  }
  tempDf <- data.frame(Period = intslab[iy], Females = sxNum[1], 
                       Males = sxNum[2])
  if (iy == 1) {
    datDf <- tempDf
  } else {
    datDf <- rbind(datDf, tempDf)
  }
}

# Eliminate row names:
rownames(datDf) <- NULL

# Save summary data table:
write.csv(x = datDf, file = "datSummary.csv", row.names = FALSE)

# =========================================== #
# ==== LIFE EXPECTANCY AND AGEING RATES: ====
# =========================================== #
# ----------------------------------------- #
# ---- Kullback Leibler discrepancies: ----
# ----------------------------------------- #
# Output list for Kullback-Leibler discrepancies:
klList <- list()

# Mortality parameter names:
parNames <- c("a0", "a1", "c", "b0", "b1")

# Number of parameters:
nPars <- length(parNames)

# Empty Kullback-Leibler data frame:
klMat <- data.frame(Periods = paste(intslab[c(1, 1, 2)], 
                                    intslab[c(2, 3, 3)], sep = "/"),
                    matrix(0, 3, nPars, dimnames = list(NULL, parNames)))

# Calculate Kullback-Leibler discrepancies to compare parameters per period:
for (sx in c("Sexf", "Sexm")) {
  klList[[sx]] <- list()
  klMat0 <- klMat
  for (ip in 1:nPars) {
    parSx <- sprintf("%s.%s", parNames[ip], sx)
    ii <- 0
    for (yi in 1:(nyints - 1)) {
      for (yj in (yi + 1):nyints) {
        ii <- ii + 1
        kli <- CalcKL(x1 = outList[[yi]]$params[, parSx],
                      x2 = outList[[yj]]$params[, parSx], 
                      length = 10000)
        klMat0[[parNames[ip]]][ii] <- kli["mqKl"]
      }
    }
    klList[[sx]] <- klMat0
  }  
}

# Create Table 1:
kltab <- rbind(data.frame(sex = rep("Female", 3), klList$Sexf),
               data.frame(sex = rep("Male", 3), klList$Sexm))

for (ip in 1:nPars) {
  klcal <- (kltab[[parNames[ip]]] - 0.5) * 2
  kltab[[parNames[ip]]] <- round(klcal, 2)
}
write.csv(kltab, file = "Table1.csv", row.names = FALSE)

# Save results:
save(list = c("outList", "datList", "klList"),
     file = "dormouseMortResultsByPeriods.RData")

# -------------------------- #
# ---- Life expectancy: ----
# -------------------------- #
# Life expectancies per sex and period:
for (sx in 1:2) {
  ssx <- c("Sexf", "Sexm")[sx]
  sxName <- c("Female", "Male")[sx]
  for (iy in 1:nyints) {
    cnames <- colnames(outList[[iy]]$params)
    theMat <- outList[[iy]]$params[, grep(ssx, cnames)]
    colnames(theMat) <- parNames
    exv <- apply(theMat, 1, function(thei) {
      CalcRemainLifeExp(theta = thei, x = 0, dx = 0.001, xmax = 20, 
                        model = "GO", shape = "bathtub")["RemLExp"]
    })
    exstats <- c(Mean = mean(exv), SD = sd(exv), Lower = c(quantile(exv, 0.025)),
                 Upper = c(quantile(exv, 0.975)))
    names(exstats) <- c("Mean", "SE", "Lower", "Upper")
    temp <- data.frame(Sex = sxName, Period = intslab[iy], t(exstats))
    if (sx == 1 & iy == 1) {
      lifeExp <- temp
    } else {
      lifeExp <- rbind(lifeExp, temp)
    }
  }
}

# Save Table S2:
write.csv(lifeExp, file = "TableS2.csv", row.names = FALSE)

# ======================= #
# ==== PLOT RESULTS: ====
# ======================= #
# ----------------------------- #
# Mortality and survival plots:
# ----------------------------- #
# id vector of periods:
lty <- 1:nyints

# Colors per period:
cols <- brewer.pal(4, "Set1")

# Plotting x and y limits:
xlim <- c(0, 6)
ylim <- list(surv = c(0, 1), mort = c(0, 3))

# Age vector:
xv <- seq(0, 10, 0.001)

# y-axis interval lengths:
dy <- c(surv = 0.2, mort = 0.5)

# Sex labels:
sexes <- c(f = "Female", m = "Male")

# Densities for polygon plots:
dens <- seq(20, 30, length = nyints)

# Angles for density plots:
angle <- seq(45, 135, length = nyints)

# Y axis names:
yaxname <- c(surv = "Survival", mort = "Hazard rate")

# Layout matrix:
laymat <- cbind(c(0, 2, 3, 0), c(4, 5, 6, 1), c(7, 8, 9, 1))

# Layout widths and heights:
widths <- c(0.1, 1, 1)
heights <- c(0.1, 0.75, 0.75, 0.1)

# Margins:
mar <- c(3, 1, 1, 1)

# Plot width in inches:
plw <- 8

# Plot height in inches:
plh <- plw / (sum(widths) / sum(heights))

# Label and axis sizes:
cex.lab <- 2
cex.ax <- 1.5

# Produce Figure 2:
if (saveResults) {
  pdf(file = "Fig2.pdf", width = plw, height = plh)
}

layout(laymat, widths = widths, heights = heights)

# X-axis label:
par(mar = mar * c(0, 1, 0, 1))
plot(xlim, c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
text(mean(xlim), 0.85, "Age (years)", cex = cex.lab, xpd = NA)
legend("right", intslab, lty = lty, col = cols, pch = NA, bty = 'n', lwd = 4, 
       seg.len = 8)

# Y-axis label:
for (idem in c("surv", "mort")) {
  par(mar = mar * c(1, 0, 1, 0))
  plot(c(0, 1), ylim[[idem]], col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.25, mean(ylim[[idem]]), yaxname[idem], srt = 90, cex = cex.lab)
  axis(side = 4, at = seq(ylim[[idem]][1], ylim[[idem]][2], dy[idem]), 
       pos = 0.5, las = 2, lwd = NA)
  
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
         labels = sprintf("%s)", letters[ilet]), cex = 1.5, font = 2,
         adj = 0)
    for (ii in 1:nyints) {
      x <- outList[[ii]]$x[outList[[ii]]$cuts[[sxname]]]
      mux <- outList[[ii]][[idem]][[sxname]][, outList[[ii]]$cuts[[sxname]]]
      
      polygon(c(x, rev(x)), c(mux[2, ], rev(mux[3, ])), 
              col = adjustcolor(cols[ii], alpha.f = 0.25), 
              border = NA)
      lines(x, mux[1, ], lty = lty[ii], col = cols[ii], lwd = 2)
    }
    axis(side = 1, at = seq(xlim[1], xlim[2], 2), 
         pos = ylim[[idem]][1])
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
# Layout matrix:
laymat <- rbind(cbind(c(0, rep(1, nPars)), c(0, 1:nPars + 1), 
                      c(nPars + 2, 1:nPars + nPars + 2), 
                      c(2 * nPars + 3, 1:nPars + 2 * nPars + 3)),
                c(0, 0, rep(nPars * 3 + 4, 2)))

# Widths, heights and width-height ratio:
widths <- c(0.1, 0.1, 1, 1)
heights <- c(0.1, rep(0.5, nPars), 0.15)
whratio <- sum(widths) / sum(heights)

# Parameter expression:
parExpr <- expression(italic(a)[0], italic(a)[1], italic(c), 
                      italic(b)[0], italic(b)[1])

# Margins:
mar <- c(3, 1, 1, 1)

# Pdf width:
pdfw <- 6

# Produce Figure 4:
if (saveResults) {
  pdf(file = "Fig4.pdf", width = pdfw, height = pdfw / whratio)
}
layout(laymat, widths = widths, heights = heights)

# y axis label:
par(mar = mar * c(1, 0, 1, 0))
plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
text(0.5, 0.5, "Parameter posterior densities", srt = 90, cex = 1.5)

# Parameter names:
for (ip in 1:nPars) {
  plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.75, 0.85, parExpr[ip], cex = 1.5, xpd = NA)
  
}

isx <- 0
for (sx in c("Sexf", "Sexm")) {
  isx <- isx + 1
  par(mar = mar * c(0, 1, 0, 1))
  plot(c(0, 1), c(0, 1), col = NA, xlab = "", ylab = "", axes = FALSE)
  text(0.5, 0.5, c("Females", "Males")[isx], cex = 1.5)
  
  par(mar = mar)
  for (ip in 1:nPars) {
    parSx <- sprintf("%s.%s", parNames[ip], sx)
    xlim <- c(NA, NA)
    ylim <- c(0, NA)
    dens <- list()
    for (yi in 1:nyints) {
      dd <- density(outList[[yi]]$params[, parSx])
      xlim[1] <- min(xlim[1], dd$x, na.rm = TRUE)
      xlim[2] <- max(xlim[2], dd$x, na.rm = TRUE)
      ylim[2] <- max(ylim[2], dd$y, na.rm = TRUE)
      dens[[intslab[yi]]] <- dd
    }
    plot(xlim, ylim, col = NA, axes = FALSE, xlab = "", ylab = "")
    for (yi in 1:nyints) {
      xx <- dens[[intslab[yi]]]$x
      yy <- dens[[intslab[yi]]]$y
      polygon(c(xx, rev(xx)), c(yy, rep(0, length(yy))), 
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
text(0.5, 1.1, "Parameter value", cex = 1.5, xpd = NA)
if (saveResults) dev.off()
