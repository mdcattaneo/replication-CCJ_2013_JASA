######################################################################################
## R CODE FOR CATTANEO-CRUMP-JANSSON (2012, JASA): Generalized Jackknife for WAD
## DATE: 24-Jan-2013
## NOTE: Replicates empirical illustration of the paper
## FILES REQUIRED:
##    1) replic-empillus.R (this file): main file to run results
##    2) replic-empillus-functions.R  : our R functions
##    3) C_Kernel_WAD-jack.so         : C object file with WAD estimation functions
##    4) replic-finaldata.dta         : final database used in empirical illustration
######################################################################################

## SET PATH TO YOUR WORKING DIRECTORY HERE:
setwd("/mainfiles/cattaneo/research/GeneralizedJackknife/shared-code/2013_Kourtellos")

## CLEAN R ENVIRONMENT & LOAD R PACKAGES (you may need to install them first!)
rm(list=ls(all=TRUE))
library(foreign); library(MASS); library(gam); library(Hmisc)

## LOAD OUR R FUNCTIONS
source("replic-empillus-functions.R")

## LOAD OUR C FUNCTIONS (choose appropriate architecture!)
#if (is.loaded("phi", PACKAGE="C_Kernel_WAD-jack")) dyn.unload("C_Kernel_WAD-jack.so"); dyn.load("C_Kernel_WAD-jack.linux64.so")
if (is.loaded("phi", PACKAGE="C_Kernel_WAD-jack")) dyn.unload("C_Kernel_WAD-jack.so"); dyn.load("C_Kernel_WAD-jack.win64.so")

######################################################################################
## CLEAN ORIGINAL DATA
## NOTE: This portion of the code is commented out because the file aer-2008-1056.dta
##       is too big. We leave the code here for replication purposes. The original
##       database (aer-2008-1056.dta) may be downloaded from the AER website.
######################################################################################
### READ ORIGINAL DATA
#data.orig = read.dta("../../empiricalillustration/2008-0156-data/aer-2008-1056.dta")
### CLEAN DATA
#data = data.orig[data.orig[,"hisp"]==0 & data.orig[,"black"]==0 & data.orig[,"sex79"]==1 &
#                 data.orig[,"education"]>=12 & data.orig[,"education"]<=16 &
#                 data.orig[,"age79"]>=15 & data.orig[,"age79"]<=19 &
#                 !is.na(data.orig[,"age79"]) & !is.na(data.orig[,"education"]) &
#                 !is.na(data.orig[,"lwageNJ2"]) & !is.na(data.orig[,"AFQT89adfin"]) &
#                 !is.na(data.orig[,"school_enrollment"]) & !is.na(data.orig[,"teacher_salary"]),]
#data = data[,c("lwageNJ2","AFQT89adfin","school_enrollment","teacher_salary")]; dim(data)
#data = data[data[,"school_enrollment"]<=5000 & data[,"teacher_salary"]<=15000,]; dim(data)
### SCALE VARIABLES FOR INTERPRETATION (WLOG)
#data[,"school_enrollment"] = data[,"school_enrollment"] / 100
#data[,"teacher_salary"] = data[,"teacher_salary"] / 1000
### SAVE FINAL DATABASE
#write.dta(data,"replic-finaldata.dta")

########################################
## READ FINAL DATABASE & SUMMARY STATS
########################################
data = read.dta("replic-finaldata.dta"); summary(data)

########################################
## PLOT E[lwageNJ2 | AFQT89adfin]
########################################
pdf(file=paste("figure.pdf",sep=""), paper="special", width=8.5, height=3)
  par(mfrow=c(1,3))
  plot(gam(lwageNJ2 ~ s(AFQT89adfin), data=data), se=T, ylab="Log Wages", xlab="AFQT (standarized)", ylim=c(-1,.6), xlim=c(-2,2))
  plot(gam(lwageNJ2 ~ s(school_enrollment), data=data), se=T, ylab="Log Wages", xlab="School size (per 100)", ylim=c(-.5,.5), xlim=c(0,40))
  plot(gam(lwageNJ2 ~ s(teacher_salary), data=data), se=T, ylab="Log Wages", xlab="Teacher salary (per 1,000)", ylim=c(-.2,.6), xlim=c(7,14))
dev.off()

########################################
## Estimation
########################################
## CENTER & SCALE VARIABLES
data[,"school_enrollment"] = (data[,"school_enrollment"] - mean(data[,"school_enrollment"])) / sd(data[,"school_enrollment"])
data[,"teacher_salary"]    = (data[,"teacher_salary"] - mean(data[,"teacher_salary"])) / sd(data[,"teacher_salary"])

## PARAMETRIC FITS
lm.lin = lm(lwageNJ2 ~ AFQT89adfin + school_enrollment + teacher_salary, data=data); summary(lm.lin)
lm.sq  = lm(lwageNJ2 ~ AFQT89adfin + I(AFQT89adfin^2) + school_enrollment + teacher_salary, data=data); summary(lm.sq)
lm.cub = lm(lwageNJ2 ~ AFQT89adfin + I(AFQT89adfin^2) + I(AFQT89adfin^3)+ school_enrollment + teacher_salary, data=data); summary(lm.cub)

## WADE Estimation
## Grid Values, Constants & Tables for results
c.hat = matrix(rep(c(.90,.95,1,1.05,1.10),3),ncol=3)
P=4; d=3; CB = -0.160354; CS = 1/8; CH = -3/(8*pi^(3/2)); n=nrow(data)
col.names=c("s","n","d","P","c","h"           ,"theta1.hat.h","sigma1.hat.h","obs.trimmed.h",
                                "h1","h2","h3","theta1.hat.hs","sigma1.hat.hs","obs.trimmed.hs")
table.hat.1d = matrix(NA, nrow=nrow(c.hat), ncol=length(col.names), dimnames=list(NULL,col.names))
table.hat.tr = matrix(NA, nrow=nrow(c.hat), ncol=length(col.names), dimnames=list(NULL,col.names))
table.hat.1d.under = matrix(NA, nrow=nrow(c.hat), ncol=length(col.names), dimnames=list(NULL,col.names))
table.hat.tr.under = matrix(NA, nrow=nrow(c.hat), ncol=length(col.names), dimnames=list(NULL,col.names))

## ROT bandwidths
sigmas.hat = apply(data[,c("AFQT89adfin", "school_enrollment", "teacher_salary")],2,scale)
bs.hat = lm(lwageNJ2 ~ AFQT89adfin + school_enrollment + teacher_salary, data=data)$coeff[2:4]

h.rot.1d = c(ROT.h.1d(s=sigmas.hat, CB, CS, CH, P, n), rep(NA,d))
h.rot.tr = c(ROT.h.tr(s=sigmas.hat, b=bs.hat, CB, CS, CH, P, n), ROT.hs.tr(s=sigmas.hat, b=bs.hat, CB, CS, CH, P, n))

## Triming
trim = apply(data[,c("AFQT89adfin", "school_enrollment", "teacher_salary")], 2, function(x) quantile(x, c(.01,.99)))
trim = apply(trim, 2, function(x) min(abs(x)))

## Results
for (c.i in 1:nrow(c.hat)) {
    #### USING ROT BANDWIDTH CHOICES
    table.hat.1d[c.i,1:5] = c(1,n,d,P,c.hat[c.i])
    table.hat.tr[c.i,1:5] = c(1,n,d,P,c.hat[c.i])

    ## M1 - COMMON BANDWIDTH (two ROT bandwidth choices: 1d & tr)
    h = c.hat[c.i,]*h.rot.1d[1]
    out.C = C.wad(y=data[,"lwageNJ2"],
                  x=as.matrix(data[,c("AFQT89adfin","school_enrollment","teacher_salary")]),
                  d=d,trim=trim,P=P,n=n,h=h)
    table.hat.1d[c.i,6:9] = c(h[1],out.C$theta1.hat,out.C$sigma1.hat,n-sum(out.C$used))

    h = c.hat[c.i,]*h.rot.tr[1]
    out.C = C.wad(y=data[,"lwageNJ2"],
                  x=as.matrix(data[,c("AFQT89adfin","school_enrollment","teacher_salary")]),
                  d=d,trim=trim,P=P,n=n,h=h)
    table.hat.tr[c.i,6:9] = c(h[1],out.C$theta1.hat,out.C$sigma1.hat,n-sum(out.C$used))

    ## M1 - DIFFERENT BANDWIDTHS
    h = c.hat[c.i,]*h.rot.tr[2:4]
    out.C = C.wad(y=data[,"lwageNJ2"],
                  x=as.matrix(data[,c("AFQT89adfin","school_enrollment","teacher_salary")]),
                  d=d,trim=trim,P=P,n=n,h=h)
    table.hat.tr[c.i,10:15] = c(h,out.C$theta1.hat,out.C$sigma1.hat,n-sum(out.C$used))

    #### USING UNDERSMOOTHING (0.90*h) ROT BANDWIDTH CHOICES
    under = 0.90
    table.hat.1d.under[c.i,1:5] = c(1,n,d,P,c.hat[c.i])
    table.hat.tr.under[c.i,1:5] = c(1,n,d,P,c.hat[c.i])

    ## M1 - COMMON BANDWIDTH (two ROT bandwidth choices: 1d & tr)
    h = c.hat[c.i,]*under*h.rot.1d[1]
    out.C = C.wad(y=data[,"lwageNJ2"],
                  x=as.matrix(data[,c("AFQT89adfin","school_enrollment","teacher_salary")]),
                  d=d,trim=trim,P=P,n=n,h=h)
    table.hat.1d.under[c.i,6:9] = c(h[1],out.C$theta1.hat,out.C$sigma1.hat,n-sum(out.C$used))

    h = c.hat[c.i,]*under*h.rot.tr[1]
    out.C = C.wad(y=data[,"lwageNJ2"],
                  x=as.matrix(data[,c("AFQT89adfin","school_enrollment","teacher_salary")]),
                  d=d,trim=trim,P=P,n=n,h=h)
    table.hat.tr.under[c.i,6:9] = c(h[1],out.C$theta1.hat,out.C$sigma1.hat,n-sum(out.C$used))

    ## M1 - DIFFERENT BANDWIDTHS
    h = c.hat[c.i,]*under*h.rot.tr[2:4]
    out.C = C.wad(y=data[,"lwageNJ2"],
                  x=as.matrix(data[,c("AFQT89adfin","school_enrollment","teacher_salary")]),
                  d=d,trim=trim,P=P,n=n,h=h)
    table.hat.tr.under[c.i,10:15] = c(h,out.C$theta1.hat,out.C$sigma1.hat,n-sum(out.C$used))
}

## Construct GENJACK estimates
table.h.1d.genjack = genjack(data=table.hat.1d[, c("s","n","c","h","theta1.hat.h","sigma1.hat.h")],
                             d=3, nc.down=1, nc.up=0, thetalab="theta1.hat.h")

table.h.tr.genjack = genjack(data=table.hat.tr[, c("s","n","c","h","theta1.hat.h","sigma1.hat.h")],
                             d=3, nc.down=1, nc.up=0, thetalab="theta1.hat.h")

table.hs.tr.genjack = genjack(data=table.hat.tr[, c("s","n","c","h1","h2","h3","theta1.hat.hs","sigma1.hat.hs")],
                              d=3, nc.down=1, nc.up=0, thetalab="theta1.hat.hs")


table.h.1d.genjack.under = genjack(data=table.hat.1d.under[, c("s","n","c","h","theta1.hat.h","sigma1.hat.h")],
                             d=3, nc.down=1, nc.up=0, thetalab="theta1.hat.h")

table.h.tr.genjack.under = genjack(data=table.hat.tr.under[, c("s","n","c","h","theta1.hat.h","sigma1.hat.h")],
                             d=3, nc.down=1, nc.up=0, thetalab="theta1.hat.h")

table.hs.tr.genjack.under = genjack(data=table.hat.tr.under[, c("s","n","c","h1","h2","h3","theta1.hat.hs","sigma1.hat.hs")],
                              d=3, nc.down=1, nc.up=0, thetalab="theta1.hat.hs")


## Construct LATEX table results
table.tex = rbind(c(table.h.1d.genjack[3,c("theta1.hat.h","theta1.hat-jack")], table.h.1d.genjack[3,"sigma1.hat.h"]/sqrt(n)),
                  c(table.h.1d.genjack.under[3,c("theta1.hat.h","theta1.hat-jack")], table.h.1d.genjack.under[3,"sigma1.hat.h"]/sqrt(n)),
                  c(table.h.tr.genjack[3,c("theta1.hat.h","theta1.hat-jack")], table.h.tr.genjack[3,"sigma1.hat.h"]/sqrt(n)),
                  c(table.h.tr.genjack.under[3,c("theta1.hat.h","theta1.hat-jack")], table.h.tr.genjack.under[3,"sigma1.hat.h"]/sqrt(n)),
                  c(table.hs.tr.genjack[3,c("theta1.hat.hs","theta1.hat-jack")], table.hs.tr.genjack[3,"sigma1.hat.hs"]/sqrt(n)),
                  c(table.hs.tr.genjack.under[3,c("theta1.hat.hs","theta1.hat-jack")], table.hs.tr.genjack.under[3,"sigma1.hat.hs"]/sqrt(n)))

print(table.tex,3)

latex(round(table.tex,3), file=paste("table.empillus.txt", sep=""), append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=c(2,1), cgroup=c("Coef.","Std. Err."),
      colheads=c("$\\mathbf{\\hat\\theta}_n(\\mathbf{\\hat{H}}_n)$", "$\\mathbf{\\tilde\\theta}_n(\\mathbf{\\hat{H}}_n,\\mathbf{c})$",
                 "$\\mathbf{\\hat\\Sigma}_n(\\mathbf{\\hat{H}}_n)$"),
      r.cgroup=c(1,1,1), rgroup=c("Common Bandwidth: ROT-1d", "Common Bandwidth: ROT-tr", "Different Bandwidths: ROT-tr"),
      rowname=c(paste("$\\mathbf{\\hat{H}}_n = ",r(table.h.1d.genjack[3,"h"]),"\\cdot\\mathbf{I}_3$",sep=""),
                paste("$\\mathbf{\\hat{H}}_n = ",under,"\\cdot",r(table.h.1d.genjack[3,"h"]),"\\cdot\\mathbf{I}_3$",sep=""),
                paste("$\\mathbf{\\hat{H}}_n = ",r(table.h.tr.genjack[3,"h"]),"\\cdot\\mathbf{I}_3$",sep=""),
                paste("$\\mathbf{\\hat{H}}_n = ",under,"\\cdot",r(table.h.tr.genjack[3,"h"]),"\\cdot\\mathbf{I}_3$",sep=""),
                paste("$\\mathbf{\\hat{H}}_n = \\limfunc{diag}(",r(table.hs.tr.genjack[3,"h1"]),",",
                                                          r(table.hs.tr.genjack[3,"h2"]),",",
                                                          r(table.hs.tr.genjack[3,"h3"]),")$",sep=""),
                paste("$\\mathbf{\\hat{H}}_n = ",under,"\\cdot","\\limfunc{diag}(",r(table.hs.tr.genjack[3,"h1"]),",",
                                                          r(table.hs.tr.genjack[3,"h2"]),",",
                                                          r(table.hs.tr.genjack[3,"h3"]),")$",sep=""))
      )



