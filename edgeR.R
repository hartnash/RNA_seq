# code written by Sarah E. Holmes and modified by Mohit Mahey for DGE analysis

# Clearing the environment 
options(error = recover)
rm(list=ls())

# Setting the working directory
setwd("~/MSU REU/poa counts")

# loading required libraries
library(edgeR)

# reading in the raw counts table
x <- read.table("combines_counts.txt", header = TRUE, row.names = "GeneID", sep = '\t')

# checking the counts file and its layout
View(x)

# legend for samples and their respective numbers
#1 = A
#2 = B
#3 = FL
#4 = S1
#5 = S3
#6 = S5
1 = resistant
2 = susceptible

# this is when comparing each with each
group <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6))

# this is when comparing SvR
group <- factor(c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

# this is when we are comparing ALl Sus vs each R
group <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
y <- DGEList(counts = x,group = group)
y$samples

keep <- filterByExpr(y)
summary(keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group-1, data = y$samples)
design

View(design)
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

plotMDS(y)

# comparison of the counts and writing results to .txt files

qlf.AllSvA <- glmQLFTest(fit, contrast = c(-1,0,0,1))
write.table(qlf.AllSvA$table, "AllSvA.txt")
AllSvA <- read.table("AllSvA.txt")
View(AllSvA)

qlf.AllSvB <- glmQLFTest(fit, contrast = c(0,-1,0,1))
write.table(qlf.AllSvA$table, "AllSvB.txt")
AllSvB <- read.table("AllSvB.txt")
View(AllSvB)

qlf.AllSvFL <- glmQLFTest(fit, contrast = c(0,0,-1,1))
write.table(qlf.AllSvFL$table, "AllSvFL.txt")
AllSvFL <- read.table("AllSvFL.txt")
View(AllSvFL)


qlf.SvR <- glmQLFTest(fit, contrast = c(-1,1))
write.table(qlf.SvR$table, "SvR.txt")
SvR <- read.table("SvR.txt")
View(SvR)


qlf.S1vA <- glmQLFTest(fit, contrast = c(-1,0,0,1,0,0))
write.table(qlf.S1vA$table, "S1vA.txt")
S1vA <- read.table("S1vA.txt")
View(S1vA)

qlf.S1vB <- glmQLFTest(fit, contrast = c(0,-1,0,1,0,0))
write.table(qlf.S1vB$table, "S1vB.txt")
S1vB <- read.table("S1vB.txt")
View(S1vB)

qlf.S1vFL <- glmQLFTest(fit, contrast = c(0,0,-1,1,0,0))
write.table(qlf.S1vFL$table, "S1vFL.txt")
S1vFL <- read.table("S1vFL.txt")
View(S1vFL)


qlf.S3vA <- glmQLFTest(fit, contrast = c(-1,0,0,0,1,0))
write.table(qlf.S3vA$table, "S3vA.txt")
S3vA <- read.table("S3vA.txt")
View(S3vA)

qlf.S3vB <- glmQLFTest(fit, contrast = c(0,-1,0,0,1,0))
write.table(qlf.S3vB$table, "S3vB.txt")
S3vB <- read.table("S3vB.txt")
View(S3vB)

qlf.S3vFL <- glmQLFTest(fit, contrast = c(0,0,-1,0,1,0))
write.table(qlf.S3vFL$table, "S3vFL.txt")
S3vFL <- read.table("S3vFL.txt")
View(S3vFL)

qlf.S5vA <- glmQLFTest(fit, contrast = c(-1,0,0,0,0,1))
write.table(qlf.S5vA$table, "S5vA.txt")
S5vA <- read.table("S5vA.txt")
View(S5vA)

qlf.S5vB <- glmQLFTest(fit, contrast = c(0,-1,0,0,0,1))
write.table(qlf.S5vB$table, "S5vB.txt")
S5vB <- read.table("S5vB.txt")
View(S5vB)

qlf.S5vFL <- glmQLFTest(fit, contrast = c(0,0,-1,0,0,1))
write.table(qlf.S5vFL$table, "S5vFL.txt")
S5vFL <- read.table("S5vFL.txt")
View(S5vFL)

qlf.AvB <- glmQLFTest(fit, contrast = c(1,-1,0,0,0,0))
write.table(qlf.AvB$table, "AvB.txt")
AvB <- read.table("AvB.txt")
View(AvB)
