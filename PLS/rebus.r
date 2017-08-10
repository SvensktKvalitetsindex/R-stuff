library("foreign")
library("Rcpp")
library("Amelia")
library("data.table")
library(psych)
# package for pls analysis
library(plspm)
options(OutDec = ",")
#####################################################
# SET WORKING DIRECTORY
#####################################################
setwd("C:\\Users\\Johan\\OneDrive - Svenskt kvalitetsindex\\Current Work\\Sp skåne")
#####################################################

#####################################################
# LOAD FILES - LOCATED IN THE WD
#####################################################
file_name = "B2B.sav"
EPSI.spss <- read.spss(file_name, use.value.labels = FALSE)
# convert list to data frame
EPSI <- as.data.frame(EPSI.spss)
# copy all variable labels in separated list
EPSI_vars <- attr(EPSI.spss, "variable.labels")
# copy all value labels as separated list
EPSI_label <- attr(EPSI.spss, "label.table") #when used, make sure to use the reverse function,
# convert 98 to NA
EPSI[EPSI == 98] <- NA
options(stringsAsFactors = F)
manifest_modell <- read.delim(file = "measurement_model.txt", head = TRUE, sep = "\t")
manifest_label <- manifest_modell$Manifest
q1_names <- read.csv(file = "Q1names.txt", head = TRUE, sep = "\t")

#####################################################
# PREPARE DATA FOR ANALYSIS
#####################################################
# select relevant columns
EPSI_MANIFEST <- EPSI[, manifest_label]
# impute missing value
data_impute = amelia(EPSI_MANIFEST)
xx <- data.table(data_impute$imputations[[1]])
# make sure that that items take values between 1-10
xx[xx > 10] <- 10
xx[xx < 1] <- 1
# data for analysis
EPSI_MANIFEST <- cbind(xx, EPSI[, "Q1"], EPSI[, "CODERESP"])
colnames(EPSI_MANIFEST) <- c(manifest_label, c("Q1", "CODERESP"))
# split data per Q1
df_list <- split(EPSI_MANIFEST, as.factor(EPSI_MANIFEST$Q1))

#####################################################
# MODELL SETUP
#####################################################
# define inner structure
IMAGE = c(0, 0, 0, 0, 0, 0, 0)
EXPECT = c(1, 0, 0, 0, 0, 0, 0)
PRODQ = c(1, 1, 0, 0, 0, 0, 0)
SERVQ = c(1, 1, 1, 0, 0, 0, 0)
VALUE = c(0, 0, 1, 1, 0, 0, 0)
EPSI = c(1, 0, 1, 1, 1, 0, 0)
LOYAL = c(0, 0, 0, 0, 0, 1, 0)
sat_path = rbind(IMAGE, EXPECT, PRODQ, SERVQ, VALUE, EPSI, LOYAL)
sat_blocks = list(which(manifest_modell$Image == -1), which(manifest_modell$Expect == -1), which(manifest_modell$ProdQ == -1), which(manifest_modell$ServQ == -1), which(manifest_modell$Value == -1), which(manifest_modell$EPSI == -1), which(manifest_modell$Loyal == -1)) # Fitness

# vector of modes (reflective indicators)
sat_mod = rep("A", 7)
#####################################################
# RUN PLS AND WRITE TO FILE
#####################################################
sink('analysis-output.txt')

mylist <- paste('comp', 1:length((df_list)), sep = '')
indexlist = list()
for (i in 1:length(df_list)) {
    mylist[i] = list(plspm(df_list[[i]], sat_path, sat_blocks, modes = sat_mod))
    # the scores can be rescaled by (rescale(mylist[[1]])-1)*100/9

    print(q1_names[i, 1])
    print(describe((rescale(mylist[[i]]) - 1) * 100 / 9))
    print(summary(mylist[[i]]))

    index = cbind((rescale(mylist[[i]]) - 1) * 100 / 9, df_list[[i]]$CODERESP)
    index$i <- i
    indexlist[[i]] <- index
}
# Stop writing to the file
big_data = do.call(rbind, indexlist)
colnames(big_data)[8] = "CODERESP"
write.table(big_data, "index.txt", sep = "\t", row.names = FALSE, dec = ",")

sink()

sink('latent_summary.txt')
for (i in 1:length(df_list)) {
    print(q1_names[i, 1])
    print(describe((rescale(mylist[[i]]) - 1) * 100 / 9))
}
sink()

#REBUS
datatemp = df_list[[1]]
sim_global = plspm(datatemp, sat_path, sat_blocks, modes = sat_mod)


# run rebus.pls and choose the number of classes
# to be taken into account according to the displayed dendrogram.
rebus_sim = rebus.pls(sim_global, stop.crit = 0.005, iter.max = 100)
local_rebus = local.models(sim_global, rebus_sim)
sink('rebus analysis-output.txt')
for (i in 1:length(local_rebus)) {
    print(local_rebus[i])
    print(describe((rescale(local_rebus[[i]]) - 1) * 100 / 9))
    print(summary(local_rebus[[i]]))
}
sink()

datafile <- cbind(datatemp$CODERESP, rebus_sim$segments)
colnames(datafile) <- c("CODERESP", "GROUP")
write.table(datafile, "rebus.txt", sep = "\t", row.names = FALSE)
