# Import necessary modules
library(RcppCNPy)

# Define input variables
k <- 3

# Read in data
setwd("/Users/hirzi/Documents/Hirzi/Cambridge/IndonesiaTrip_Aug22/Workshop/Tutorial/")
admix_Q_raw <- read.table(paste0("GL_75inds.admix.pcangsd.K",k,".admix.",k,".Q"), stringsAsFact=F);
pop_labels <- read.table("SampleInfoAdmix_5inds.txt", header = TRUE)

# Make BarplotQ input
admix_Q_labelled <- cbind(pop_labels[,c(1,2)],admix_Q_raw)
# Rename columns
idx_shift <- 2
for (q in seq(1,(ncol(admix_Q_labelled) - idx_shift))) {
  colnames(admix_Q_labelled)[idx_shift+q] <- paste0("Q",q)
}
# Sort df by "Pop"
admix_Q_labelled_sorted <- admix_Q_labelled[with(admix_Q_labelled, order(Pop, as.integer(sub('\\D+', '', Ind)))),]

# And plot
par(mfrow=c(1,1), mar = c(12, 4, 2, 1))
barplot(t(as.matrix(admix_Q_labelled_sorted[,seq(3,2+k)])), col=c("firebrick","royalblue3", "gold"), border=NA, names.arg=admix_Q_labelled_sorted$Ind, las=2, cex.names=0.85)

