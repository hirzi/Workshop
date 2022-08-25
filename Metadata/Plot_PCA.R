# Import necessary modules
library(methods)
library(ggplot2)
library(plotly)
library(RcppCNPy)
library(scales)
library(ggpubr)

# Set working directory (change this to your own)
setwd("/Users/hirzi/Documents/Hirzi/Cambridge/IndonesiaTrip_Aug22/Workshop/Tutorial/")

# Define number of principle components to plot
components<-"1-2"

# Read input file
covar <- read.table("GL_75inds.pcangsd.cov", stringsAsFact=F);

# Read annot file
annot <- as.data.frame(read.table("SampleInfoPCA_5inds.txt", sep="\t", header=T));

# Parse components to analyze
comp <- as.numeric(strsplit(components, "-", fixed=TRUE)[[1]])

# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
PC$Region <- factor(annot$IID)
PC$Metaregion <- factor(annot$REGION)
PC$ID <- factor(annot$FID)
PC$COV <- as.numeric(as.character(annot$DEPTH))
PC$LAT <- as.numeric(as.character(annot$LAT))
PC$LON <- as.numeric(as.character(annot$LON))
cols_regions <- c("Balkans" = 16, "Julian_Alps" = 3, "Apennines" = 17, "Central_Alps" = 15, "French_Alps+Jura" = 8, "French_Alps_Longicaulis" = 4, "Monte_Baldo+Dolomites" = 11, "Calabria" = 5)
#shape_regions <- c("Balkan" = 19, "Apennine" = 15, "Alpine" = 17)
shape_regions <- c("Balkan" = 21, "Apennine" = 22, "Alpine" = 24)
col_regions <- c("Balkan" = "darkseagreen2", "Apennine" = "orange", "Alpine" = "darkred")
col2_regions <- c("Balkan" = "grey5", "Apennine" = "grey5", "Alpine" = "grey5")

par(mar=c(5.1,4.1,8,2.1))
par(mfrow=c(1,1))

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

pl <- ggplot() + geom_point(data=PC, size = 6, aes_string(x=x_axis, y=y_axis, color=PC$Pop, shape = PC$Region)) + theme_bw() + scale_colour_hue(name = "POPULATION") + scale_shape_manual(name = "REGION", values = cols_regions) + ggtitle(title)
plot(pl)

# Make an interactive plotly plot
pl_plotly <- ggplotly(pl)
pl_plotly
