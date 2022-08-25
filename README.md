# Inferring population structure & demography from low-coverage whole genome sequencing data

  - [Introduction](#introduction)
  - [Workshop - initial preparation](#workshop---initial-preparation)
  - [Workshop - data](#workshop---data)
  - [Workshop - programs we will use](#workshop---programs-we-will-use)
  - [Instructions - preparing our Docker containers](#instructions---preparing-our-docker-containers)
      - [1. Open Docker](#step-1-open-docker)
      - [2. Running Ubuntu on Docker](#step-2-running-ubuntu-on-docker)
      - [3. Pulling Docker images](#step-3-pulling-docker-images)
      - [4. Creating and mounting a volume](#step-4-creating-and-mounting-a-volume)
      - [5. Download data](#step-5-download-data)
      - [6. Getting a hand with command-line in Docker](#step-6-getting-a-hand-with-command\-line-in-docker)
      - [7. Index files](#step-7-index-files)
  - [Principle component analysis (PCA)](#principle-component-analysis)
  - [Admixture and ancestry analysis](#Admixture-and-ancestry-analysis)
  - [Site frequency spectrum and summary statistics](#site-frequency-spectrum-and-summary-statistics)
      - [One population](#one-population)
      - [Two populations](#two-populations)
      - [Three populations](#three-populations)

<br>

# Introduction

## Why go low-coverage?

**Experimental design** All scientific lines of inquiry start with a question. From this question, we (the researcher) try to come up with an experimental design to best address said question. If money and time were no object, we may e.g. plan for an experiment with 10 treatments, 10 replicates per treatment, and 1000 samples per replicate. However, in the real world, the experimental design is not determined purely by how best to address the biological question at hand, but also by *cost, time and technical feasibility*.

Typically, we are faced with the following trade-off; of having either **1)** more samples at the expense of data per sample or **2)** less samples but with more data per sample.

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Experimental_design.png" width="650"> 

With respect to sequencing and genetic data, if we start with the perfect or complete representation of a unit data as the whole genome sequenced at high coverage (e.g. 50x), less data can imply one of two things: **1)** sequencing a reduced or sub-representation of the genome, i.e. using genetic markers like microsatellites or limited sets of single nucleotide polymorphisms (SNPs) or **2)** sequencing the whole genome but at low coverage (depth). I.e. a trade-off of breadth vs depth. To give a concrete example, imagine you had enough money to sequence 1 million reads, and that this is sufficient to sequence your whole genome at 2x or 10% of your genome at 20X, *which would you choose*? 

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/breadth_vs_depth.png" width="650"> 

The answer can be difficult, and lies in weighing the respective advantages and disadvantages of these two alternatives. Briefly, both approaches have their caveats. For the former (i.e. subset of genetic markers), we **1)** assume the sub-selection of the genome to be representative of the whole genome, **2)** are prone to ascertainment bias, and **3)** lack data in unsequenced parts of the genome which prohibits detection of new, potentially interesting genetic variants. For the latter (low-coverage sequencing), our certainty in the genotype call (i.e. whether something is **A**, **C**, **T** or **G**) is low, because we read each position fewer times, and hence are prone to more sequencing errors in our genotype calls. Furthermore, the power to detect rare variants is strongly diminished.

That said, in the last years (for both cases), there has been notable advances in alleviating these respective downsides. For the former, modern sequencing techniques provide more markers at increasingly cheaper prices, while for the latter, we now have methods that explicitly accommodate the uncertainty in genotype calls in analyses (i.e., rather than having to explicitly call genotypes, we consider their genotype likelihoods instead, aka the probability of the data given a specific genotype).

Low-coverage methods also lend themselves well to the sequencing and analysis of ancient DNA, where high-coverage, high-quality DNA sequences may not be easily attainable.


### Summary
Through working with genotype likelihoods rather than relying on discrete (lossy) genotype calls, low-coverage methods are useful for they allow for the statistical propagation of uncertainty from raw sequencing data to downstream analysis. The effect may be huge for data at low-coverage or negligable for data at high-coverage (where analytical results will tend to converge between low-coverage methods and classical genotype-call based methods). In addition, working directly with genotype likelihoods typically involves fewer (potentially lossy) processing steps, e.g. genotype calling and various filtering.

<br>

# Workshop - initial preparation

In this session you will learn how to use low-coverage whole genome data to infer population structure and admixture (ancestry), via two methods:

  - **Principal Components Analysis (PCA)**
  - **Admixture analysis**

Population genetic analyses of NGS data is typically run on large linux-based computing clusters. For this workshop, since we do not have access to this, we will be running population genetic analyses in **Docker**. 

### What is Docker?
Docker is an open source container-based technology that runs in an isolated, self-contained package that can be efficiently distributed and executed in a portable manner across a wide range of computing platforms. Containerisation in concept is very similar to virtualisation, i.e. a method of isolating an application from the underlying machine operating system. The difference between virtual machines and containers is that containers do not require a full operating system to operate but rather the application and dependencies, means they are much smaller in size (Gharib 2018). 

## Step 1. Make sure you have Docker Desktop installed on your computer
Instructions for intalling Docker can be found [here](https://www.docker.com/get-started/), with OS-specific information (including minimum system requirements) and download links for [Windows](https://docs.docker.com/desktop/install/windows-install/), [Mac](https://docs.docker.com/desktop/install/mac-install/) and [Linux](https://docs.docker.com/desktop/linux/).

## Step 2. Make sure you are familiar with the basics of bash shell scripting
We will be working almost exclusively through the command-line in docker, so if you have not used shell scripting before or are a little rusty with it, it may be helpful to have a look at a tutorial like [this](https://linuxconfig.org/bash-scripting-tutorial-for-beginners) or a cheat sheet like [this](https://bioinformaticsworkbook.org/Appendix/Unix/UnixCheatSheet.html#gsc.tab=0) before proceeding with this workshop.

<br>

# Workshop - data

We will be working with low-coverage whole-genome sequencing (WGS) data (average 2x) of a plant species *Dianthus sylvetris* (Wood pink). This is a perennial plant species that grows throughout the mountain ranges of Europe (inc. the Alps, Apennines & Dinarides). Given the European mountain ranges experienced repeated bouts expansion and recession of glacial ice sheets during the last ca. 2 million years (the Quaternary glaciations or "ice ages"), this species likely experienced a complex and dynamic demographic history. Additionally, this species inhabits a large elevational range (0 - 2500 meters), with the consequence that contemporary populations exhibit a remarkable degree of local adaptation in phenotypic and life history traits.

Add map!

The populations and data that we will use in this workshop represents a small subset from a larger study that covered the geographic and ecological range of the species (https://www.biorxiv.org/content/10.1101/2022.06.07.495159v1). Our data is in BAM format (i.e. mapped sequencing data in binary format) and span a ~2MB region from four scaffolds selected at random across the genome. (may want to extract an interesting region for FST/PBS analysis, 15 inds pops. alternatively, consider taking data directly from Simons's data and do FST scans). We will use this data to infer this species' population structure first via principle component analysis (PCA) and then via admixture (ancestry) analysis.

<br>

# Workshop - programs we will use
For this practical, we will be using:

[SAMtools](https://www.htslib.org/)

[ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next
Generation Sequencing Data)

[PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd)

We will run these programs in Docker containers so there is no need to install these programs locally on your computer.

<br>

# Instructions - preparing our Docker containers

## Step 1\. Open Docker
First, open Docker Desktop on your computer. Then, open Windows Terminal (Powershell), Mac Terminal or Linux terminal. If Docker has ben succesfully installed, you should see the Docker help menu after typing "docker" in the terminal. Docker itself has a whole syntax for usage ([Docker manual](https://docs.docker.com/engine/reference/commandline/docker/), but today we will focus more on the bash-based command-lines as that is what is mainly used to work with NGS data (again here Docker is only used as container to run NGS programs without installing them locally).

## Step 2\. Running Ubuntu on Docker
We first need a Linux distribution (e.g. Ubuntu) running in our Docker container. Use "docker pull" to pull a Docker image or repository from a registry. 

	docker pull ubuntu

This will pull the latest version of Ubuntu to Docker.

To run Ubuntu, we use "docker run"

	docker run -it --rm ubuntu

The -it instructs Docker to allocate a pseudo-TTY connected to the container’s stdin; creating an interactive bash shell in the container. I.e. this allows us to run Ubuntu interactively in Docker. --rm automatically remove the container when it exits. To exit the bash interactive shell, enter: exit 13.

## Step 3\. Pulling Docker images

Let's now "pull" Docker images for the other programs we will be using.

	docker pull biocontainers/samtools:v1.9-4-deb_cv1
	docker pull zjnolen/angsd
	docker pull didillysquat/pcangsd

## Step 4\. Creating and mounting a volume
We have now "pulled" the necessary programs into Docker. We now need the data. First we'll need to create a volume in Docker to store downloaded and generated data. Here, we will run an Ubuntu container with a named volume via the --mount argument (i.e. we will name this created volume "myvol" and define the associated container path as "/data". The -w allows the command to be executed inside the defined working directory. Before we download data, we will first need to install some other programs within the running container, namely git and wget. 

	docker run --name base --mount source=myvol,target=/data -w /data -it --rm ubuntu
	apt update
	apt install git
	apt-get install wget

While we are at it, let's install a text editor too (vim nano).

	apt-get install vim nano

## Step 5\. Download data
Then we download data. All the BAM files as well as population metadata are deposited in the github: https://github.com/hirzi/Workshop (make public!). Let's pull this directory to our Docker container and our named volume.

	git clone https://github.com/hirzi/Workshop.git

Once downloaded, let's navigate to the /data directory containing the samples. Once there and before we continue, let's download the reference sequence. The reference is deposited in a separate repository (because it's large file).

	wget https://www.dropbox.com/s/1nggjjhrcjseuwx/assembly_homozygous.fa?dl=1

Let's rename this reference sequence file.
	
	mv 'assembly_homozygous.fa?dl=1' assembly_homozygous.fa

## Step 6\. Getting a hand with command-line in Docker

Now that we have all the data downloaded, let's try a few bash commands to see what we have.

	cd /data/Workshop/Data/
	ls
	ls -lrth *bam

How many BAM files (i.e. samples) are there? (hint: ls \*bam | wc -l). In total we have 95 samples (13 populatiosn with 5 individuals each and 2 populations with 15 individuals each).

Now, let's check directory permissions. Do we have permission to write to the /data/Workshop/Data/ directory? You can check this with the command "ls -l". If not, let's change permissions so that you can write to this directory.

	chmod 777 /data/Workshop/Data

## Step 7\. Index files
We now have the data stored in a named volume. Even when we close the Docker container, the volume is maintained. We can always re-access the volume by mounting it it when running a new container. Before we perform any analysis with BAM and fasta files (or really any large sequencing file), we need to index them so that programs can parse through them efficiently. To do this, let's open another terminal tab, in which we will run Samtools. Note the named volume that we mount to.

	docker run --name samtools --mount source=myvol,target=/data -w /data/ -it --rm biocontainers/samtools:v1.9-4-deb_cv1
	
Index the fasta file (i.e. the reference sequence):

	samtools faidx assembly_homozygous.fa

Index all the BAM files (here we'll do this in a for loop):

	for i in $(ls *bam); do samtools index ${i}; done

!!!!!!!!!!!!!!! Make sample list (let's actually have this in the Github, in case htere are difference with file ordering)
!!!!!!!!!!!!!!!#ls *bam > samples.list

<br>

# Principle component analysis

Principle component analysis (PCA) is a method that reduces the dimensionality of a dataset to emphasize variation and thus bring out strong (and more easily) interpretable patterns in a dataset. It does so by transforming a large set of (potentially correlated) variables into a smaller set of uncorrelated variables while minimising information loss.

Have a look [here](https://setosa.io/ev/principal-component-analysis) to intuitively see how this is done.

For NGS data, each variant site (e.g. single nucleotide polymorphism (SNP)) represents one dimesion (or axis) of variation. Thus, if we have 1 million SNPs, we have 1 million dimensions (or axes) of variation. Obviously, this would be extremely difficult to visualise without some sort of transformation applied (we can normally only visualise 2-3 axes at once). Additionally, many SNPs will be highly correlated (linked), meaning that PCA can be highly effective at reducing the dimensionality of the dataset while retaining most of the variance.

To perform PCA, we want to find a rotation and scaling of the data to find new set of orthogonal axes that maximuse variation within the data. We can do this easily on a genotype table, e.g. in R, SNPRelate or other tools, however, recall that for low-coverage data, we potentially have a large uncertainty in genotype calls and consequently in variant calls, in addition to a large amount of missing data. The tool that we will use, PCAngsd performs PCA by first estimating individual allele frequencies (in an iterative approach), and then using the estimated individual allele frequencies as prior information for unobserved genotypes (in the data) to estimate a genetic covariance matrix [(Meisner & Albrechtsen 2018)](https://academic.oup.com/genetics/article/210/2/719/5931101).

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Pcangsd_admix.gif" width="400"> 

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Pcangsd_pca.png" width="400"> 


## Step 1\. PCAngsd input

Let's first open a Docker container in a new tab (running ANGSD):

	docker run --name angsd --mount source=myvol,target=/data -w /data/ -it --rm zjnolen/angsd

Here, we first assign the reference sequence fasta file to a variable, so that it is easy to refer to in the next steps (we can expand a variable *i* via ${*i*}.

	cd /data/Workshop/Data
	REF=/data/Workshop/Data/assembly_homozygous.fa
	metadata_dir=/data/Workshop/Metadata

PCAngsd takes as input genotype likelihoods in beagle format, which we generated in the step before using the `-doGLF 2`option.

	angsd -GL 2 -out GL_75inds -ref ${REF} -nThreads 4 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -only_proper_pairs 1 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -bam ${metadata_dir}/samples_5inds.list -rf ${metadata_dir}/scaffolds.list

## Step 2\. Run PCA via PCAngsd

Let's open a Docker container in a new tab (running PCAngsd) and first look at PCAngsd's options:
	
	docker run --name pcangsd --mount source=myvol,target=/data -w /data/ -it --rm didillysquat/pcangsd
	pcangsd.py -h

Then run the PCA

	cd /data/Workshop/Data
	prefix="GL_75inds"
	pcangsd.py -beagle ${prefix}.beagle.gz -threads 2 -o ${prefix}.pcangsd

We can then plot the results (PC1 vs PC2) in R as follows (don't forget to upload pop metadata files for PCA, admixture, FST to your github!);
	R code

## Step 3\. Plot PCA results

Let's plot the results of the PCA. Let's first copy the output (GL_75inds.pcangsd.cov) from the Docker container (volume) to your local computer (here we create a temporary container (named temp) with our named volume mounted). Let's also download the population metada and some plotting code while we are it.

	docker run --name temp --mount source=myvol,target=/data -w /data ubuntu
	docker cp temp:/data/Workshop/Data/GL_75inds.pcangsd.cov ./Desktop/
	docker cp temp:/data/Workshop/Metadata ./Desktop/
	docker stop temp

The output file together with some plotting scripts should now have been downloaded to your Desktop (if not, please check the paths in the code above).

Open Plot_PCA.R in RStudio and run the code. It will be necessary to change the path on line 10 (after "setwd") of the code to the path of your Desktop (or the path where you downloaded to). If you have time, see if you can follow some of the code.

<br>

<details>

<summary>Click here to expand Plot_PCA.R code</summary>

<br>

	# Import necessary modules
	library(methods)
	library(ggplot2)
	library(plotly)
	library(RcppCNPy)
	library(scales)
	library(ggpubr)

	# Set working directory (change this to your own)
	setwd("/Users/hirzi/Desktop/")

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

</details>

<br>

<details>

<summary>Click here to expand and see PCA result</summary>

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/PCA_result.png" width="800"> 

</details>

<br>

How do we interpet the results? How much variance is explained by the first two principle components (PCs)? 

In certain cases, geographic clines, admixture, bottlenecks and complex demography in general can complicate interpretations of PCA. The following papers provide some good pointers on how to deduce more complex patterns from PCA: [Novembre & Stephens 2008](https://www.nature.com/articles/ng.139), [François et al. 2010](https://academic.oup.com/mbe/article/27/6/1257/1109324), and [Gompert & Buerkle 2016](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12380)

<br>

# Admixture and ancestry analysis

Another way to investigate population structure is by looking at admixture (or ancestry) proportions in populations. This can be achieved through a model-based clustering method where we assume ***K*** populations (***K*** may be unknown), each of which is characterised by a set of allele frequencies. Individuals are then assigned probabilistically to 1 - ***K*** populations (where assignment to > 1 population indicates population admixture). Here, we will use a genotype likelihood implementation of admixture inference using PCAngsd. 

## Step 1\. Infer admixture proportions via PCAngsd. 

Similar to before, PCAngsd takes as input genotype likelihoods in beagle format, which we generated before in the previous section. To estimate admixture proportions in PCangsd, we add the -admix argument. When inferring admixture, it is always advisable to do so over different values of ***K***, so we will iterate the command in a *for* loop. Here, we do so via the -e argument which defines the number of eigenvalues, rather than via -admix_K (as the latter is not recommended). We then define an alpha (sparseness regularisation) parameter, which can be specified manually (-admix_alpha) or here, automatically (-admix_auto) specifying only a soft upper bound.

Let's go back to the terminal running our PCAngsd container and run the following command:

	cd /data/Workshop/Data
	for k in $(seq 1 2); do
		pcangsd.py -beagle ${prefix}.beagle.gz -threads 2 -e ${k} -admix -admix_auto 10000 -o ${prefix}.admix.pcangsd.K$((${k}+1))
	done

## Step 2\. Plot admixture results

Let's plot the results from the admixture analysis. Let's first copy the output (GL_75inds.admix.pcangsd.K3.admix.2.Q, GL_75inds.admix.pcangsd.K3.admix.3.Q) from the Docker container (volume) to your local computer.

	docker run --name temp --mount source=myvol,target=/data -w /data ubuntu
	docker cp temp:/data/Workshop/Data/GL_75inds.admix.pcangsd.K2.admix.2.Q ./Desktop/
	docker cp temp:/data/Workshop/Data/GL_75inds.admix.pcangsd.K3.admix.3.Q ./Desktop/
	docker stop temp

Then, open Plot_admix.R in RStudio (this will be in the same folder as PCA.R used in the previous section). It will be necessary to change the path on line 8 (after "setwd") of the code to the path of your Desktop (or the path where you downloaded to). Run the code for ***K***= 2 and 3 (line 5).

<br>

<details>

<summary>Click here to expand Plot_admix.R code</summary>

<br>

	# Import necessary modules
	library(RcppCNPy)

	# Define number of Ks
	k <- 3

	# Set working directory (change this to your own) and read in data
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

</details>

<br>

<details>

<summary>Click here to expand and see admixture results</summary>

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Admixture_result_K2.png" width="800"> 

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Admixture_result_K3.png" width="800"> 

</details>

<br>

How do we interpet the results? Is there a lot of admixture between populations? If so, between which populations?

Note: it is important to note that inferences of admixture can easily be misinterpreted, e.g. if the dataset cannot biologically be delimited into discrete ***K*** populations (as in the case of a continous geographic cline), [if the dataset is imbalanced](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12512) (with regard to sample sizes of each cluster), or under complex demography (e.g. [recent bottleneck](https://www.nature.com/articles/s41467-018-05257-7)).

<br>

# Site frequency spectrum and summary statistics

To interogate the genealogy of a set of samples, population geneticists typically rely on summary statistics that contain information of the underlying genealogical tree of the data.  Among the most informative (and commonly used) statistics for this is the site-frequency spectrum (SFS). To consider the relation between genealogy (which is not directly observable) and the SFS, consider a genealogical tree.

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/genealogies_demos.png" width="800"> 

Mutations can occur anywhere on the genealogical tree, and we can assume they appear randomly at a relatively fixed rate; hence, mutation events are proportional to branch length. Mutations on terminal (or external) branches are denominated as singletons, since they are unique (private) to one branch (population). Mutations on internal branches are classified e.g. as doubletons, tripletons, etc. depending on how many terminal branches (populations) the mutation is present. The SFS (of a population) is simply the distribution of these frequency classes (singletons, doubletons, tripletons, etc.) in the population.

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/sfs_demos.png" width="800"> 

From the above two figures, we can see that different demographic processes, e.g. constant size, expansion and decline, are expected to effect the genealogy, and hence the SFS, of a population in particular ways. These expectations can be directly derived from coalecent theory (covered in the workshop yesterday, see Wakely's 2009 book "Coalescent theory" for a good overview).

Because of this information held in the SFS, many diversity statistics (e.g. nucleotide diversity, Watterson's theta) and neutrality statistics (e.g. Tajima's D, Fay & Wu's H, Zeng's E) are based on functions (statistical summaries) of the SFS. These summary statistics are relatively easy to compute (easily calculated from the SFS and do not require phasing of genotypes into haplotypes) and can be quite effective in detecting selection on intermediate to long evolutionary timescales. It is important to note that because both selection and demography can affect genealogies (and hence the SFS) in similar ways, their signals can sometimes be confounded. For more recent selection, haplotype-based selection inference methods (such as those based on extended haplotye homozygosity (EHH) and derivatives, see [selscan](https://github.com/szpiech/selscan) for a modern implementation) are more appropriate.

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/sfs_selection.png" width="800">

(*from Nielsen 2005*)

## One population

Here, we will ANGSD to calculate the SFS for single populations, from which we will estimate various thetas and neutrality statistics.

First, let's go back to the terminal running our ANGSD container (or load it again via: docker run --name angsd --mount source=myvol,target=/data -w /data/ -it --rm zjnolen/angsd)

Let's define some path and file variables:
	
	cd /data/Workshop/Data
	metadata_dir=/data/Workshop/Metadata
	REF=/data/Workshop/Data/assembly_homozygous.fa
	out_dir=/data/Workshop/Data/1pop_sumstats
	mkdir ${out_dir}

Here, we will calculate the SFS and summary statistics for three populations comprising 15 individuals each. *Note, in case the below analyses runs slow, consider increasing the number of threads via the -P argument in the angsd command (limited by number of CPU cores/threads on computer).*

	for i in $(seq 1 3); do
		# assign population list
		pop=`sed -n ${i}p < ${metadata_dir}/pop15inds.list`
		#pop=`sed -n ${i}p < ${metadata_dir}/pop5inds.list`
		# assign individual BAM files that constitute population 
		pop_bamlist=${metadata_dir}/Population_lists_15inds/${pop}_15inds.list
		#pop_bamlist=${metadata_dir}/Population_lists_5inds/${pop}.list
		# First, we calculate the site allele frequency likelihoods (SAFs).
		angsd -b ${pop_bamlist} -doSaf 1 -anc ${REF} -ref ${REF} -GL 2 -P 4 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -only_proper_pairs 1 -doCounts 1 -setMinDepth 10 -setMaxDepth 75 -minInd 3 -out ${out_dir}/${pop}
		# Then, we obtain the maximum likelihood estimate of the SFS (here for the folded spectrum since we do not have information of ancestral states)
		realSFS ${out_dir}/${pop}.saf.idx -P 2 -fold 1 > ${out_dir}/${pop}.sfs
		# We calculate the thetas for each site
		realSFS saf2theta ${out_dir}/${pop}.saf.idx -sfs ${out_dir}/${pop}.sfs -outname ${out_dir}/${pop}
		# Anf finally estimate various thetas and neutrality statistics
		thetaStat do_stat ${out_dir}/${pop}.thetas.idx
	done

Let's have a look at the results. 

First, let's look at the SFS (here e.g. for the population "Val_da_la_Stura")

	cd ${out_dir}
	cat Val_da_la_Stura.sfs

Let's plot this in R.

	barplot(sfs[-c(1,length(sfs))]) # where sfs is the vector of frequencies

How does the SFS look like? Is it what we expect?

Let's also have a look at the per-site thetas (again here for the population "Val_da_la_Stura"; but also do so for the other populations). The file can be very big so let's open it with *less* (which we will first need to install)

	apt-get install less
	thetaStat print Val_da_la_Stura.thetas.idx | less -S

Finally, we can have a look at various neutrality statistics for the different scaffold

	cat Val_da_la_Stura.thetas.idx.pestPG

For detecting genetic loci under selection, we can plot the various thetas and neutrality statistics along a sliding window to find loci whose statistics are distinct to that of the rest of the genome (for help with interpreting some the output neutrality statistics, refer to [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1667063) and [here](https://academic.oup.com/genetics/article/207/1/229/5930677)). This typically requires generating a null distribution for the summary statistic (either from the empirical genome-wide distribution or simulated under an appropriate demographic model, as demonstrated in yesterday's workshop) to obtain a measure of statistical significance. 

We can also use the SFS to more explicitly infer the demography of the population, e.g. by leveraging the expected waiting times between coalescent events to assess the variation over time in the effective population size affects e.g. via methods like Bayesian skyline or [stairway](https://github.com/xiaoming-liu/stairway-plot-v2) plots, or by comparing simulated SFS (generated under particular demographic models) against the empirical (observed) SFS (e.g. via [dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master), [moments](https://bitbucket.org/simongravel/moments/src/main), [momi](https://github.com/popgenmethods/momi2), [fastsimcoal](http://cmpg.unibe.ch/software/fastsimcoal27)). 

Note: We could also calculate the SFS for all populations in our dataset in Docker, but note that for the other 12 populations, we only have 5 individuals per population. This means that estimations of allele frequencies and summary statistics will be much more coarse (as we can detect fewer frequency classes). You can try to calculate the SFS and summary statstics based on 5 individuals per population to see this effect [OPTIONAL]. *The trade-off of having more populations but with less accurate data per population vs. less populations with more accurate data per population is an important consideration when designing population genetic studies.*

## Two populations

A useful statistic to calculate between pairs of populations is their population genetic differentiation (***FST***). Here, we will use ANGSD to calculate the ***FST*** between our population pairs. To do this, ANGSD first calculates the two population or joint (2D) SFS, which it uses as a prior (together with the site allele frequency likelihoods of the single populations, calculated in the previous step) to calculate the ***FST***.

To calculate ***FST*** between all our population pairs, let's first make a list of all population pairs.

	cd ${metadata_dir}
	chmod u+x make_list_pairwise.sh
	./make_list_pairwise.sh

We can then calculate the ***FST*** for each population pair by looping over each line of the pairwise list. *Note, this may take a few minutes to complete. In case it takes too long, consider increasing the number of threads via the -P argument in the angsd command (limited by number of CPU cores/threads on computer).*

	cd /data/Workshop/Data
	for i in $(seq 1 3); do
		# read in population pairs list
		pop=`sed -n ${i}p < ${metadata_dir}/pop_name_pairs`
		# Note, set allows you to define the elements of your list as variables, according to their order
		set -- $pop
		# Calculate the 2D SFS prior
		realSFS ${out_dir}/${1}.saf.idx ${out_dir}/${2}.saf.idx -P 4 -fold 1 > ${out_dir}/${1}.${2}.sfs
		# Prepare FSTs for easy window analysis. The option -whichFst 1 is preferable for small sample sizes.
		realSFS fst index ${out_dir}/${1}.saf.idx ${out_dir}/${2}.saf.idx -sfs ${out_dir}/${1}.${2}.sfs -fstout ${out_dir}/${1}.${2}.stats -whichFst 1
		# Get the global estimate
		FST=`realSFS fst stats ${out_dir}/${1}.${2}.stats.fst.idx`
		printf "${1}\t${2}\t${FST}\n" >> ${out_dir}/FST_summary.txt
		cat ${out_dir}/FST_summary.txt
		# Get sliding window estimates
		realSFS fst stats2 ${out_dir}/${1}.${2}.stats.fst.idx -win 10000 -step 1000 > ${1}.${2}.stats.slidingwindow
	done

Let's have a look at the ***FST*** results.

	cat FST_summary.txt
	Kapetanovo_Jezero.Legn_Mar_Scuol.stats.slidingwindow | less -S
	Kapetanovo_Jezero.Val_da_la_Stura.stats.slidingwindow  | less -S
	Kapetanovo_Jezero.Legn_Mar_Scuol.stats.slidingwindow  | less -S

Which population pairs are genetically closest to each other? Which are furthest apart? A common way to detect selection leveraging population pairs is to calculate ***FST*** in a sliding window across the genome to find loci whose ***FST*** is significantly higher (indicative of positive selection in one population but not the other) or lower (potentially indicative of balancing selection) than the genome-wide average. As in the case of single-population outlier stastistics, this typically requires generating a null distribution for ***FST*** to obtain a measure of statistical significance.

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Manhattan_plot.png" width="600">

In addition to using the 2D SFS as a prior to calculate ***FST***, we can use it directly to infer the demography of the population pair, by comparing simulated SFS (generated under particular demographic models) against the empirical (observed) SFS (e.g. via [dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master), [moments](https://bitbucket.org/simongravel/moments/src/main), [momi](https://github.com/popgenmethods/momi2), [fastsimcoal](http://cmpg.unibe.ch/software/fastsimcoal27)). Here, migration rates and time of population divergence, in addition to effective population sizes, are typical demographic parameters of interest.

## Three populations

When we have three populations, we can measure the branch lengths between the three populations (employing one as an outgroup) to detect extreme allele frequency change (in a locus) in one population (relative to the other non-outgroup population). This statistic, called the population branch statistic (***PBS***) represents another way to infer selection.

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/PBS.png" width="300">

(*EPAS 1 gene; Yi et al. 2010*)

We can use ANGSD to calculate the PBS via extending the previous ***FST*** command to three populations.

	cd /data/Workshop/Data/1pop_sumstats
	realSFS fst index Kapetanovo_Jezero.saf.idx Legn_Mar_Scuol.saf.idx Val_da_la_Stura.saf.idx -sfs Kapetanovo_Jezero.Legn_Mar_Scuol.sfs -sfs Kapetanovo_Jezero.Val_da_la_Stura.sfs -sfs Legn_Mar_Scuol.Val_da_la_Stura.sfs -fstout 3pops_KJ_LMS_VDLS.pbs -whichFst 1
	# Get global estimates
	PBS=`realSFS fst stats ${out_dir}/3pops_KJ_LMS_VDLS.pbs.fst.idx`
	printf "Kapetanovo_Jezero\tLegn_Mar_Scuol\tVal_da_la_Stura\t${PBS}\n" > ${out_dir}/PBS_summary.txt
	# Get sliding window estimates
	realSFS fst stats2 3pops_KJ_LMS_VDLS.pbs.fst.idx -win 10000 -step 1000 > 3pops_KJ_LMS_VDLS.pbs.slidingwindow

Let's have a look at the results.

	cat PBS_summary.txt
	realSFS fst print 3pops_KJ_LMS_VDLS.pbs.fst.idx | less -S

Here, the columns are: region	chr	midPos	Nsites	Fst01	Fst02	Fst12	PBS0	PBS1	PBS2, and populations 0, 1 and 2 refer to Kapetanovo_Jezero, Legn_Mar_Scuol and Val_da_la_Stura respectively. Note that sometimes we may get negative PBS and FST values in the output; these are equivalent to 0.

Separately, we can calculate the three population (3D) SFS as follows (note that calculating the joint SFS for >/=3 populations can be very slow, *so we won't run this*):

	realSFS ${out_dir}/Kapetanovo_Jezero.saf.idx ${out_dir}/Legn_Mar_Scuol.saf.idx ${out_dir}/Val_da_la_Stura.saf.idx -P 8 > ${out_dir}/3pops_KJ_LMS_VDLS.sfs

Similar to above, the 3D-SFS can be used infer the demography of the population trio, by comparing simulated SFS (generated under particular demographic models) against the empirical (observed) SFS (e.g. via [dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master), [moments](https://bitbucket.org/simongravel/moments/src/main), [momi](https://github.com/popgenmethods/momi2), [fastsimcoal](http://cmpg.unibe.ch/software/fastsimcoal27)).

<br>

===================================================================

<br>

**You have now learnt how to infer population structure (including admixture), estimate various diversity metrics and perform tests for selection on low-coverage WGS data - well done! :-)**

<br>
