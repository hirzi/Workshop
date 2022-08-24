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
      - [For single populations](#for-single-populations)
      - [For two populations; population genetic differentiation](#for-two-populations\;-population-genetic-differentiation)
      - [For three populations](#for-three-populations)

<br>

# Introduction

## Why go low-coverage?

**Experimental design** All scientific lines of inquiry start with a question. From this question, we (the researcher) try to come up with an experimental design to best address said question. If money and time were no object, we may e.g. plan for an experiment with 10 treatments, 10 replicates per treatment, and 1000 samples per replicate. However, in the real world, the experimental design is not determined purely by how best to address the biological question at hand, but also by *cost, time and technical feasibility*.

Typically, we are faced with the following trade-off; of having either **1)** more samples at the expense of data per sample or **2)** less samples but with more data per sample.

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Experimental_design.png" width="650"> 

With respect to sequencing and genetic data, if we start with the perfect or complete representation of a unit data as the whole genome sequenced at high coverage (e.g. 50x), less data can imply one of two things: **1)** sequencing a reduced or sub-representation of the genome, i.e. using genetic markers like microsatellites and sets of single nucleotide polymorphisms (SNPs) or **2)** sequencing the whole genome but at low coverage (depth). I.e. a trade-off of breadth vs depth. To give a concrete example, imagine you had enough money to sequence 1 million reads, and that this is sufficient to sequence your whole genome at 2x or 10% of your genome at 20X, *which would you choose*? 

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/breadth_vs_depth.png" width="650"> 

The answer can be difficult, and lies in weighing the respective advantages and disadvantages of these two alternatives, in the context of the biological question at hand. Briefly, both approaches have their caveats. For the former (i.e. subset of genetic markers), we **1)** assume the sub-selection of the genome to be representative of the whole genome, **2)** are prone to ascertainment bias, and **3)** lack data in unsequenced parts of the genome which prohibits detection of new, potentially interesting genetic variants. For the latter (low-coverage sequencing), our certainty in the genotype call (i.e. whether something is **A**, **C**, **T** or **G**) is much lower, due to the fact that we are reading each position fewer times, and hence are prone to more sequencing errors in our genotype calls. Further, the power to detect rare variants is strongly diminished.

That said, in the last years (for both cases), there has been notable advances in alleviating these respective downsides. For the former, modern sequencing techniques provide more markers (at increasingly cheaper prices), while for the latter, we now have methods that explicitly accommodate the uncertainty in genotype calls in analyses (i.e., rather than having to explicitly call genotypes, we consider their genotype likelihoods instead, aka the probability of the data given a specific genotype).

Low-coverage methods also lend themselves well to the sequencing and analysis of ancient DNA, where high-coverage, high-quality DNA sequences may not be attainable.


### Summary
Through working with genotype likelihoods rather than relying on discrete (lossy) genotype calls, low-coverage methods are useful for they allow for the statistical propagation of uncertainty from raw sequencing data to downstream analysis. The effect may be huge for data at low-coverage or minimal for data at high-coverage (where analytical results will tend to converge between low-coverage methods and classical genotype-call based methods). In addition, working directly with genotype likelihoods generally involve fewer (potentially lossy) processing steps, e.g. genotype calling and various filtering.

<br>

# Workshop - initial preparation

In this session you will learn how to use low-coverage whole genome data to infer population structure and admixture (ancestry), via two methods:

  - **Principal Components Analysis (PCA)**
  - **Admixture analysis**

Population genetic analyses of NGS data is typically run on large linux-based computing clusters. For this workshop, since we do not have access to this, we will be running population genetic analyses in **Docker**. 

### What is Docker?
Docker is an open source container based technology that runs in an isolated, self-contained package that can be efficiently distributed and executed in a portable manner across a wide range of computing platforms. Containerization in concept is very similar to virtualization, i.e. a method of isolating an application from the underlying machine operating system. The difference between virtual machines and containers is that containers do not require a full operating system to operate but rather the application and dependencies, means they are much smaller in size (Gharib 2018). 

## Step 1. Make sure you have Docker Desktop installed on your computer
Instructions for intalling Docker can be found [here](https://www.docker.com/get-started/), with OS-specific information (including minimum system requirements) and download links for [Windows](https://docs.docker.com/desktop/install/windows-install/), [Mac](https://docs.docker.com/desktop/install/mac-install/) and [Linux](https://docs.docker.com/desktop/linux/).

## Step 2. Make sure you are familiar with the basics of bash shell scripting
We will be working almost exclusively through the command-line in docker, so if you have not used shell scripting before or are a little rusty with it, it may be helpful to have a look at a tutorial like [this](https://linuxconfig.org/bash-scripting-tutorial-for-beginners) or a cheat sheet like [this](https://bioinformaticsworkbook.org/Appendix/Unix/UnixCheatSheet.html#gsc.tab=0) before proceeding with this workshop.

<br>

# Workshop - data

We will be working with low-coverage whole-genome sequencing (WGS) data (average 2x) of a plant species *Dianthus sylvetris* (Wood pink). This is a perennial plant species that grows throughout the mountain ranges of Europe (e.g. Alps, Apennines & Dinarides). Given the European mountain ranges experienced repeated bouts expansion and recession of glacial ice sheets during the last ca. 2 million years (the Quaternary glaciations or "ice ages"), this species likely experienced a complex and dynamic demographic history. Additionally, this species inhabits a large elevational range (0 - 2500 meters), with the consequence that contemporary populations exhibit a remarkable degree of local adaptation in phenotypic and life history traits. Thee genetic bases of these observed adaptations, however, remain unknown.

Add map!

The populations and data that we will use in this workshop represents a small subset from a larger study that covered the geographic and ecological range of the species (https://www.biorxiv.org/content/10.1101/2022.06.07.495159v1). Our data is in BAM format (i.e. mapped sequencing data) and span a ~2MB region from four scaffolds selected at random across the genome. (may want to extract an interesting region for FST/PBS analysis, 15 inds pops. alternatively, consider taking data directly from Simons's data and do FST scans. Oruse selscan on vcf). We will use this data to interogate this species' population structure via principle component analysis (PCA) and admixtre/ancestry analysis.

<br>

# Workshop - programs we will use
For this practical, we will be using:

[SAMtools](https://www.htslib.org/)

[ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next
Generation Sequencing Data)

[PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd)

[ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) (think if want to include!)

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
We have now "pulled" the necessary programs into Docker. We now need the data. First we'll need to create a volume in Docker to store downloaded and generated data. Here, we will run an Ubuntu container with a named volume via the --mount argument (i.e. we will name this created volume "myvol" and define the associated container path as "/data". The -w allows the command to be executed inside the defined working directory. Before we download data, we will first need to install some other programs within the running container, namely git. 

	docker run --name base --mount source=myvol,target=/data -w /data -it --rm ubuntu
	apt update
	apt install git

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

To perform PCA, we want to find a rotation and scaling of the data to find new set of orthogonal axes that maximuse variation within the data. We can do this easily on a genotype table, e.g. in R, SNPRelate or other tools, however, recall that for low-coverage data, we potentially have a large uncertainty in genotype calls and consequently in variant calls. The tool that we will use, PCAngsd performs PCA by first estimating individual allele frequencies (in an iterative approach), and then using the estimated individual allele frequencies as prior information for unobserved genotypes (in the data) to estimate a genetic covariance matrix [(Meisner & Albrechtsen 2018)](https://academic.oup.com/genetics/article/210/2/719/5931101).

<br>

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Pcangsd_admix.gif" width="400"> 

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/Pcangsd_pca.png" width="400"> 


## Step 1\. PCAngsd input

Here, we first assign the reference sequence fasta file to a variable, so that it is easy to refer to in the next steps (we can expand a variable *i* via ${*i*}.

	REF=/data/Workshop/Data/assembly_homozygous.fa

PCAngsd takes as input genotype likelihoods in beagle format, which we generated in the step before using the `-doGLF 2`option.

	angsd -GL 2 -out GL_95inds -ref ${REF} -nThreads 4 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -only_proper_pairs 1 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -bam samples.list -rf scaffolds.list
	angsd -GL 2 -out GL_75inds -ref ${REF} -nThreads 4 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -only_proper_pairs 1 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -bam samples_5inds.list -rf scaffolds.list

## Step 2\. Run PCA via PCAngsd. 

Let's first look at PCAngsd's options:
	
	pcangsd.py -h

Then run the PCA (adjust the number of threads accordingly)

	prefix="GL_95inds"
	prefix="GL_75inds"
	pcangsd.py -beagle ${prefix}.beagle.gz -threads 2 -o ${prefix}.pcangsd

We can then plot the results (PC1 vs PC2) in R as follows (don't forget to upload pop metadata files for PCA, admixture, FST to your github!);
	R code

## Step 3\. Copy file to/from the Docker container (volume) to local computer. Alternatively, let's write R script and upload on Github repo. Then easy to run as Rscript, output as html. Then copy to local computer.

Copy files to/from container (volume) to local filesystem (here we create a temporary container (named temp) with our named volume mounted)

	docker run --name temp --mount source=myvol3,target=/data -w /data ubuntu
	docker cp temp:/data/Workshop/Data/GL_95inds.pcangsd.cov ./Desktop/Test/
	#remove/detach docker contained after copying

## Step 4\. Plot PCA results

Let's add the plotly PCA results (embed html)!  Let's add map to githubs totorial.

How do we interpet the results?

We find three distinct clusters. How much variance is explained by the first two principle components (PCs)? 

Here, the PCA is relatively straightforward to interpret, however, in many cases, geographic clines, admixture, bottlenecks and complex demography in general can make interpretations of PCA more difficult to interpret. The following papers provide some good pointers on how to deduce more complex patterns from PCA: [(Novembre & Stephens 2008)] (https://www.nature.com/articles/ng.139), [(François et al. 2010)] (https://academic.oup.com/mbe/article/27/6/1257/1109324), and [(Gompert & Buerkle 2016)] (https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12380)

# Admixture and ancestry analysis

# Add refs everywhere!

Another way to investigate population structure is by looking at admixture (or ancestry) proportions in populations. This can be achieved through a model-based clustering method where we assume ***K*** populations (***K*** may be unknown), each of which is characterised by a set of allele frequencies. Individuals are then assigned probabilistically to 1 - ***K*** populations (where assignment to > 1 population indicates population admixture). Here, we will use a genotype likelihood implementation of admixture inference using PCAngsd. 


## Step 1\. Infer admixture proportions via PCAngsd. 

Similar to before, PCAngsd takes as input genotype likelihoods in beagle format, which we generated before in the previous section. To estimate admixture proportions in PCangsd, we add the -admix argument. When inferring admixture, it is always advisable to do so over different values of ***K***, so we will iterate the command in a *for* loop. Here, we do so via the -e argument which defines the number of eigenvalues, rather than via -admix_K (as the latter is not recommended). We then define an alpha (sparseness regularisation) parameter, which can be specified manually (-admix_alpha) or here, automatically (-admix_auto) specifying only a soft upper bound.

	for k in $(seq 1 2); do
		pcangsd.py -beagle ${prefix}.beagle.gz -threads 2 -e ${k} -admix -admix_auto 10000 -o ${prefix}.admix.pcangsd.K$((${k}+1))
	done

## Step 2\. Copy file to/from the Docker container (volume) to local computer. Alternatively, let's write R script and upload on Github repo. Then easy to run as Rscript, output as html. Then copy to local computer.

Either write R code, or upload R code on Github repo to allow easy execution of Rscript in Docker

Copy files to/from container (volume) to local filesystem (here we create a temporary container (named temp) with our named volume mounted)

	docker run --name temp --mount source=myvol3,target=/data -w /data ubuntu
	docker cp temp:/data/Workshop/Data/admix_results ./Desktop/Test/
	#remove/detach docker contained after copying

## Step 3\. Plot admixture results

Let's add the plotly/interactive admixture results (embed html)!

Similar to our PCA analysis, we find three distinct clusters, with minimal admixture between clusters. It is important to note that inferences of admixture can easily be misinterpreted, e.g. if the dataset cannot biologically be delimited into discrete K populations (as in the case of continous geographic clines), if the dataset is imbalanced (with regard to sample sizes of each cluster), or under complex demography (e.g. [(bottleneck)] (https://www.nature.com/articles/s41467-018-05257-7)).

[(Puechmaille 2016)] (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12512)
[(Gilbert 2016)] (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12521)
[(Meirmans 2015)] (https://onlinelibrary.wiley.com/doi/10.1111/mec.13243)

# Site frequency spectrum and summary statistics

To interogate the genealogy of a set of samples, population geneticists typically rely on summary statistics that contain information of the underlying genealogical tree of the data.  Among the most informative (and commonly used) statistics for this is the site-frequency spectrum (SFS). To consider the relation between genealogy (which is not directly observable) and the SFS, consider a genealogical tree.

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/genealogies_demos.png" width="650"> 

Mutations can occur anywhere on the genealogical tree, and we can assume they do so randomly at a relatively fixed rate; hence, mutation events are proportional to branch length. Mutations on terminal (or external) branches are denominated as singletons, since they are unique (private) to one branch (population). Mutations on internal branches are classified e.g. as doubletons, tripletons, etc. depending on how many terminal branches (populations) the mutation is present. A histogram of these frequency classes (singletons, doubletons, tripletons, etc.) then defines the SFS.

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/sfs_demos.png" width="650"> 

From the above two figures, we can see that different demographic processes, e.g. constant size, expansion and decline, are expected to effect the genealogy, and hence the SFS, of a population in particular ways. These expectations can be directly derived from coalecent theory (covered in the workshop yesterday, see Wakely's 2009 book "Coalescent theory" for a good overview).

Because of this information held in the SFS, many diversity statistics (e.g. nucleotide diversity, Watterson's theta) and neutrality statistics (e.g. Tajima's D, Fay & Wu's H, Zeng's E) are based on functions (statistical summaries) of the SFS. These summary statistics are relatively easy to compute (easily calculated from the SFS and do not require phasing of genotypes into haplotypes) and can be quite effective in detecting selection on intermediate to long evolutionary timescales (see [(here)] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1667063) and [(here)] (https://academic.oup.com/genetics/article/207/1/229/5930677) for good references). It is important to note that because both selection and demography can affect genealogies (and hence the SFS) in similar ways, their signals can sometimes be confounded. For more recent selection, haplotype-based selection inference methods (such as those based on extended haplotye homozygosity (EHH) and derivatives, see [(selscan)] (https://github.com/szpiech/selscan) for a modern implementation) are more appropriate.

<img src="https://github.com/hirzi/Workshop/blob/main/Example_figures/sfs_selection.png" width="650">

## Single populations

Here, we will ANGSD to calculate the SFS for single populations, from which we will estimate various thetas and neutrality statistics.

	for i in $(seq 1 15); do
		pop=`sed -n ${i}p < pop5inds.list`
		pop_bamlist='./Population_lists_5inds/'${pop}.5inds
		# First, we calculate the site allele frequency likelihoods (SAFs).
		angsd -b ${pop_bamlist} -doSaf 1 -anc ${REF} -ref ${REF} -GL 2 -P 2 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -only_proper_pairs 1 -doCounts 1 -setMinDepth 10 -setMaxDepth 75 -minInd 3 -out ${OUT}/${pop}
		# Then, we obtain the maximum likelihood estimate of the SFS (here for the folded spectrum since we do not have information of ancestral states)
		realSFS ${OUT}/${pop}.saf.idx -P 2 -fold 1 > ${OUT}/${pop}.sfs
		# We calculate the thetas for each site
		realSFS saf2theta ${OUT}/${pop}.saf.idx -sfs ${OUT}/${pop}.sfs -outname ${OUT}/${pop}
		# Anf finally estimate various thetas and neutrality statistics
		thetaStat do_stat ${OUT}/${pop}.thetas.idx
	done

Let's have a look at the results. For detecting genetic loci under selection, we can plot the various thetas and neutrality statistics along a sliding window to find loci whose statistics are distinct to that of the rest of the genome. This typically requires generating a null distribution for the summary statistic (either from the empirical distribution or simulated under an appropriate demographic model) to obtain a measure of statistical significance. We can also use the single population (1D)-SFS to more explicitly interogate the demography of the population, e.g. by leveraging the expected waiting times between coalescent events to assess the variation over time in the effective population size affects e.g. via stairway plot, or by comparing simulated SFS (generated under particular demographic models) against the empirical SFS (e.g. dadi, moments, momi2, fastsimcoal). 

## Two populations

A common statistic to calculate between pairs of populations is their population genetic differentiation (***FST***). Here, we will use ANGSD to calculate ***FST*** between our population pairs. ANGSD first calculates the two population (2D)-SFS, which it then uses as a prior (joint with the site allele frequency likelihoods of the single populations, calculated in the previous step) to calculate the ***FST***.

To calculate ***FST*** between all our population pairs, we will first make a list of all population pairwise comparisons. We do this this with the make_list_pairwise.sh script in the xxx/Data/xxx folder.

	./make_list_pairwise.sh

We can then calculate ***FST***s for each population pair by looping over each line of the pairwise list.

	for i in $(seq 1 15); do
		pop=`sed -n ${i}p < ${working_dir}/list_pop_name_pairs/pop_name_pairs`
		# Note, set allows you to define the elements of your list as variables, according to their order
		set -- $pop
		# Calculate the 2D SFS prior
		realSFS ${input_dir}/${1}.saf.idx ${input_dir}/${2}.saf.idx -P 2 -fold 1 > ${output_dir}/${1}.${2}.ml
		# Prepare FSTs for easy window analysis. The option -whichFst 1 should be preferable for small sample sizes.
		realSFS fst index ${input_dir}/${1}.saf.idx ${input_dir}/${2}.saf.idx -sfs ${output_dir}/${1}.${2}.ml -fstout ${output_dir}/${1}.${2}.stats -whichFst 1
		# Get the global estimate
		FST=`/cluster/project/gdc/shared/tools/angsd0_933/angsd/misc/realSFS fst stats ${output_dir}/${1}.${2}.stats.fst.idx`
		printf "${1}\t${2}\t${FST}\n" >> ${output_dir}/FST_summary.txt
		# Get sliding window estimates
		realSFS fst stats2 ${output_dir}/${1}.${2}.stats.fst.idx -win 50000 -step 10000 > slidingwindow
	done

Let's have a look at the ***FST*** results. Which population pairs are genetically closest to each other? Which are furthest apart? You can use the included plotting script (xxx.R) to plot a heatmap of the ***FST*** results (PREPARE R PLOTTING SCRIPT FOR FST, EMBED THE RESULTS HERE!!). Similar to the 1D-SFS, the 2D-SFS can be used infer the demography of the population pair, by comparing the simulated SFS (generated under particular demographic models) against the empirical SFS (e.g. dadi, moments, momi2, fastsimcoal). Here, additional demographic parameters, such as/namely migration and population divergence,....


## Three populations

For three populations, we can consider the population branch statistic (***PBS***). We can use ANGSD to calculate the PBS via a similar command to that for ***FST*** (Select 3 populations, i.e. one population each from each lineage. Code is similar to that above for FST (see ANGSD and low-cov tutorial website)_

	realSFS fst index pop1.saf.idx pop2.saf.idx pop3.saf.idx -sfs pop1.pop2.ml -sfs pop1.pop3.ml -sfs pop2.pop3.ml -fstout out.pbs -whichFst 1

Separately, we can calculate the three population (3D)-SFS in ANGSD as follows (not that joint SFS for >3 populations can be very slow):

	realSFS ${input_dir}/${1}.saf.idx ${input_dir}/${2}.saf.idx ${input_dir}/${3}.saf.idx -P 8 > ${output_dir}/${1}.${2}.${3}.sfs

, PBS, 3D-SFS, input for dadi/moments/momi2

Let's have a look at the ***PBS*** results. 

Similar to above, the 3D-SFS can be used infer the demography of the population trio, by comparing the simulated SFS (generated under particular demographic models) against the empirical SFS (e.g. dadi, moments, momi2, fastsimcoal).

<br>

===================================================================
