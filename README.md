# Workshop 4 - Population structure and demography of low-coverage whole genome sequencing data

Why go low?


Experimental design.
All scientific lines of inquiry start with a question. From this question, the researcher comes up with an experimental design that can best address said question. If money and time were no object, the researcher would e.g. plan for an experiment with 10 treatments, 10 replicates per treatment, and 1000 samples per replicate. However, in the real world, the experimental design is not determined purely by how best to address the biological question at hand, but also by cost, time and technical feasibility.

Pic

Generally, one is faced with the following trade-off; of having either 1) more samples at the expense of data per sample or 2) less samples but with more data per sample.

With respect to  genetics and sequencing, if we start with the perfect or complete representation of a unit data as the whole genome sequenced at high coverage (e.g. 50x), less data can imply 1 of two things: 1) sequencing a reduced or sub- representation of the genome, i.e. using genetic markers like microsats and SNPs or 2) sequencing the whole genome but at low coverage. I.e. a trade-off of breadth vs depth. To give a concrete example, imagine you had enough money to sequence 1 million reads, and that this is sufficient to sequence your whole genome at 2x or 10% of your genome at 20X, which would you choose? 

Pic

The answer can be difficult, and lies in weighing the respective advantages and disadvantages of these 2 alternatives, in the context of the biological question at hand.
Briefly, both approaches have their caveats. For the former (i.e. genetic markers), you make the assumption that your sub-selection of the genome is representative of the whole genome, you are prone to ascertainment bias, and you lack data in unsequenced parts of the genome, hence making it inappropriate for e.g. when looking for new genetic variants. For the latter (low-coverage), you’re certainty in the genotype call (whether something is A,C,T,G) is much lower, due to the fact that you’re reading each position fewer times, and hence you’re prone to more sequencing errors in your genotype calls.

That said, in the last years for both cases, there has been notable advances in alleviating these respective downsides. For the former, we’ve steadily been developing techniques which provide more and more markers, while for the latter, we now have methods that explicitly accommodate the uncertainty in genotype calls in our analysis, or put another way, we don’t have to explicitly call genotypes but rather we can consider their genotype likelihoods instead (aka the probability of the data given a specific genotype).

Low-coverage methods also lend themselves well to the sequencing and analysis of ancient DNA, where high-coverage, high-quality DNA sequences may not be easily attainable.

To be concise, it is the statistic propagation of uncertainty from raw sequenving data to downstream analysis (via working with genotyp likelihoods rather than discrete (lossy) genotype calls) that make low-coverage methods useful. The effect may be huge at low coverage or minimal/negligable at high coverage (where results will tend to converge to classical genotype-call based methods). By working directly with genotype likelihoods, less lossy) steps (e.g. genotype calling, various filtering) need to be performed, leading to less loss of potentially informative data.

# Workshop 4 (morning session) goals.

<br> <br>

In this session you will learn how to use low-coverage whole genome data
to do:

  - Principal Components Analysis (PCA)
  - Admixture analysis

# Initial preparation

Population genetic analyses of NGS data is typically run on large linux-based computing clusters. For this workshop, since we don't have access to this, we will be running/performing population genetic analyses in Docker. 

## What is Docker?
Docker is an open source container based technology that runs in an isolated, self-contained package that can be efficiently distributed and executed in a portable manner across a wide range of computing platforms. Containerization in concept is very similar to virtualization, i.e. a method of isolating an application from the underlying machine operating system. The difference between virtual machines and containers is that containers do not require a full operating system to operate but rather the application and dependencies, means they are much smaller in size (Gharib 2018). 

## 1. Make sure you have Docker Desktop installed on your computer (choose installation depending on your computer's OS).
https://www.docker.com/get-started/
https://docs.docker.com/desktop/install/windows-install/
https://docs.docker.com/desktop/install/mac-install/
https://docs.docker.com/desktop/linux/

Note that Docker has the following minimum system requirements:
Windows 10 64-bit: Home or Pro 21H1 (build 19043) or higher, or Enterprise or Education 20H2 (build 19042) or higher.
macOS must be version 10.15 or newer. That is, Catalina, Big Sur, or Monterey.
64-bit processor
4GB system RAM

## 2. Make sure you are familiar with the basics of bash shell scripting.
We will be working almost exclusively through the command line in docker, so if you have not used shell scripting before or are getting rusty on it, it may be helpful to have a look at a tutorial like [this one](https://linuxconfig.org/bash-scripting-tutorial-for-beginners) or a cheat sheet like [this one](https://bioinformaticsworkbook.org/Appendix/Unix/UnixCheatSheet.html#gsc.tab=0) before proceeding to the next step.

# Workshop data

We will be working with low-coverage NGS date (average 2x) of a plant species Dianthus sylvetris (Wood pink). This is a perennial plant species that grows throughout the mountain ranges of Europe (Alps, Apennines, Dinarides). Because the European mountain ranges experienced repeated bouts expansion and recession of glacial ice sheets during the last ca. 2 million years (the Quaternary glaciations or "ice ages"), this species likely experienced a complex demographic history. Additionally, the species inhabits a large elevational range (0-2500m), with the consequence that contemporary populations exhibit a remarkable degree of local adaptation in phenotypic and life history traits. We are interested to investigate the genetic bases of these observed adaptations.

The populations and data that we will use in this workshop represents a small subset from a larger study that covered the geograpihc and ecologcal range of the species (https://www.biorxiv.org/content/10.1101/2022.06.07.495159v1). Our NGS data is in BAM format (i.e. mapped sequencing data) and span a ~2MB region from 4 scaffolds selected at random across genome. (may want to extract an interesting region for FST/PBS analysis, 15 inds pops. alternatively, consider taking data directly from Simons's data and do FST scans. Oruse selscan on vcf). We will use this data to interogate this species' population structure via principle component analysis (PCA) and admixtre/ancestry analysis.


# Programs
For this practical, we will be using:
[SAMtools] (https://www.htslib.org/),
[ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next
Generation Sequencing Data),
[ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and
[PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd).
We will run the programs in a Docker container so there is no need to install these programs locally on your computer.

# Instructions - preparing our Docker containers

## Step 1. Open Docker
First, open Docker Desktop on your computer. Then, open Windows Terminal (Powershell), Mac Terminal or Linux terminal. If Docker has ben succesfully installed, you should see the Docker help menu after typing "docker" in the terminal. Docker itself has a whole syntax for usage ([Docker manual](https://docs.docker.com/engine/reference/commandline/docker/), but today we will focus more on the bash-based command-lines as that is what is mainly used to work with NGS data (again here Docker is only used as container to run NGS programs without installing them locally).

## Step 2. Running Ubuntu on Docker
We first need a Linux distribution (e.g. Ubuntu) running in our Docker container. Use "docker pull" to pull a Docker image or repository from a registry. 

	docker pull ubuntu

This will pull the latest version of Ubuntu to Docker.

To run Ubuntu, we use "docker run"

	docker run -it --rm ubuntu

The -it instructs Docker to allocate a pseudo-TTY connected to the container’s stdin; creating an interactive bash shell in the container. I.e. this allows us to run Ubuntu interactively in Docker. --rm automatically remove the container when it exits. To exit the bash interactive shell, enter: exit 13.

## Step 3. Pull Docker images

Let's now "pull" Docker images for the other programs we will be using.

	docker pull biocontainers/samtools:v1.9-4-deb_cv1
	docker pull zjnolen/angsd
	docker pull didillysquat/pcangsd

## Step 4. Creating and mounting a volume
We have now "pulled" the necessary programs into Docker. We now need the data. First we'll need to create a volume in Docker to store downloaded and generated data. Here, we will run an Ubuntu container with a named volume via the --mount argument (i.e. we will name this created volume "myvol" and define the associated container path as "/data". The -w allows the command to be executed inside the defined working directory. Before we download data, we will first need to install some other programs within the running container, namely git. 

	docker run --name base --mount source=myvol,target=/data -w /data -it --rm ubuntu
	apt update
	apt install git

While we're at it, let's install a text editor too (vim nano).

	apt-get install vim nano

## Step 5. Download data
Then we download data. All the BAM files as well as population metadata are deposited in the github: https://github.com/hirzi/Workshop (make public!). Let's pull this directory to our Docker container and our named volume.

	git clone https://github.com/hirzi/Workshop.git

Once downloaded, let's navigate to the /data directory containing the samples. Once there and before we continue, let's download the reference sequence. The reference is deposited in a separate repository (because it's large file).

	wget https://www.dropbox.com/s/1nggjjhrcjseuwx/assembly_homozygous.fa?dl=1

Let's rename this reference sequence file.
	
	mv 'assembly_homozygous.fa?dl=1' assembly_homozygous.fa

## Step 6. Getting a hand of the command-line in Docker

Now that we have all the data downloaded, let's try a few bash commands to see what we have.

	cd /data/Workshop/Data/
	ls
	ls -lrth *bam

How many BAM files (i.e. samples) are there? (hint: ls \*bam | wc -l). In total we have 95 samples (13 populatiosn with 5 individuals each and 2 populations with 15 individuals each).

Now, let's check directory permissions. Do we have permission to write to the /data/Workshop/Data/ directory? You can check this with the command "ls -l". If not, let's change permissions so that you can write to this directory.
chmod 777 /data/Workshop/Data

# Step 7. Index files
We now have the data stored in a named volume. Even when we close the Docker container, the volume is maintained. We can always re-access the volume by mounting it it when running a new container. Before we perform any analysis with BAM and fasta files (or really any large sequencing file), we need to index them so that programs can parse through them efficiently. To do this, let's open another terminal tab, in which we will run Samtools. Note the named volume that we mount to.

	docker run --name samtools --mount source=myvol,target=/data -w /data/ -it --rm biocontainers/samtools:v1.9-4-deb_cv1
	
Index the fasta file (i.e. the reference sequence):

	samtools faidx assembly_homozygous.fa

Index all the BAM files (here we'll do this in a for loop):

	for i in $(ls *bam); do samtools index ${i}; done


!!!!!!!!!!!!!!! Make sample list (let's actually have this in the Github, in case htere are difference with file ordering)
!!!!!!!!!!!!!!!#ls *bam > samples.list



# Principle component analysis

What is principle component analysis (PCA)? PCA is a method that reduces the dimensionality of a dataset to emphasize variation and bring out strong (and more easily) interpretable) patterns in a dataset. It does this by transforming a large set of (potentially correlated) variables into a smaller set of uncorrelated variables while minimising information loss.

We can visualise how a PCA works through this helpful, interactive link:

https://setosa.io/ev/principal-component-analysis/ - link or embed!


For NGS data, each variant site (e.g. single nucleotide polymorphism (SNP)) represents one dimesion (axis) of variation. So if we have 1 million SNPs, we have 1 million dimensions (axes) of variation. Obviously, this would be extremely difficult to visualise without some sort of transformation (we can normally only visualise 2-3 axes at once). Additionally, many SNPs will be highly correlated (linked), e.g. think height and weight, implying the a dimension-reduction method like PCA can be highly effective at reducing the dimensionality of the dataset while retaining most of the variance.

So we want to find a a rotation and scaling of the data to find axes that maximuse variation within the data. We can do this simply, e.g. in R or other tools (eg. SNP relate), however recall we have low coverage data, and our certainty in genotype calls and consequently variant calls will be high. So we want to work with genoytpe likelihoods.

Also (refer to PCAngsd website), missing data in these normal algorithms typically impute missing data as the the mean/zero, which leads to particular behaviour.

Add pic from PCAngsd

To perform PCA on lowcov data,  we have to infer the genetic covariance matrix, which can be estimated in different ways. Here we'll use PCANGSD. See lcs and pcangsd tutorials.

## Step 1: PCAngsd takes as input genotype likelihoods in beagle format, which we generated in the step before using the `-doGLF 2`option.

	REF=/data/Workshop/Data/assembly_homozygous.fa

Here, we have assigned the reference sequence fasta file to a variable, so that it is easy to refer to in the next steps (we can expand a variable *i* via ${*i*).

	angsd -GL 2 -out GL_95inds -ref ${REF} -nThreads 4 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -only_proper_pairs 1 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -bam samples.list -rf scaffolds.list
	angsd -GL 2 -out GL_75inds -ref ${REF} -nThreads 4 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -only_proper_pairs 1 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -bam samples_5inds.list -rf scaffolds.list


# Step 2: Run PCA via pcangsd. We can first look at PCAngsd's options via:
	
	pcangsd.py -h

And then perform the PCA (adjusting the number of threads accordingly)

	prefix="GL_95inds"
	prefix="GL_75inds"
	pcangsd.py -beagle ${prefix}.beagle.gz -threads 2 -o ${prefix}.pcangsd


# We can then plot the resutls (PC1 vs PC2) in R as follows (don't forget to upload pop metadata files for PCA, admixture, FST to your github!);
	R code

# How do we interpet the results (add html results to Github tutorial!). We find three distinct clusters. Let's add map to githubs totorial. How much variance is explained by the first two PCs? Provide links/references on how to interpret/not interpret PCA results

# Admixture/ancestry analysis

In some cases we also want to infer genome-wide admixture proportions
for each individuals. Similar to PCA, there are different ways of
inferring admixture proportions from genotype likelihoods. Here, we will
use [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and
will also present a way of inferring admixture proportions with PCAngsd.
Similar to the PCA, we want to use an LD-pruned SNP dataset for this
analysis to reduce the impact of non-independent SNPs on the ancestry
inference.

ngsAdmix uses a genotype likelihood file in beagle format (same as for
PCAngsd) as input, which is specified using the `-likes` option. In
addition, there are a range of parameters that can be adjusted. Here we
only set the number of ancestry clusters using the `-K` option to K=2.
In reality, it is advisable to compare different numbers of ancestry
clusters by iterating over different values of K.

In case the analysis is not converging, one can also increase the
maximum number of EM iterations using the `-maxiter` option.

ngsAdmix produces three different outputs files:

  - A run log containing the log likelihood of the admixture estimates:
    `.log file`
  - A file containing the allele frequencies of all ancestral
    populations (in this case two ancestral clusters): `.fopt file`
  - A file containing the admixture proportions for each individual:
    `.qopt file`

We are mostly interested in the admixture proportions and the log
likelihood of the estimates. The log likelihoods can be compared between
runs with different values of K to select the most likely number of
ancestral clusters (However, this should always be interpreted in a
biologically meaningful context)

In addition to ngsAdmix, PCAngsd can also be used to infer admixture
proportions. The command is similar to the one used for the PCA with the
addition of to a `-admix` option. The input file for ngsAdmix and
PCAngsd is the same, making comparisons relatively easy.

Other than ngsAdmix, PCAngsd can automatically infer the most likely
number of ancestry cluster. However, one can also set the number of
clusters using the `-admix_K` option.

Here, we have PCangsd output individual admixture proportions (-admix), and also output population specific allele frequencies (-admix_save). We iterate over K (here via the -e argument which defines the number of eigenvalues, rather than via -admix _K as the latter is not recommended).

To estimate admixture proportions in PCangsd, you need to define an alpha parameter (sparseness regularisation parameter). This can be specified manually (-admix_alpha) or automatically searching for the optimal alpha (-admix_auto), specifiying only a soft upper bound.
	for k in $(seq 1 2); do
		pcangsd.py -beagle ${prefix}.beagle.gz -threads 2 -e ${k} -admix -admix_auto 10000 -o ${prefix}.admix.pcangsd.K$((${k}+1))
	done

Alternatively, we can use NGSAdmix (slower)
Remember the LSB_JOBINDEX goes from 1,2,3,....
To change, add constant to variable e.g. $((${i}+c))
	for k in $(seq 1 2); do
		NGSadmix -likes ${prefix}.beagle.gz -K $((${k}+1)) -P 2 -o ${prefix}_$((${k}+1))
	done


Copy files to/from container (volume) to local filesystem (here we create a temporary container (named temp) with our named volume mounted)

	docker run --name temp --mount source=myvol3,target=/data -w /data ubuntu
	docker cp ./Desktop/Test/samples_5inds.list temp:/data/Workshop/Data/
	docker cp temp:/data/Workshop/Data/GL_95inds.pcangsd.cov ./Desktop/Test/
	docker cp temp:/data/Workshop/Data/admix_results ./Desktop/Test/
	#remove/detach docker contained after copying


# Site frequency spectrum and summary statistics

## For single populations.
	for i in $(seq 1 15); do
		pop=`sed -n ${i}p < pop5inds.list`
		pop_bamlist='./Population_lists_5inds/'${pop}.5inds
		# Note 1: BAQ adjusts (downgrades) quality scores arround indels, to reflect greater uncertainty around indel regions. Since we already perform GATK IndelRealigner, this is less necessary here. Otherwise, baq 2 (extended baq) would be preferable over baq 1. See: https://github.com/ANGSD/angsd/issues/97, http://bioinformatics.oxfordjournals.org/content/early/2011/02/13/bioinformatics.btr076, https://www.biostars.org/p/440490/, http://www.htslib.org/doc/samtools-mpileup.html
		# Note 2: Use GL 2! For BAMs which have been recalibrated and have had overlapping reads merged/clipped, GL 1 results in weirdly shaped SFS (see: http://www.popgen.dk/angsd/index.php/Glcomparison and SFS_troubleshooting.xlsx). 
		angsd -b ${pop_bamlist} -doSaf 1 -anc ${REF} -ref ${REF} -GL 2 -P 2 -minMapQ 1 -minQ 1 -C 50 -remove_bads 1 -only_proper_pairs 1 -doCounts 1 -setMinDepth 10 -setMaxDepth 75 -minInd 3 -out ${OUT}/${pop}
		# Obtain the maximum likelihood estimate of the SFS (here for the folded spectrum)
		realSFS ${OUT}/${pop}.saf.idx -P 2 -fold 1 > ${OUT}/${pop}.sfs
		# Calculate the thetas for each site
		realSFS saf2theta ${OUT}/${pop}.saf.idx -sfs ${OUT}/${pop}.sfs -outname ${OUT}/${pop}
		# Estimate Tajimas D and other statistics
		thetaStat do_stat ${OUT}/${pop}.thetas.idx
	done

## For two populations; population genetic differentiation g2 populations, FST, 2D-SFS, input for dadi/moments. PREPARE R PLOTTING SCRIPT FOR FST!
	# The following requires an input file which lists all pairwise comparisons, e.g. Airolo Bayasse
	# This can be produced by 1.Make_list_pairwise.sh
	for i in $(seq 1 15); do
		pop=`sed -n ${i}p < ${working_dir}/list_pop_name_pairs/pop_name_pairs`
		# Remember, set allows you to define the elements of your list as variables, according to their order
		set -- $pop
		##calculate the 2dsfs prior (choose folded or unfolded)
		realSFS ${input_dir}/${1}.saf.idx ${input_dir}/${2}.saf.idx -P 2 -fold 1 > ${output_dir}/${1}.${2}.ml
		##prepare the fst for easy window analysis etc. Option -whichFst 1 should be preferable for small sample sizes.
		realSFS fst index ${input_dir}/${1}.saf.idx ${input_dir}/${2}.saf.idx -sfs ${output_dir}/${1}.${2}.ml -fstout ${output_dir}/${1}.${2}.stats -whichFst 1
		#get the global estimate
		FST=`/cluster/project/gdc/shared/tools/angsd0_933/angsd/misc/realSFS fst stats ${output_dir}/${1}.${2}.stats.fst.idx`
		printf "${1}\t${2}\t${FST}\n" >> ${output_dir}/FST_summary.txt
		#get sliding window estimates
		realSFS fst stats2 ${output_dir}/${1}.${2}.stats.fst.idx -win 50000 -step 10000 > slidingwindow
	done

## 3 populations, PBS, 3D-SFS, input for dadi/moments
# Select 3 populations, i.e. one population each from each lineage. Code is similar to that above for FST (see ANGSD and low-cov tutorial website)
	realSFS fst index pop1.saf.idx pop2.saf.idx pop3.saf.idx -sfs pop1.pop2.ml -sfs pop1.pop3.ml -sfs pop2.pop3.ml -fstout out.pbs -whichFst 1

# Note that calling e.g. -n 4 and -P 16 allows better CPU usage as it allows for hyperthreading
	realSFS ${input_dir}/${1}.saf.idx ${input_dir}/${2}.saf.idx ${input_dir}/${3}.saf.idx -P 8 > ${output_dir}/${1}.${2}.${3}.sfs

## optional: dadi
## optional: selscan (requires vcfs)

