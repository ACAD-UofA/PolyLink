
# `PolyLink: gene-based pathway enrichment` <img src="inst/sticker/polylinkr_150px.png" align="right" />
Is an R package that performs gene-based pathway enrichment, which can also be used as evidence for polygenic selection in case the software is used with selection signals evidence. The package explicitly also accounts for linkage desiquilibrium between adjacent loci belonging on the same pathway.

While we use the same input format, and a good number of functions from [POLYSEL](https://github.com/CMPG/polysel), the core algorithms are substantially different in PolyLink

## Software requirements
The following scripts expect that [SweepFinder2](http://www.personal.psu.edu/mxd60/sf2.html), and [Plink](https://www.cog-genomics.org/plink/1.9) (>=v1.9) are in the path, as well as a fully working version of R and the called packages:

- Matrix
- RColorBrewer
- data.table
- doParallel
- foreach
- ggplot2
- igraph
- matrixStats
- qvalue
- stringi

## Example: Using the software with selection scans
PolyLink can work with a range of analyses. For example, we can use selection signal scores to check if we have a pathway enrichment and evidence for polygenic selection at the pathway level. Here we illustrate using selection scan example using [SweepFinder2]() with ancient DNA dataset (see Examples). 

### Notes on data format
The bash wrapper script that runs creates SweepFinder2 input (SFS) expects the data in Plink format. You can get that by using `convertf` utility available as part of the AdmixTools suite. However, make sure that the `.fam` is properly formated with the population IDs in the first column, and that the phenotype column is setup to either `1`, or `0` and not `-9`. Also, sex chromososmes are ignored and can be removed from the dataset.

Also, needed is a refernce file (fasta). Plink has the upleasant habit of polarising each SNP based on the minor allele frequency whithin the dataset. To remedy that prior to the generation of the SFS, we repolarise each SNP using the reference allele using bcftools (load the module if necessary).

In the shared repository there's an example file for a subset of the 1240k SNP set. As the polarisation step could be accomplished with Plink 2.x if used (--update-ref-allele refAllele.txt).

### Running SweepFinder2
The warapper (`sweepfinder2_splitByPop.sh`  under `scripts` folder)script is a simple bash that creates an SFS file for each chromosome for each population in the plink file. Also, make sure to edit the file to point to the righ location of the `vcf2SF.py` python script (`scripts` folder), which depend on ***python3*** The script could be used as follows:

```
$ bash sweepfinder2_splitByPop.sh <plink file basename> <ref.fasta>
```

**Note:** Please note that sweepfinder jobs (per chromosome, per population) takes a fair amount of time to run, and scales up principally with the number of loci. Each job takes for 1M SNPs between _8-12h walltime_. So it is advised to replace the SweepFinder2 calls (in the bash script) with a job submission. And/or replace the loop inside the script with `parallel` or `xargs` to activate parallelisation. Each job requires 2 cores and 2GB of ram.


### Analysing SweepFinder results
#### Score genes an determine outliers
The SweepFinder outputs are located in SweepFinder folder. There is one folder for each population. Each population folder contains a folder called 'out', which contains the SweepFinder CLR scores used in later analyses. The other folder contains allele frequencies used to estimate the effective population size.

The Raw_files folder contains two files that are used in the first step – assigning the SweepFinder CLR scores to individual genes for each population `Entrez_Gene_IDs_positions_HG19.txt` – and the second step – creating the files for the PolyLink pathway enrichment analysis `PolyLink_GO_pathways_HG19.txt`. Both files were created using custom R scripts that can be made available, but feel free to generate your set of gene and pathway files.

To perform the first step, you will need to run the script `Determine_outlier_genes_and_sweeps.R`. All R scripts can be found in the Scripts folder. Running this script will generate two new files with CLR scores assigned to genes `GeneScores_max_LR.txt`, and scores corrected for gene length `GeneScores_max_LR_GeneLengthCorrected.txt`. These output files will end up in the `GeneScores` folder.

Also, the `Determine_outlier_genes_and_sweeps.R` will create a set of sweeps by clustering sets of 'outlier' genes. Outlier genes are determined by assigning p and q values to each gene for each population, and using a q-value threshold to identify outliers. The gene scores have an approximately standard normal distribution (i.e. Z scores), allowing simple calculation of p values. The q values calculation uses Storey's method, but note that we only calculate this for the upper tail of the distribution (i.e. Z.0), since the p-value distribution is U-shaped violating one of the assumptions of the q-value estimation. Once the outlier genes are identified, these are combined into sweep regions if the distance between the two outlier genes is less than a pre-specified value. We have used a distance of 1Mb for our sweeps, and it seems to work fairly well.

## Pathway analysis (PolyLink)
Once the gene scores have been generated, the pathway analysis can be performed. The pathway analysis is based on a modified version of PolyLink that accounts for genomic clustering amongst genes in a pathway using a permutation method that we have devised. To create the input files for the pathway analysis, use the script named `Create_PolyLink_input_files.R` in the Scripts folder. This will generate a set of three files for each of the populations you are testing, containing information about the gene scores, the pathway names, and how genes are allocated to each pathway. The input files can be found in the `PolyLink` folder, in the `Input` subfolder.

Once the input files are generated, you can run the pathway enrichment analysis using the `PolyLink_CMD.R` in the Scripts folder. This is implemented to run from the command line, and is reasonably fast (~1-2hrs on a 4 core laptop). The default settings run 1M permutations, which always converge for the populations we have worked with so far. The resulting output files can be found in the 'Output' subfolder in the `PolyLink` folder.

Finally, the pathway enrichment output files can be analysed using the `Analyse_PolyLink_output.R` script in the Scripts folder. This will provide a set of enriched pathways defined at a specific *q-value* threshold. We have been using a cutoff of **q<0.05** for our populations.

Note, that I have included a set of population data in this SweepFinder folder that can be used to run the scripts and make sure the everything is working correctly. I have already run the scripts locally, and the outputs can be found in the relevant folders. So you may find it useful to rerun the whole process to see if you get the same outputs (noting that the pathway *p-values* will be slightly different due to the permutation based method). If you choose to do this, make a copy of the output files first, since the current setup will overwrite the existing files.
