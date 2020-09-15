
# the PRODIGY experiment on Mycobacterium tuberculosis (MTB) data to find rare mutations confering drug ressistance. 
The original R package was designed to prioritize driver genes for individual cancer patients. 

The details of the ORIGINAL method are described in
Dinstag G. & Shamir R. PRODIGY: personalized prioritization of driver genes. bioRxiv (2018), 
https://www.biorxiv.org/content/early/2018/10/30/456723

# Experiment rational and motivation (Oren Tzfadia, ITM Sept 2020)
We wish to currently focus our efforts on searching for drug-resistance rare mutations in the genome of Mycobacterium tuberculosis. We aim to do so, by combining whole genome sequencing (WGS), and gene expression profiles (RNAseq) from clinical samples. 
Reading the PRODIGY publication - made me think we can maybe try to harness cancer-tailored-methods such as PRODIGY (https://github.com/Shamir-Lab/PRODIGY), for hunting mutations in a handful of M. tuberculosis clinical samples. PRODIGY (written in R code), was designed originally to be cancer specific to find the best mutation combined with differential expression model and conclude that this model best explains (drives) cancer tumorigenesis. In case of MTB the effect of mutation is not if the patient will be developing the disease, but if mutation yes / not yields resistance to drug (antibiotics). My hypothesis is, if resistance is reflected by differential expression then by applying the same logic, we can argue that mutation combined with gene expression change (to suggest resistance), so the model should be working the same way to highlight such mutations. Hence, I think PRODIGY can highlight driver mutation in the bug (causing diff expression in the bug) leading to drug resistance (like driver mutations to develop cancer which is reflected by differential gene expression). Our basic logic is that the mutation causes cancer which is reflected by differential gene expression. However, if substantial changes in expression in the bug are expected to accompany resistance, we can use Prodigy. If not - another approach is needed).

## Workflow suggested: 
First we need to select mother - daughter samples from the BCCM collection, and cultur them for up to 6 weeks. Followed by genomic DNA extraction and sequencing. The resulting sequencing output (fastq data files), is then run through the MTBseq computationl pipeline which first map reads to a reference genome (often M. tuberculosis strain H37Rv) and then call genomic variants, creating a table of SNPs. The resulting SNP lists can then be used for a variety of analyses, such as strain typing, transmission clustering and drug resistance profiling. The results of these tasks are then reported to the end user (for example, a clinician or researcher; for more info see figure 3 at Meehan et al, 2019; https://www.nature.com/articles/s41579-019-0214-5).

## Required MTB data:
In sense of their input requirements, PRODIGY needs to be adjusted for MTB genome as follows: As input we need to collect SNPs and gene expression data per sample, and PPI networks and known pathways for MTB, with the interaction links between them (see paper Mei et al 2008 BMC Genomics). The reasoning for these data sets is to find the path/subnetwork from the tested mutation to the genes in the effected pathway, so we will need to have the underling networks and pathways to make the connection.

We have the pathways (Mycobrowser) and PPI networks (STRING https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4873-9), for MTB.
As for gene expression profiles, we need to keep in mind that proper control for differential expression is expected (it doesn't have to be the same patients of course as long as everything is normalized properly). We have samples from wild type ‘mother’ strains. These were subjected to sub-lethal antibiotics lives, and grew on plate. Then the colonies survived were analysed and found to have resistance mutations. We can extract RNA from both the ‘mother’ strains and the resistance ones (‘daughter’), in addition to the WGS (which already exists in BCCM collection). 

For now, we can use the examples data sets in the PRODIGY git. Then we can use DNA/RAN data from Inaki’s paper (REF). 

## Package installation
```r
library(devtools)
install_github("Shamir-Lab/PRODIGY")
```

PRODIGY was developed using and is dependent on the following packages (minimal version required is specified):

- MASS_7.3-50
- DESeq2_1.16.1
- igraph_1.2.2
- graphite_1.22
- ff_2.2-14
- plyr_1.8.4
- biomaRt_2.32.1
- PCSF_0.99.1
- mixtools_1.1.0
- ggplot2_3.0.0
- cowplot_0.9.3

## Simple run example
```r
library(PRODIGY)
# Load SNP+expression data derived from TCGA
data(COAD_SNV)
data(COAD_Expression)
# Load STRING network data 
data(STRING_network)
network = STRING_network
# Take samples for which SNP and expression is available 
samples = intersect(colnames(expression_matrix),colnames(snv_matrix))[1:5]
# Get differentially expressed genes (DEGs) for all samples
expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]
library(DESeq2)
DEGs = get_DEGs(expression_matrix,samples,sample_origins=NULL,beta=2,gamma=0.05)
# Identify sample origins (tumor or normal)
sample_origins = rep("tumor",ncol(expression_matrix))
sample_origins[substr(colnames(expression_matrix),nchar(colnames(expression_matrix)[1])-1,nchar(colnames(expression_matrix)[1]))=="11"] = "normal"	
# Run PRODIGY
all_patients_scores = PRODIGY_cohort(snv_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,
			pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins,write_results = F, results_folder = "./",beta=2,gamma=0.05,delta=0.05)
# Get driver gene rankings for all samples 
results = analyze_PRODIGY_results(all_patients_scores) 

# If files were stored in results_folder, read them to all_patients_scores first
all_patients_scores = list()
for(sample in samples)
{
	all_patients_scores[[sample]] = as.matrix(read.csv(file = paste(results_folder,sample,"_influence_scores.txt",sep=""),sep="\t",header=T,row.names=1))
	
}
#Load gold standard drivers from CGC. Here we used only genes that were annotated with a driver SNP by CGC.
data(gold_standard_drivers)
check_performances(results,snv_matrix,gold_standard_drivers)
```
