README for “Quantification of dynamic couplings between distal chromatin regions reveals functional correlations between genes”

First, please download the required data as described below. The script "Master_GNM_analysis_script.m" assumes that all of the following files are in the folder called "InputData" in the "MatlabCode" folder. The files needed are:

chr22_5kb.RAWobserved
chr22_5kb.VCnorm
GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed
ENCFF000SKV.pk
ENCFF002EMO.tsv

The first two are Hi-C data files for the GM12878 cell line from Rao et al., available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525 as part of the GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz file.
The bed file can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1155959.
The .pk file can be downlaoded from ENCODE (https://www.encodeproject.org/experiments/ENCSR000EJD/), with accession code ENCFF000SKV. 
The .tsv file is also found on ENCODE (https://www.encodeproject.org/experiments/ENCSR752QCX/), with accession code ENCFF002EMO.


Once the above files have been downloaded and are located in the "InputData" folder, simply running the “MatlabCode/Master_GNM_analysis_script.m” in MATLAB will go through all of the analysis on chromosome 22, and should take around 8 minutes to complete. First, the Kirchhoff matrix and its spectral decomposition will be computed, then the comparisons with experimental data and structural domains are carried out. In order to run the same analysis on other chromosomes, the Hi-C data from Rao et al. can be found at http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525. In order to compute TADs, use the Armatus tool which can be found at https://github.com/kingsfordgroup/armatus. Our implementation of the compartment computation method described in Lieberman-Aiden et al. (2009) can be found in the CompartmentAnalysis folder.

The co-expression analysis of CCDDs was done outside of MATLAB, so all of the files and code for this portion are found in the CoExpAnalysis folder, with descriptions below.


Co-expression analysis of CCDDs:

1. Download the experiments listed in the "SRRids.txt" file (also listed in Table S1 of Supp. Info.) from the SRA (http://www.ncbi.nlm.nih.gov/sra)

2. Quantify using Salmon (available at https://combine-lab.github.io/salmon/) with the Python script "script-runsalmon.py".

3. Using the "refFlat.txt" file (from the UCSC Genome Browser, http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) as a full reference of transcripts and genes from hg38, "writeTranscriptRefList.py" writes reference files for each gene and transcript, assuming the expression files written by Salmon are contained in a folder called “transcriptexpfiles”.

4. To convert the quantification files from transcript to gene level, use "transcript2genes.py". 

5. To calculate the co-expression and distance between each every pair, run “calcGeneDistances.py”.

6. Assuming you have already used the MATLAB function “findCCDDs.m” to write a file containing the locations of the CCDDs for every chromosome, the script “geneCorrelationsbyRegion.py” will then calculate the Pearson correlation of expression values for each gene pair in each CCDD. The input to the Python script should be the prefix of the CCDD files, i.e. ‘../MatlabCode/OutputFiles/ccddLocs_chr’

7. Finally, to plot the histograms comparing background gene pairs to those within the CCDDs, use “geneCorrelationsComparison.py”
