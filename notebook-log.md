Project Description: The genus Prunus, which includes important species such as almonds, cherries, peaches, and plums, plays a critical role in agriculture, ecology, and human history. ITS1 (Internal Transcribed Spacer 1) regions of nuclear ribosomal DNA are widely used in phylogenetic studies due to their high variability and suitability for resolving evolutionary relationships at the genus and species level.

Data resources:ITS1 sequences were obtained from GenBank, including Prunus sp. BSP-2004-2 partial ITS1, isolates A1, A2, and A6. Additional ITS1 sequences from Prunus serrulata var. spontane T-26 (Accession: LC600963.1) and Prunus dulcis cultivar Lauranne (Accession: HF969277.1) were included to provide broader taxonomic context.

Alignment installation:
1. Installed miniconda3
2. Installed Anaconda
3. conda install -c bioconda clustalw
conda install -c bioconda t-coffee (Not working for some reason)
conda install -c bioconda muscle

Align it with clustalw:
Cd to where dataset is, then run
To count: grep ">" Prunus-Dulcis-ITS1.fasta | wc -l
clustalw2 -ALIGN -INFILE=primatesAA.fasta -OUTFILE=primatesAA-aligned.fasta -OUTPUT=FASTA

clustalw2 -ALIGN -INFILE=Prunus-Dulcis-ITS1.fasta -OUTFILE=Prunus-Dulcis-ITS1-aligned.fasta -OUTPUT=FASTA