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
To count: grep ">" ITS1-region.fasta | wc -l
clustalw2 -ALIGN -INFILE=primatesAA.fasta -OUTFILE=primatesAA-aligned.fasta -OUTPUT=FASTA

clustalw2 -ALIGN -INFILE=ITS1-region.fasta -OUTFILE=ITS1-region-aligned.fasta -OUTPUT=FASTA

RStudio Parsimony:
# Load necessary libraries
library(ape)       # For phylogenetic tree construction and DNA handling
library(phangorn)  # For additional phylogenetic methods

# Load aligned FASTA file directly from GitHub
dna <- fasta2DNAbin(file = "https://raw.githubusercontent.com/tyyxuu/Phylogenetic-project/main/Phylogenetic%20project%20dataset/ITS1-region-aligned.fasta")

# Print basic information about the loaded sequence data
print(dna)  # This helps verify that the sequences were read correctly

# Compute a genetic distance matrix using the TN93 model
D <- dist.dna(dna, model="TN93")  

# Display nucleotide composition to check for gaps or ambiguous bases
table(as.character(dna))  

# Remove sequences that have more than 50% gaps ('-') to clean the dataset
dna <- dna[apply(as.character(dna), 1, function(x) sum(x == "-") / length(x)) < 0.5, ]

# Recalculate the distance matrix using the "raw" model with pairwise deletion
# - "raw" model treats each site independently without assuming an evolutionary model
# - "pairwise.deletion=TRUE" ensures that missing data (gaps) are ignored when computing distances
D <- dist.dna(dna, model="raw", pairwise.deletion=TRUE)

# Construct a Neighbor-Joining (NJ) tree
# - `nj(D)` is the standard NJ method
# - If `nj(D)` fails (due to missing values or data issues), use `njs(D)` as a fallback
tre <- tryCatch(nj(D), error = function(e) njs(D))

# Plot the phylogenetic tree
plot(tre, cex=0.6)  # Adjust text size for readability

# Add a title to the tree plot
title("Neighbor-Joining Tree")  
