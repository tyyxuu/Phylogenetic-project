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

IQ-Tree alignment 
# Install the 'ape' package (only needed once)
install.packages("ape")  

# Load the 'ape' package for phylogenetic analysis
library(ape)  

# Read the Newick tree from the tree file
tree <- read.tree("ITS1-region-aligned.fasta.treefile")  

# Plot the phylogenetic tree as a standard phylogram
plot.phylo(tree, type = "phylogram", cex = 0.8)  

# Add a title to the tree plot
title("Phylogenetic Tree")  

# Save the tree as a PNG file
png("phylogenetic_tree.png", width = 800, height = 600)  
plot.phylo(tree, type = "phylogram", cex = 0.8)  
dev.off()  # Close the PNG graphics device and save the file

# Display edge labels (branch lengths, bootstrap values if present)
edgelabels()  

# Display node labels (internal nodes of the tree)
nodelabels()  

# Attempt to root the tree at node 22 (node numbers are internal tree structure identifiers)
tree2 = root(tree, node= 22)  

# Plot the newly rooted tree
plot(tree2)  

# Replot the original tree for comparison
plot(tree)  

# Display edge and node labels again to check their placements
edgelabels()  
nodelabels()  

# Attempt to root the tree at node 1 (incorrect: node number should be greater than the number of taxa)
tree2 = root(tree, node= 1)  # This produces an error

# Attempt to root the tree at node 13 (valid node number)
tree2 = root(tree, node= 13)  
plot(tree2)  # Plot the new rooted tree

# Replot the original tree for comparison
plot(tree)  
edgelabels()  
nodelabels()  

# Attempt to root the tree at node 20
tree2 = root(tree, node = 20)  
plot(tree2)  # Plot the re-rooted tree

# Replot the original tree one last time
plot(tree)  
edgelabels()  
nodelabels()  

# Final plot of the re-rooted tree
plot(tree2)  


---Maximum Likelihood

## Step 1: Load Required Packages and Data

```r
library(ape)
library(phangorn)

# Load aligned ITS1 sequences
dna <- fasta2DNAbin(file = "https://raw.githubusercontent.com/tyyxuu/Phylogenetic-project/main/Phylogenetic%20project%20dataset/ITS1-region-aligned.fasta")
print(dna)

# Initial distance matrix using TN93 model
D <- dist.dna(dna, model = "TN93")

# Visualize character composition (for gap inspection)
table(as.character(dna))

# Filter sequences with >50% gaps
dna <- dna[apply(as.character(dna), 1, function(x) sum(x == "-") / length(x)) < 0.5, ]

# Recalculate distance using raw model with pairwise deletion
D <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE)

# Build NJ tree
tre <- tryCatch(nj(D), error = function(e) njs(D))

# Plot NJ tree
plot(tre, cex = 0.6)
title("Neighbor-Joining Tree")

# Convert DNA alignment to phyDat format
dna_phy <- phyDat(dna, type = "DNA")

# Use NJ tree (JC model) as starting tree
tree_init <- nj(dist.dna(dna, model = "JC"))

# Build ML tree using GTR+G+I model
fit <- pml(tree_init, data = dna_phy)
fit_ml <- optim.pml(fit, model = "GTR", optGamma = TRUE, optInv = TRUE, rearrangement = "stochastic")

# Plot ML tree
plot(fit_ml$tree, main = "Maximum Likelihood Tree (GTR+G+I)", cex = 0.6)

# Perform 100 bootstrap replicates
bs <- bootstrap.pml(fit_ml, bs = 100, optNni = TRUE)

# Plot ML tree with bootstrap support (threshold = 50%)
plotBS(fit_ml$tree, bs, p = 50, main = "ML Tree with Bootstrap Support")


