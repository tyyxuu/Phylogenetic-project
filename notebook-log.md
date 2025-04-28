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

**MrBayes Phylogenetic Analysis Notebook**

**Objective:**
To perform Bayesian phylogenetic inference using MrBayes on an aligned ITS1 region DNA dataset.

---

**Step 1: Prepare the Data**
- Aligned sequences were initially in FASTA format: `ITS1-region-aligned.fasta`
- Converted to NEXUS format using Biopython:

```python
from Bio import AlignIO
alignment = AlignIO.read("ITS1-region-aligned.fasta", "fasta")
for record in alignment:
    record.annotations["molecule_type"] = "DNA"
AlignIO.write(alignment, "ITS1-region-aligned.nex", "nexus")
```

---

**Step 2: Add MrBayes Command Block**
The following block was appended to the end of the `.nex` file:

```nexus
Begin mrbayes;
   set autoclose=yes nowarn=yes;
   lset nst=6 rates=gamma;
   outgroup KX192120.1;
   mcmc ngen=1000000 samplefreq=100 printfreq=1000 diagnfreq=1000 nchains=4 temp=0.2 savebrlens=yes;
   sumt burnin=2500;
End;
```

- `nst=6` uses the GTR model
- `rates=gamma` accounts for rate variation among sites
- `outgroup` designates the sequence used to root the tree
- `ngen=1000000` specifies 1 million MCMC generations
- `samplefreq=100` samples every 100 generations
- `burnin=2500` discards the first 25% of samples

---

**Step 3: Run MrBayes**
In the terminal, the MrBayes program was run from the compiled source directory:

```bash
~/MrBayes/src/mb
```

At the MrBayes prompt:

```mrbayes
execute ITS1-region-aligned.nex
```

---

**Step 4: Output Files**
After completion, MrBayes produces the following:

- `ITS1-region-aligned.nex.run1.p` and `.run2.p`: parameter traces
- `ITS1-region-aligned.nex.run1.t` and `.run2.t`: tree samples
- `ITS1-region-aligned.nex.con.tre`: **consensus tree** with posterior probabilities

---

**Next Steps:**
- Visualize the consensus tree using [FigTree](https://github.com/rambaut/figtree/releases)
- Check convergence and posterior summaries
- Optionally re-run analysis with more generations or model refinements

---

**Notes:**
- MrBayes was run on macOS, M1 architecture
- Execution time: approximately 10–30 minutes depending on sequence length and hardware
- Future extensions: partitioned models, clock models, or larger datasets
# Reproducible Note for Bayesian Phylogenetic Analysis (ITS1 Region)

## Project Information
- **Analysis Target:** ITS1-region aligned sequences
- **Goal:** Construct a Bayesian phylogenetic tree
- **Tools Used:**
  - BEAST 2.7.7
  - BEAUTi 2.7.7
  - TreeAnnotator 2.7.7
  - FigTree v1.4.5_pre

---

## 1. Sequence Data Preparation
- Input file: `ITS1-region-aligned.nex`
- Format: Aligned NEXUS file
- Aligned 12 taxa across 701 sites.

---

## 2. XML Configuration Using BEAUTi
- Imported `.nex` file into BEAUTi 2.
- Partition Settings:
  - **Data Type:** nucleotide
  - **Site Model:** GTR + Gamma (4 gamma categories, shape estimated)
  - **Clock Model:** Strict clock (rate estimated)
  - **Tree Prior:** Coalescent Constant Population
- Prior Distributions:
  - Base frequencies: Dirichlet[4.0, 4.0, 4.0, 4.0]
  - Gamma shape parameter: Exponential[1.0]
  - Population size: 1/X prior
  - Substitution rates: Gamma distributions with customized parameters
- MCMC Setup:
  - Chain length: 10,000,000
  - Pre-burnin: 0
- Output file: `ITS1-region-aligned.xml`

---

## 3. BEAST MCMC Analysis
- Command: Run BEAST 2.7.7 with the generated XML.
- Summary:
  - MCMC chain length: 10 million generations
  - Total runtime: ~278 seconds
  - End likelihood: `-1559.1203600913755`

---

## 4. Maximum Clade Credibility Tree Extraction
- Tool: TreeAnnotator 2.7.7
- Settings:
  - Burn-in: 10% (first 1000 trees discarded)
  - Tree target: Maximum clade credibility (MCC) tree
- Output file: `ITS1-region-aligned_annotated.tree`

---

## 5. Visualization and Plotting
- Tool: FigTree v1.4.5_pre
- Visualization Steps:
  - Layout: Rectangular tree
  - Tip Labels: Enabled (species names shown)
  - Branch Labels: Posterior probabilities displayed
  - Scale Bar: Enabled
  - Appearance: Line width increased, font size adjusted
- Export: Final tree exported as high-resolution PDF/PNG.

---

## 6. Notes on Posterior Probabilities
- Posterior probabilities close to 1 indicate strong branch support.
- Posterior values < 0.5 indicate weakly supported branches.
- Some branches had very low posterior (~0.03) and were treated cautiously during interpretation.

---

## ✅ Project Status
- Full analysis pipeline completed successfully.
- Final outputs:
  - MCC tree file
  - Annotated tree figure (PDF/PNG)

---

## Quick Summary
> This project reproducibly built a Bayesian phylogenetic tree based on ITS1-region sequences using BEAST2.7.7, and visualized the final tree with FigTree v1.4.5_pre.


