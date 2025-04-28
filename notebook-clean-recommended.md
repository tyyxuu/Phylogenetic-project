# Phylogenetic Analysis Project Notebook

## Project Description
The genus *Prunus*, which includes almonds, cherries, peaches, and plums, plays a critical role in agriculture, ecology, and human history.  
The ITS1 (Internal Transcribed Spacer 1) regions of nuclear ribosomal DNA are widely used in phylogenetic studies due to their high variability and suitability for resolving evolutionary relationships at the genus and species level.

---

## Data Resources
- ITS1 sequences obtained from GenBank, including:
  - *Prunus sp.* BSP-2004-2 isolates A1, A2, A6
  - *Prunus serrulata* var. *spontanea* T-26 (Accession: LC600963.1)
  - *Prunus dulcis* cultivar Lauranne (Accession: HF969277.1)

---

## Alignment Preparation

### 1. Software Setup
- Installed miniconda3
- Installed Anaconda
- Installed alignment tools:
  ```bash
  conda install -c bioconda clustalw
  conda install -c bioconda muscle

2. Multiple Sequence Alignment
Using ClustalW:

clustalw2 -ALIGN -INFILE=ITS1-region.fasta -OUTFILE=ITS1-region-aligned.fasta -OUTPUT=FASTA

Parsimony Tree Construction (RStudio)
library(ape)
library(phangorn)

dna <- fasta2DNAbin(file = "https://raw.githubusercontent.com/tyyxuu/Phylogenetic-project/main/Phylogenetic%20project%20dataset/ITS1-region-aligned.fasta")
print(dna)

D <- dist.dna(dna, model="TN93")
dna <- dna[apply(as.character(dna), 1, function(x) sum(x == "-")/length(x)) < 0.5, ]

tre <- tryCatch(nj(D), error = function(e) njs(D))

plot(tre, cex=0.6)
title("Neighbor-Joining Tree")

IQ-TREE Maximum Likelihood Tree (Visualization in R)
library(ape)

tree <- read.tree("ITS1-region-aligned.fasta.treefile")

plot.phylo(tree, type = "phylogram", cex = 0.8)
title("IQ-TREE Phylogenetic Tree")

png("phylogenetic_tree.png", width = 800, height = 600)
plot.phylo(tree, type = "phylogram", cex = 0.8)
dev.off()

Maximum Likelihood Analysis (phangorn R package)
library(ape)
library(phangorn)

dna <- fasta2DNAbin(file = "https://raw.githubusercontent.com/tyyxuu/Phylogenetic-project/main/Phylogenetic%20project%20dataset/ITS1-region-aligned.fasta")
dna <- dna[apply(as.character(dna), 1, function(x) sum(x == "-")/length(x)) < 0.5, ]

tree_init <- nj(dist.dna(dna, model = "JC"))

dna_phy <- phyDat(dna, type = "DNA")

fit <- pml(tree_init, data = dna_phy)
fit_ml <- optim.pml(fit, model="GTR", optGamma=TRUE, optInv=TRUE, rearrangement="stochastic")

plot(fit_ml$tree, main="Maximum Likelihood Tree (GTR+G+I)", cex=0.6)

bs <- bootstrap.pml(fit_ml, bs=100, optNni=TRUE)
plotBS(fit_ml$tree, bs, p=50, main="ML Tree with Bootstrap Support")

MrBayes Bayesian Phylogenetic Inference
Step 1: Data Preparation
Convert aligned FASTA to NEXUS:
from Bio import AlignIO

alignment = AlignIO.read("ITS1-region-aligned.fasta", "fasta")
for record in alignment:
    record.annotations["molecule_type"] = "DNA"
AlignIO.write(alignment, "ITS1-region-aligned.nex", "nexus")
Step 2: Add MrBayes Block
Begin mrbayes;
   set autoclose=yes nowarn=yes;
   lset nst=6 rates=gamma;
   outgroup KX192120.1;
   mcmc ngen=1000000 samplefreq=100 printfreq=1000 diagnfreq=1000 nchains=4 temp=0.2 savebrlens=yes;
   sumt burnin=2500;
End;
Step 3: Run MrBayes
~/MrBayes/src/mb
execute ITS1-region-aligned.nex
Step 4: Output Files
.run1.p and .run2.p: MCMC traces

.run1.t and .run2.t: Sampled trees

.con.tre: Consensus tree with posterior probabilities

BEAST2 Bayesian Phylogenetic Analysis
Step 1: XML Preparation (BEAUTi 2)
Imported ITS1-region-aligned.nex.

Configured:

Site Model: GTR + Gamma

Clock Model: Strict

Tree Prior: Coalescent Constant Population

MCMC:

Chain Length: 10,000,000

Output file: ITS1-region-aligned.xml

Step 2: BEAST2 MCMC Run
Run completed successfully.

Total runtime: ~278 seconds

End likelihood: ~ -1559.12

Step 3: Summarizing with TreeAnnotator
Burn-in: 10% (1000 trees discarded)

Maximum Clade Credibility Tree (MCC) selected.

Output file: ITS1-region-aligned_annotated.tree

Step 4: Visualization with FigTree
Layout: Rectangular

Tip Labels: Species names shown

Branch Labels: Posterior probabilities displayed

Scale Bar: Enabled

Exported final tree figure as PDF and PNG.

Final Notes
Posterior probabilities near 1 indicate strong clade support.

Posterior < 0.5 indicate low support and were interpreted cautiously.

The BEAST2 Bayesian pipeline completed reproducibly.

Project Status: Completed Successfully
