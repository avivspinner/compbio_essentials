# Living Document Instructions

This repository is paired with the [Open Source Tools in Computational Biology article](https://www.pillar.vc/playlist/article/open-source-tools-for-computational-biology/?utm_source=tldrnewsletter) that I wrote for Pillar VC. 
The .yml file will contain all essential packages needed to hit the (computational biology) ground running! 

## Contributing
Fork the repo, edit the README, and then make a PR.

# Open Source Tools for Computational Biology

The movement towards open source tools and knowledge for building biotech is important for the next wave of biotech innovation. We rounded up free open source software that founders, computational biologists, and data scientists are using that could help you get started on building products.

As a PhD student in computational biology, these are the basic development environments that I’m working in almost every day:

### Python

- **Jupyter Notebooks/Jupyter Lab**: very common, can split python code into code chunks. Jupyter Lab is especially useful for cloud computing because directories can be viewed and it creates different tabs for multiple notebooks/files.

- **Google Collab Notebooks**: extremely easy to collaborate with others and supports python and R.

### R
- **R Studio**: very easy-to-use development environment for programming in R
- **R Studio Cloud**: Like R studio, but all in the cloud.
### Bash
- **Terminal** for Macs and **Command Prompt** for Windows
- **Screen**: widely used terminal multiplexer that allows the user to view multiple sessions at a time.

## The Essential Packages
Regardless of specialty, there are several packages that are universally relevant. I will split these between python and R as well:

### Python
- [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [Matplotlib](https://matplotlib.org/) (and specifically matplotlib.pyplot), [Seaborn](https://seaborn.pydata.org/): fundamental package for numerical computing, analysis, and visualization
- [scikit-learn](https://scikit-learn.org/stable/): everything machine learning!

### R
- [tidyverse](https://www.tidyverse.org/): a powerful and well-documented set of packages for all general data manipulation and visualization. This includes: dplyr, ggplot2, tidyr, purrr, tibble, stringr, and forcats.
- [keras](https://keras.io/): everything deep learning!

Some disciplines have tools built for their specific field. We rounded up the most popular software packages in protein modeling, CRISPR, single cell RNA-seq, sequence alignments and computational chemistry. We created an environment (with a corresponding [.yml file](https://github.com/hanspinner/compbio_essentials) for easy activation! Simply download this file and create the environment with `conda env create -f /path/to/compbio_essentials.yml`. To activate the environment with `conda activate compbio_essentials`.) that contains all python packages listed below. Here is a [tutorial on how to connect it to JupyterLab](https://softwarejargon.com/jupyterlab-and-conda-environment-installation-and-setup/) and get coding!

### Protein Modeling
[EV Couplings](https://evcouplings.org/): open-source python package and web server for modeling proteins based on evolutionary couplings.
- Try running EV Couplings on HIV capsid protein (or any protein of your choice) by entering the UniProt id “F5BS86_9HIV1” into the [web server](https://v2.evcouplings.org/).

[EVE](https://github.com/debbiemarkslab/EVE): open-source python tool that leverages variational autoencoders to predict pathogenicity of variation in proteins.

[AlphaFold](https://alphafold.ebi.ac.uk/): open-source protein structure prediction server
- Try folding GFP (sequence linked [here](https://www.uniprot.org/uniprot/P42212.fasta)) in their [collab notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)!

[Rosetta Fold](): open-source protein structure prediction server

### CRISPR
[CRISPResso2](http://crispresso2.pinellolab.org/submission): analysis of genome editing outcomes from deep sequencing data with a fantastic [tutorial page](https://crispresso.pinellolab.partners.org/help).

[mali-dual-crispr-pipeline](https://github.com/ucsd-ccbb/mali-dual-crispr-pipeline): code for the dual-CRISPR screen analysis pipeline developed to analyze results from the dual-CRISPR screening system set up by the lab of Dr. Prashant Mali. The software is developed by Amanda Birmingham and Roman Sasik of the Center for Computational Biology and Bioinformatics at the University of California, San Diego.

### Single Cell RNA-seq
[scanpy](https://scanpy.readthedocs.io/en/stable/): Python-based toolkit for analyzing single cell transcriptomic data

[seurat](https://satijalab.org/seurat/): R-based toolkit for analyzing single cell transcriptomic data

[scVAE](https://github.com/scvae/scvae): Machine learning (Variational Autoencoder) to model single cell transcriptomics. They have really [helpful documentation here](https://scvae.readthedocs.io/en/stable/guide.html) that walks the user through how to use it (and you can skip the first step if you use the environment we already built!)

[kallisto](https://pachterlab.github.io/kallisto/):  ultra-fast single cell RNA-seq data analysis tool that doesn’t require full alignments.

[bustools](https://github.com/BUStools/bustools): tool for working with Barcode, UMI, and Set (BUS) files that can be used with kallisto.

### Genomics
[Hail](https://hail.is): provides accessible interfaces for exploring genomic data

[bedtools](https://bedtools.readthedocs.io/en/latest/): ‘swiss-army knife’ of genomic data analysis

[samtools](http://www.htslib.org/): general suite of programs for analyzing high-throughput sequencing data

[Anvi’o](https://anvio.org/): software suite for a wide array of disciplines within genomics including metagenomics, metatranscriptomics, pangenomics, metapangenomics, phylogenomics, and microbial population genetics. They also have a large swathe of tutorials and resources!

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi): basic alignment search tool that can be run through python or through their web interface.

[GATK](https://gatk.broadinstitute.org/hc/en-us): toolkit for identifying and discovering new single nucleotide polymorphisms (SNP) in DNA and RNAseq data.

[STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html): alignment tool specifically designed for for RNAseq data


### Image Analysis
[ImajeJ](https://imagej.net/software/fiji/)

### Synthetic Biology
[jecalles/synbio](https://github.com/jecalles/synbio): A simple python package for doing synthetic biology

[jecalles/z3helpers](https://github.com/jecalles/z3helpers): ergonomics around z3

[TimothyStiles/poly](https://github.com/TimothyStiles/poly): A Go package for engineering organisms



### General Biology
[Biopython](https://biopython.org/): broadly applicable python package for all things computational biology, especially handling and parsing different bio-specific file formats.

[Bioconductor](https://www.bioconductor.org/): broadly applicable R package, similar functional role as biopython.

[Snakemake](https://snakemake.readthedocs.io/en/stable/): workflow management system to create a robust computational pipeline.

[grailbio/reflow](https://github.com/grailbio/reflow): A language and runtime for distributed, incremental data processing in the cloud

[drc](https://cran.r-project.org/web/packages/drc/drc.pdf): R package to easily generate dose response curves, especially useful for ELISAs.

[GrowthCurver](https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html):  R package that summarizes and creates visualizations of optical density-based growth curves

### Computational Chemistry
[OpenMM](https://openmm.org/): toolkit for molecular simulation

[Amber](https://ambermd.org/): molecular dynamics software for DNA and proteins

[DeepChem](https://deepchem.io/): python package on deep learning for computational chemistry and biology

[TorchDrug](https://torchdrug.ai/): python package centered around deep learning platform for drug discovery

We would love to get your feedback on which tools you like, especially in fields we didn’t cover. You can reach out to Han at spinner@pillar.vc or check out their twitter thread for suggestions from the community.
