# Investigating-the-presence-of-missing-peptides-across-the-central-dogma-

# The Final_codes_for_CDS.ipynb file is a Google colab notebook containing the entire process pipeline of missing peptide identification and CDS mapping. This file can be uploaded in google colab and is ready to run. This notebook requires only one input file at the start of the code. The input file is the complete_uniprot_prot.fasta Supplementary file 1. 

# For the ncRNA_and_promoters.ipynb file, the input files are provided in supplementary files: Supplementary file 7 and supplementary file 15. The RNAcentral input file can be accessed via http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz. 

# The Secondary structure final code requires Supplementary file 4 or supplementary file 5 as inputs. 
# Before runnning the code, Vienna RNA fold must be installed in the linux environments. For users working on windows operating system, below are the steps. 

# wsl --install    Restart the system after installation.

# wsl --status     his installs Ubuntu by default on Windows systems

# LAUNCH UBUNTU IN WINDOWS

# sudo apt update

# sudo apt upgrade -y

# sudo apt install build-essential wget curl git -y    Install basic build utilities

# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh    Download Miniconda

# bash Miniconda3-latest-Linux-x86_64.sh    Run installer

# During installation:
# Press Enter to continue
  # Type yes to accept license
  # Choose installation location (default recommended)
  # Type yes when asked to initialize conda

# source ~/.bashrc     Restart terminal

# conda --version      Verify conda installation

# conda create -n rnafold_env python=3.13 -y     Create a new conda environment

# conda activate rnafold_env       Activate conda

# ViennaRNA is distributed through Bioconda. Configure channels in this order:
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

# conda install viennarna -y      Install ViennaRNA

# RNAfold --version    Check installation

# conda install pandas numpy scipy tqdm openpyxl matplotlib seaborn -y      Install necessary python packages

# Run SECONDARY STRUCTURE_FINAL_CODE.PY





# FOR THE CODE in Stop_codon_Reverse_peptide_ncORF.ipynb 
# Input files are Supplementary file 5: final_pep_frames_tranls.csv  &  Supplementary file 4: CQWW_frames.csv  for Stop codon analysis (Code cell 1 & 2) 
# Input file Supplementary file 5: final_pep_frames_tranls.csv for Reverse-peptide analysis (code cell 5)

# ncORF analysis
  # Input file is 41587_2022_1369_MOESM2_ESM.xlsx for ncORF analysis (code cell 8)
  # Input file: Supplementary file 5 for code cell 9
  # Input files for code cell 10:  1. 41586_2026_10459_MOESM13_ESM.xlsx  -> ORF-level protein evidence
  # 2. 41586_2026_10459_MOESM14_ESM.xlsx  -> peptide-level MS evidence

