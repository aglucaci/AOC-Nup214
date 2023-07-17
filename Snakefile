"""

2023 - Analysis of Orthologous Collections (AOC).
@Author: Alexander G. Lucaci

"""

# =============================================================================
# Imports
# =============================================================================

import itertools
import os
import sys
import csv
import json
from pathlib import Path
from snakemake.utils import min_version
min_version("7.20")
from Bio import Entrez
from ete3 import NCBITaxa
import pandas as pd
from ete3 import Tree
import glob
#from tqdm import tqdm

# =============================================================================
# Configuration
# =============================================================================

configfile: 'config.yml'

print("# Loaded configuration YAML file from", 'config.yml')

cluster_json = "cluster.json"

with open(cluster_json, "r") as fh:
  cluster = json.load(fh)
  fh.close()
#end with

print("# Loaded cluster configuration JSON file:", cluster_json)
 
Label = config["Label"]

Taxon = config["Taxon"]

# Get working directory
BASEDIR = os.getcwd()

TEST_MODE = False

if TEST_MODE == True:
    Nucleotide_file = os.path.join(BASEDIR, 
                                   "test", 
                                   config["TEST_NUCLEOTIDE"])
    
    Protein_file    = os.path.join(BASEDIR, 
                                   "test", 
                                   config["TEST_PROTEIN"])  
else:
    Nucleotide_file = os.path.join(BASEDIR, 
                                   "data",  
                                   Taxon,
                                   config["Nucleotide"])
    
    Protein_file    = os.path.join(BASEDIR, 
                                   "data", 
                                   Taxon, 
                                   config["Protein"])
#end if

CSV = os.path.join(BASEDIR, 
                   "data", 
                   Taxon, 
                   config["CSV"])

print("# We are operating out of base directory:", 
      BASEDIR)

print("# Using nucleotide data from:", 
      Nucleotide_file)

print("# Using protein data from:", 
      Protein_file)

print("# Using the analysis label:", 
      Label)

print("# Using the CSV Ortholog file:", 
      CSV)

# Create output directories
os.makedirs(os.path.join(BASEDIR, 
                         "results"), 
                         exist_ok = True)

OUTPUT_DIR = os.path.join(BASEDIR, 
                          "results", 
                          Taxon)

os.makedirs(OUTPUT_DIR, 
            exist_ok = True)

#print("# Directory '% s' created" % OUTPUT_DIR)

print("# We will create and store results in:", 
      OUTPUT_DIR)

# =============================================================================
# Compute settings
# =============================================================================

PPN = cluster["__default__"]["ppn"]

# =============================================================================
# Batch files
# =============================================================================

PREMSA = os.path.join(BASEDIR, 
                      config["PREMSA"])

POSTMSA = os.path.join(BASEDIR, 
                       config["POSTMSA"])

FILTER_OUTLIERS_BF = os.path.join(BASEDIR, 
                                  "hyphy-analyses", 
                                  "find-outliers", 
                                  "find-outliers-slac.bf")

FITMG94 = os.path.join(BASEDIR, 
                       "hyphy-analyses", 
                       "FitMG94", 
                       "FitMG94.bf")

# =============================================================================
# Hard-coded HyPhy settings
# =============================================================================

HYPHY = config["HYPHY"]
HYPHYMPI = config["HYPHYMPI"]
#RES = os.path.join(BASEDIR, config["RES"])

# =============================================================================
# Clustering threshold hyperparameter
# =============================================================================

TN93_T = config["TN93_Threshold"]

# =============================================================================
# Rule all
# =============================================================================

rule all:
    input:
        os.path.join(OUTPUT_DIR, "Codons.fa"),
        os.path.join(OUTPUT_DIR, "Codons.fa.log"),
        os.path.join(OUTPUT_DIR, "UTR5.fa"),
        os.path.join(OUTPUT_DIR, "UTR3.fa"), 
        os.path.join(OUTPUT_DIR, "Codons.fa_protein.fas"),
        os.path.join(OUTPUT_DIR, "Codons.fa_nuc.fas"),
        os.path.join(OUTPUT_DIR, "Codons.fa_protein.aln"),
        os.path.join(OUTPUT_DIR, "Codons.fa_codons.fasta"),
        os.path.join(OUTPUT_DIR, "Codons.fa_codons_duplicates.json"),
        os.path.join(OUTPUT_DIR, 
                     Label + "_" + Taxon + "_codons.SA.fasta"),
        os.path.join(OUTPUT_DIR, 
                     Label + "_" + Taxon + "_codons.SA.fasta.treefile")
#end rule all

print("# Moving on to processing rules")

# =============================================================================
# Rules
# =============================================================================

rule get_codons:
    output:
        codons  = os.path.join(OUTPUT_DIR,
                              "Codons.fa"),
        UTR5    = os.path.join(OUTPUT_DIR, 
                              "UTR5.fa"),
        UTR3    = os.path.join(OUTPUT_DIR, 
                              "UTR3.fa"),
        Logfile = os.path.join(OUTPUT_DIR,
                              "Codons.fa.log"),
        #ProteinLengthsCSV = os.path.join(OUTPUT_DIR, 
        #                                 "ProteinLengths.csv")
    params:
        Nucleotides = Nucleotide_file,
        Proteins    = Protein_file,
        Codons      = os.path.join(OUTPUT_DIR,
                              "Codons.fa"),
        Logfile     = os.path.join(OUTPUT_DIR,
                              "Codons.fa.log"),
        UTR5        = os.path.join(OUTPUT_DIR, 
                              "UTR5.fa"),
        UTR3        = os.path.join(OUTPUT_DIR, 
                              "UTR3.fa"),
        ProteinLengthsCSV = os.path.join(OUTPUT_DIR, 
                                         "ProteinLengths.csv")
    script:
        "scripts/codons.py"
        #"scripts/CodonStats.py"
#end rule

# =============================================================================
# Alignment
# =============================================================================

rule pre_msa:
    input: 
        codons = rules.get_codons.output.codons
    output: 
        protein_fas    = os.path.join(OUTPUT_DIR, 
                                      "Codons.fa_protein.fas"),
        nucleotide_fas = os.path.join(OUTPUT_DIR, 
                                      "Codons.fa_nuc.fas")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} {PREMSA} --input {input.codons}"
#end rule 

rule mafft:
    input:
        protein = rules.pre_msa.output.protein_fas
    output:
        protein_aln = os.path.join(OUTPUT_DIR, 
                                   "Codons.fa_protein.aln"),
    shell:
        "mafft --auto {input.protein} > {output.protein_aln}"
#end rule

rule post_msa:
    input: 
        protein_aln = rules.mafft.output.protein_aln,
        nucleotide_seqs = rules.pre_msa.output.nucleotide_fas  
    output: 
        codons_fas = os.path.join(OUTPUT_DIR, 
                                  "Codons.fa_codons.fasta"),
        duplicates_json = os.path.join(OUTPUT_DIR, 
                                       "Codons.fa_codons_duplicates.json"),

    shell: 
        "mpirun -np {PPN} {HYPHYMPI} {POSTMSA} --protein-msa {input.protein_aln} --nucleotide-sequences {input.nucleotide_seqs} --output {output.codons_fas} --duplicates {output.duplicates_json}"
#end rule 

# =============================================================================
# Remove ambiguous codons, a source of noise.
# =============================================================================

rule strike_ambigs:
   input:
       in_msa = rules.post_msa.output.codons_fas
   output:
       out_strike_ambigs = os.path.join(OUTPUT_DIR, 
                                        Label + "_" + Taxon + "_codons.SA.fasta")
   shell:
       "{HYPHY} scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule 

# =============================================================================
# IQ-TREE for ML tree inference
# =============================================================================

rule iqtree: # Unrooted
    input:
        codons_fas = rules.strike_ambigs.output.out_strike_ambigs
    output:
        tree = os.path.join(OUTPUT_DIR, 
                     Label + "_" + Taxon + "_codons.SA.fasta.treefile")
    shell:
        "iqtree -s {input.codons_fas} -T AUTO -B 1000"
    #end shell
#end rule iqtree

# =============================================================================
# End of file
# =============================================================================

