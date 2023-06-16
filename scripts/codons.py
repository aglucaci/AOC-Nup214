#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
codons.py

@Author: Alexander G. Lucaci
 
"""

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
import os
import sys
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from tqdm import tqdm

# =============================================================================
# Declares
# =============================================================================

PROTEIN_FASTA     = snakemake.params.Proteins
TRANSCRIPTS_FASTA = snakemake.params.Nucleotides
CODONS            = snakemake.params.Codons
Logfile           = snakemake.params.Logfile
UTR5              = snakemake.params.UTR5
UTR3              = snakemake.params.UTR3

results = []
results_UTR5 = []
results_UTR3 = []
no_match = []
successful_count = 0
num_errors = 0
errors_IDs = []

# =============================================================================
# Helper functions
# =============================================================================

def already_in_results(transcript_desc, results):
    Found = False
    for record in results:  # results stores transcript records that passed
        if transcript_desc == record.description:  # already exists?
            Found = True
            break
        # end if
    # end for
    return Found
# end method

def log(msg, Logfile):
    with open(Logfile, 
              "a") as fh2:
        print(msg, 
              file=fh2)
        fh2.close()
    # end with
# end method

def Process(protein_desc, 
            protein_seq, 
            TRANSCRIPTS_FASTA, 
            species, 
            seq_threshold=None):

    global results, no_match, Logfile
    
    # Loop over all of the TRANSCRIPTS_FASTA sequences
    with open(TRANSCRIPTS_FASTA, "r") as transcript_handle:
        for m, transcript_record in enumerate(SeqIO.parse(transcript_handle, 
                                                          "fasta")):
            DONE = False
            
            # Grab Transcript Data
            transcript_id = transcript_record.id
            transcript_desc = transcript_record.description
            transcript_seq = transcript_record.seq
            
            exists = already_in_results(transcript_desc, 
                                        results)

            if species not in transcript_desc:
                # print("# Mismatch between species") 
                # move on to the next one
                # only look at sequences from your species, 
                # not something similar.
                continue
            # end if

            start = 0
            NT_SEQ_LENGTH = len(protein_seq) * 3
            # It should not get less than NT_SEQ_LENGTH
            # Doesn't need to run to zero.
            
            #print("## Starting to translate sequence subsets", transcript_desc)
            
            while start < len(str(transcript_seq)):
                coding_dna = ""
                try:
                    coding_seq = transcript_seq[start: start + NT_SEQ_LENGTH]
                    coding_dna = coding_seq.translate()  # translated, 
                                                         # universal code
                except:
                    pass
                # end try
                
                # Exit upon first match, 
                # may be useful to see how many matches.
                #print(len(coding_dna), len(protein_seq))
                if coding_dna == str(protein_seq) and exists == False:
                    DONE = True
                    break
                else:
                    start += 1
                # end if
            # end while
            if DONE == True:
                break
            #end if
        # end for
    # end with

    if DONE == True:
        # return transcript_id, transcript_desc, coding_seq
        transcript_record.seq = coding_seq
        
        """
        record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
                        IUPAC.protein),
                   id="YP_025292.1", name="HokC",
                   description="toxic membrane protein")
        """
        UTR5 = transcript_seq[0: start]
        tx_record_UTR5 = SeqRecord(Seq(UTR5),
                                   id=transcript_id, 
                                   description=transcript_desc)
        
        UTR3 = transcript_seq[start + NT_SEQ_LENGTH:]
        tx_record_UTR3 = SeqRecord(Seq(UTR3),
                                   id=transcript_id, 
                                   description=transcript_desc)
        
        #tx_record_UTR5 = transcript_record
        #UTR5 = transcript_seq[0: start + 1]
        #transcript_record.seq = UTR5
        
        #tx_record_UTR3 = transcript_record
        #UTR3 = transcript_seq[start + NT_SEQ_LENGTH + 1:]
        #tx_record_UTR3.seq = UTR3
        
        return transcript_record, tx_record_UTR5, tx_record_UTR3
    else:
        return "NO_MATCH", "", ""
    # end if
# end method

def CheckCounts(PROTEIN, 
         TRANSCRIPTS):  
    # Really to verify things number of seqs match
    
    print("# TRANSCRIPT INPUT FILE:", TRANSCRIPTS)
    print("# PROTEIN INPUT FILE:", PROTEIN)
    protein_list = []
    transcript_list = []

    with open(TRANSCRIPTS, "r") as handle:
        trans_count = 0
        for record in SeqIO.parse(handle, "fasta"):
            trans_count += 1
            transcript_list.append(record.description)
        #end for
        print("# Number of transcripts:", trans_count)
        handle.close()
    #end with

    with open(PROTEIN, "r") as handle:
        prot_count = 0
        for record in SeqIO.parse(handle, "fasta"):
            prot_count += 1
            protein_list.append(record.description)
        #end for
        print("# Number of proteins:", prot_count)
        handle.close()
    #end with

    return trans_count, prot_count
# end method

# =============================================================================
# Main subroutine
# =============================================================================

#print("# Processing... ", "\n")

# Get counts of Transcripts and Proteins
#transcript_count, protein_count = CheckCounts(PROTEIN_FASTA,
#                                             TRANSCRIPTS_FASTA)

# We make an assumption here, based on the NCBI Orthologs DB 
# that all species match up in transcript and protein fasta 

# This is exceptional, will need to look for species name 
# (from protein desc.) in transcript desc.

# =============================================================================
# Empty all output files
# =============================================================================

# Empty the log file
with open(Logfile, "w") as fh:
    print("", file=fh)
    fh.close()
# end with

# Create empty output file.
with open(CODONS, "w") as fh:
    fh.write("")
fh.close()

with open(UTR5, "w") as fh:
    fh.write("")
fh.close()

with open(UTR3, "w") as fh:
    fh.write("")
fh.close()

# =============================================================================
# Log
# =============================================================================
with open(Logfile, "a") as fh2:
    print("# Writing CDS data to:", CODONS, file=fh2)
    fh2.close()
# end with

# =============================================================================
# Get average sequence length.
# =============================================================================
"""
prot_seq_lengths = []

with open(PROTEIN_FASTA, "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, "fasta")):
        prot_seq_lengths.append(len(record.seq.ungap("-")))
    # end for
# end with

prot_handle.close()

avg_sequence_length = sum(prot_seq_lengths) / len(prot_seq_lengths)

avg_sequence_length_nt = avg_sequence_length * 3

with open(logfile, "a") as fh2:
    # print()
    # print("# Processing:", protein_desc, file=fh2)
    print("# Average sequence length is (PROTEIN AA ungapped):",
          avg_sequence_length, file=fh2)
    print("# Average sequence length is (NUCLEOTIDE NUC):",
          avg_sequence_length_nt, file=fh2)
    fh2.close()
# end with
"""

# =============================================================================
# Grab the protein, can I find a match in any of the mRNA transcript?
# =============================================================================
with open(PROTEIN_FASTA, 
          "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, 
                                           "fasta")):
        
        # Grab protein data.
        protein_id = record.id
        protein_desc = record.description
        protein_seq = record.seq
        
        print("#", 
              n+1, 
              " Processing:", 
              protein_desc)

        """
        # Send to log
        with open(logfile, "a") as fh2:
            print("\n", 
                  file=fh2)
            print("# Processing:", 
                  protein_desc, 
                  file=fh2)
            fh2.close()
        # end with
        """
        
        # Quality filter
        FILTERS = ["LOW QUALITY PROTEIN", "PARTIAL"]

        for item in FILTERS: 
            if item in str(protein_desc).upper():
                with open(Logfile, 
                          "a") as fh2:
                    
                    print("# Skipping this sequence due to quality issues",
                          record.description, file=fh2)
                # end with
                fh2.close()
                continue
            # end if
        #end for

        species = ""

        if "[" in protein_desc:
            species = record.description.split("[")[1].replace("]", 
                                                               "")
            #print(species)
        # end if

        """
        with open(logfile, "a") as fh2:
            print("# Species:", species, file=fh2)
            fh2.close()
        # end with
        """
        
        # Heavy lifting here.
        #tx_record = Process(protein_desc, 
        #                    protein_seq.ungap("-"), 
        #                    TRANSCRIPTS_FASTA, 
        #                    species, 
        #                    int(avg_sequence_length_nt))
        
        tx_record, tx_UTR5, tx_UTR3 = Process(protein_desc, 
                                        protein_seq.ungap("-"), 
                                        TRANSCRIPTS_FASTA, 
                                        species)

        if type(tx_record) != str:
            with open(Logfile, 
                      "a") as fh2:
                print("# Match:", 
                      tx_record.description, 
                      "\n", 
                      file=fh2)
                fh2.close()
            # end with
            
            results.append(tx_record)
            results_UTR5.append(tx_UTR5)
            results_UTR3.append(tx_UTR3)
        else:
            """
            # No match...
            with open(logfile, 
                      "a") as fh2:
                print("# -- NO Match -- \n", 
                      file=fh2)
            # end with
            fh2.close()
            """
            no_match.append(protein_desc)
    # end for
# end with

# =============================================================================
# Report on no matches, this needs to be written to a log
# =============================================================================
with open(Logfile, "a") as fh2:
    print("--- The following had no matches", file=fh2)
    for item in no_match:
        print(item, file=fh2)
    # end for
    fh2.close()
# end with

# =============================================================================
# Write out Codons
# =============================================================================
SeqIO.write(results, 
            CODONS, 
            "fasta")

# =============================================================================
# Write out 5UTR seqs
# =============================================================================
SeqIO.write(results_UTR5, 
            UTR5, 
            "fasta")

# =============================================================================
# Write out 3UTR seqs
# =============================================================================
SeqIO.write(results_UTR3, 
            UTR3, 
            "fasta")


# =============================================================================
# End of file
# =============================================================================
