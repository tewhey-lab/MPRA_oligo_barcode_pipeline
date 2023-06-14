###go through fastq file to pull barcodes

import pandas as pd
import os
import sys
import re
import time
import pathlib
import gzip
argv = sys.argv

#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

fastqfile = argv[1]
#dictfile = argv[2]
out_id = argv[2]
read_number = int(argv[3])
bc_len = int(argv[4])

current_path = os.getcwd()

def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return filename

#print("Reading Parsed File...")
#BC_list = pd.read_table(dictfile, header=None, index_col=False, usecols=[0,1])
#print("Setting up Dictionary...")
#BC_dict = dict(zip(BC_list[0], BC_list[1]))

with open("%s/%s.match" % (current_path, out_id), "w") as match_oligo:
    with open("%s/%s.reject.fastq" % (current_path, out_id), "w") as reject_fastq:
        with open("%s/%s.reject.bc" % (current_path, out_id), "w") as reject_bc:
            with gzip.open(fastqfile, "rt") as handle:
                print("Reading Records...")
                for record in SeqIO.parse(handle, "fastq"):
                    seq_only = record.seq
                    if read_number != 2:
                        seq_only = seq_only.reverse_complement()
                    seq_only = str(seq_only)
                    # Grab the first 20 bases of the sequence
                    if read_number == 2:
                        bc_seq = seq_only[0:bc_len]
                    if read_number != 2:
                        bc_seq = seq_only[-bc_len:]

                    # Check for barcode presence in the dictionary
                    # if bc_seq in BC_dict:
                    match_oligo.write("%s\t%s\n" % (record.name,bc_seq))
                    #if bc_seq not in BC_dict:
                    #    SeqIO.write(record, reject_fastq, "fastq")
                    #    reject_bc.write("%s\t%s\n" % (record.name,bc_seq))
