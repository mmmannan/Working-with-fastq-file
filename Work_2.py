# Calculating fastq stats
# Calculate id, GC content, A count, T count, C count, G count, average quality for each reads
# write those parameters for each reads as columns to a csv file

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd

handle = open("reads.fq")
seq_list = []
quality_scores_list = []
for seq_record in SeqIO.parse(handle, "fastq"):
    sequence = seq_record.seq
    sequence = str(sequence)
    seq_list.append(sequence)
    quality_scores = seq_record.letter_annotations["phred_quality"]
    quality_scores_list.append(quality_scores)
handle.close()

def average(List):
    '''Average of a list of numbers'''
    sum = 0
    for i in range(len(List)):
        sum = sum + List[i]
        avg = sum / len(List)
    return avg


read_id = []
length_bp = []
gc_content = []
A_count = []
T_count = []
C_count = []
G_count = []
average_quality = []

for i in range(len(seq_list)):
    length = len(seq_list[i])
    gc = round(gc_fraction(seq_list[i]) * 100, 2)
    A = seq_list[i].count("A")
    T = seq_list[i].count("T")
    C = seq_list[i].count("C")
    G = seq_list[i].count("G")
    avg_quality = round(average(quality_scores_list[i]), 2)
    read_id.append(f"read {i + 1}")
    length_bp.append(length)
    gc_content.append(gc)
    A_count.append(A)
    T_count.append(T)
    C_count.append(C)
    G_count.append(G)
    average_quality.append(avg_quality)


dict = {"read ID" : read_id, "length (bp)": length_bp, "GC content (%)": gc_content,
        "A count" : A_count, "T count" : T_count, "C count" : C_count,
        "G count" : G_count, "Average Quality" : average_quality}

df = pd.DataFrame(dict)

df.to_csv('reads_fq.csv', header = True, index = False)

