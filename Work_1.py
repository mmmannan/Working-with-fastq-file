# Filtering DNA sequence from sample1.fq
# Create an output file sample1_filtered.fq only with whose reads
# which are longer or equal to 200 and do not contain any unambiguous bases

from Bio import SeqIO

handle = open("sample1.fq")
seq_list = []
for seq_record in SeqIO.parse(handle, "fastq"):
    sequence = seq_record.seq
    sequence = str(sequence)
    seq_list.append(sequence)
handle.close()

Nucleotides = ["A", "T", "C", "G"]
def validate_DNA_seq(dna_seq):
    '''Check the read to make sure it does not contain uncalled bases'''
    tmpseq = dna_seq.upper()
    for base in tmpseq:
        if base not in Nucleotides:
            return False
    return True

filtered_reads = []
for i in range(len(seq_list)):
    if validate_DNA_seq(seq_list[i]) == True:
        if len(seq_list[i]) >= 200:
            read = seq_list[i]
            filtered_reads.append(read)

with open("sample1_filtered.fq", "w+") as filtered_file:
    for line in filtered_reads:
        filtered_file.write(F"{line}\n")



