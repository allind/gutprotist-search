#! /usr/bin/env python
#usage: script.py [bamfile] [buscos]
from pysam import AlignmentFile
from Bio import SeqIO
import sys
def main(argv):

	#parse reference file..?
	#parse samtools header
	bam = AlignmentFile(sys.argv[1])

	observed = []
	for read in bam.fetch():
		if read.reference_name not in observed:
			observed.append(read.reference_name)

	ref_seqs = {}
	for seq in SeqIO.parse(sys.argv[2], 'fasta'):
		if seq.id in observed:
			ref_seqs[seq.id] = str(seq.seq)

	#count coverage
	print("Organism\tGene\tReadcount\tCorrect_bases\tIncorrect_bases\tTotal_bases\tSubjlen\tCoverage\tPercent_ID")
	for o in observed:
		species = "_".join(o.split('_')[1:3])
		contig_counts = bam.count(o, start = 0, end = len(ref_seqs[o]))

		counts = bam.count_coverage(o, start = 0, end = len(ref_seqs[o]))
		pos_ids = []
		trues = 0
		falses = 0
		total = 0
		for ref_pos in range(0, len(ref_seqs[o])):
			total += sum(counts[nt][ref_pos] for nt in range(4))
		
		if total == 0:
			continue

		for ref_pos in range(0, len(ref_seqs[o])):
			ref_allele = ref_seqs[o][ref_pos]
			depth = sum(counts[nt][ref_pos] for nt in range(4))
			count_a = counts[0][ref_pos]
			count_c = counts[1][ref_pos]
			count_g = counts[2][ref_pos]
			count_t = counts[3][ref_pos]
			values = [o, ref_pos + 1, ref_allele, depth, count_a, count_c, count_g, count_t]
			if depth > 0:
				#now we calculate the percentage
				not_n = True
				if ref_allele == "A":
					true = count_a
					false = count_c + count_g + count_t
				elif ref_allele == "C":
					true = count_c
					false = count_a + count_g + count_t
				elif ref_allele == "G":
					true = count_g
					false = count_a + count_c + count_t
				elif ref_allele == "T":
					true = count_t
					false = count_a + count_c + count_g
				else:
					#it's an n, skip it
					not_n = False
				#maybe just have it as an absolute. if there's one mismatch it's all wrong.
				if not_n:
					if false > 0:
						falses += 1
					else:
						trues += 1

		seqlen = len(ref_seqs[o])
		coverage = round(((trues + falses)/seqlen) * 100, 2)
		if trues == 0 and falses == 0:
			pid = 0
		else:
			pid = round((trues/ (trues + falses)) * 100, 2)
		print(species + '\t' + o + '\t' + str(contig_counts) + '\t' 
			+ str(trues) + '\t' 
			+ str(falses) + '\t' 
			+ str(trues + falses) + '\t' 
			+ str(seqlen) + '\t'
			+ str(coverage) + '\t'
			+ str(pid))

if __name__ == "__main__":
  main(sys.argv)
