import os
from Bio import SeqIO

markers = ["ndhF","rpL16","rps16-trnQ","trnC-rpoB","trnD-trnT","trnT-trnL"]  # "rpS16" is forward only
contigpath = "seq/contig-adjust/"
assembledpath = "seq/assembled/"

if not os.path.exists(assembledpath):
	os.makedirs(assembledpath)
for marker in markers:
	seqs = list(SeqIO.parse(contigpath + marker + "_aligned_contigs_adjusted.fas", "fasta"))

	# Assumes contigs are sequencial and the reverse is first
	# loop util no more bases; build new (consensus) seq

	numSeqs = len(seqs)
	j = 1
	consensus = {}

	for n in seqs:
		if n.id.endswith("R"):
			seqLength = len(n)
			newSeq = ""
			for i in range(seqLength):
				rev = n.seq[i]
				forw = seqs[j].seq[i]
				if forw == "A" or rev == "A":
					if forw == "C" or rev == "C": newSeq += "M"
					if forw == "G" or rev == "G": newSeq += "R"
					if forw == "T" or rev == "T": newSeq += "W"
				if forw == "C" or rev == "C":
					if forw == "G" or rev == "G": newSeq += "S"
					if forw == "T" or rev == "T": newSeq += "Y"
				if forw == "G" and rev == "T" or forw == "T" and rev == "G": newSeq += "K"
				# The order below matters: if the N and gap tests both succeed, ensures N
				if forw == "N" and rev != "N": newSeq += rev
				if rev == "N" and forw != "N": newSeq += forw
				if forw == "-" and rev != "-": newSeq += rev
				if rev == "-" and forw != "-": newSeq += forw
				if forw == rev: newSeq += forw  # Same basepair overrides degenerate bps
				# Not implemented:
				# A and C and T = H
				# A and C and G = V
				# A and T and G = D
				# C and T and G = B
			consensus[n.id[:-2]] = newSeq
		j += 1

	# write a file with ">" + consensus.key + "\n" + consensus.value + "\n" <loop>
	f = open(assembledpath + marker + "_assembled.fas", "a")
	for key, value in consensus.items():
		f.write(">" + key + "\n" + value + "\n")
	f.close()
