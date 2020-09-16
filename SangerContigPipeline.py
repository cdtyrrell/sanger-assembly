#!/usr/bin/env python3

import os
import math
from Bio import SeqIO
from Bio import Seq

markers = ["ndhF","rpL16","rpS16","rps16-trnQ","trnC-rpoB","trnD-trnT","trnT-trnL"]
rawpath = "seq/raw/"

def find_overlap(forward, reverse, percent = 0.5, tolerance = 8):
	minLen = min(len(reverse.seq), len(forward.seq))
	overlap = math.ceil(minLen * percent)
	addTo = minLen - overlap
	pairwise = tolerance + 1  # needed to get the while started
	lastMinPairwise = minLen  # sets this very high
	while pairwise > tolerance:
		pairwise = 0
		for i in range(overlap):
			if forward.seq[i] != reverse.seq[i + addTo]: pairwise += 1
		po = pairwise/overlap
		if po < lastMinPairwise and overlap > 2 * tolerance:  # prevent trivial values from overwriting last_min_pairwise by adjusting the tolerance
			lastMinPairwise = po
			lastMinAddTo = addTo
		overlap -= 1
		addTo += 1
		if overlap < 2: break
	if lastMinPairwise < po: offset = lastMinAddTo
	else: offset = addTo - 1
	return offset


def make_seq_array(marker, path = rawpath):
	dir = path + marker + "/"
	seq2fas = []
	dirObject = sorted(os.scandir(dir), key = lambda file: file.name, reverse=True)  # reverse puts the 'R' files first
	for entry in dirObject:
		if entry.path.endswith(".ab1"):
			try:
				seq_file = open(entry.path, "rb")
			except FileNotFoundError:
				print("File not found.")
				exit()

			seq_tmp = SeqIO.read(seq_file, "abi-trim")
			n = seq_tmp.name.split("_")
			if entry.path.lower().endswith("r.ab1"):
				seq_tmp.seq = seq_tmp.seq.reverse_complement()
				seq_tmp.id = n[0] + "_" + n[1] + "_R"
			else: seq_tmp.id = n[0] + "_" + n[1] + "_F"
			seq2fas.append(seq_tmp)
	return seq2fas


for marker in markers:
	seqArray = make_seq_array(marker)
	arrLen = len(seqArray)
	unusedSeq = [True] * arrLen
	if not os.path.exists("seq/contig"):
		os.makedirs("seq/contig")
	for i in range(arrLen):
		if unusedSeq[i]:
			for j in range(i+1, arrLen):
				if seqArray[i].id[:-2] == seqArray[j].id[:-2]:
					unusedSeq[i] = False
					unusedSeq[j] = False
					if seqArray[i].id.endswith("R"):
						slideDist = find_overlap(seqArray[j], seqArray[i], 1, 8)
						gapSeq = "-" * slideDist
						seqArray[j] = Seq.Seq(gapSeq) + seqArray[j]
					else:
						slideDist = find_overlap(seqArray[i], seqArray[j], 1, 8)
						gapSeq = "-" * slideDist
						seqArray[i] = Seq.Seq(gapSeq) + seqArray[i]
					print(seqArray[i].id + " <-> " + seqArray[j].id + " = {0}".format(slideDist))
					break
		SeqIO.write(seqArray, "seq/contig/" + marker + "_aligned_contigs.fas", "fasta")
