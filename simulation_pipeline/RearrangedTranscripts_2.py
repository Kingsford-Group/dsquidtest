#!/bin/python

'''
This script retrieves the transcript sequences post-rearrangements.

genome rearrangement simulation: RSVsim
used outputs: deletions.csv inversions.csv translocations.csv

The corresponding effect on transcripts:
	- deletion.csv: all exon regions in the deleted region are also deleted from all transcripts:
	- inversion.csv: assume the after the inversion, the exon region is the same but reverse complemented. 
		* inversion with a gene: all transcripts with an inverted region has a reverse complated sequence corresponding to that region.
		* inversion spanning multiple genes: substitute the within inversion region with the within inversion region of the genes at the other end
		* inversion one side hits genes: the within inversion region is lost
	- translocation.csv (same as inversion): 
		* both breakpoint hit genes: directly connect the transcript regions around each breakpoint
		* one breakpoint hits genes: transcript terminates at the breakpoint
'''

import sys
import copy
import collections
import numpy as np
import pandas as pd
from SpecialTranscriptClass import *

DEL = collections.namedtuple("DEL", ["id", "chr", "startpos", "endpos"])
INV = collections.namedtuple("INV", ["id", "chr", "startpos", "endpos"])
TRA = collections.namedtuple("TRA", ["id", "chr1", "startpos1", "endpos1", "invert1", "chr2", "startpos2", "endpos2", "invert2"])
ADJACENCY = collections.namedtuple("ADJACENCY", ["var_id", "chr1", "startpos1", "endpos1", "isleft1", "chr2", "startpos2", "endpos2", "isleft2"])


def FindStopCodon(string):
	pos = -1
	if "TAG" in string.upper():
		pos = string.upper().index("TAG")
	if "TAA" in string.upper():
		tmppos = string.upper().index("TAA")
		if pos >= 0 and tmppos < pos:
			pos = tmppos
	if "TGA" in string.upper():
		tmppos = string.upper().index("TGA")
		if pos >= 0 and tmppos < pos:
			pos = tmppos
	return pos


def ReadDeletions(filename):
	deletions = []
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
#			assert(line.strip() == "Name\tChr\tStart\tEnd\tSize\tBpSeq")
			continue
		strs = line.strip().split("\t")
		id = strs[0]
		chr = strs[1]
		startpos = int(strs[2]) - 1
		endpos = int(strs[3])
		# add to deletions
		deletions.append( DEL(id, chr, startpos, endpos) )
	fp.close()
	return deletions


def ReadInversions(filename):
	inversions = []
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
#			assert(line.strip() == "Name\tChr\tStart\tEnd\tSize\tBpSeq_3prime\tBpSeq_5prime")
			continue
		strs = line.strip().split("\t")
		id = strs[0]
		chr = strs[1]
		startpos = int(strs[2]) - 1
		endpos = int(strs[3])
		# add to inversions
		inversions.append( INV(id, chr, startpos, endpos) )
	fp.close()
	return inversions


def ReadTranslocations_insertion(filename, ChrLengths):
	deletions = []
	translocations = []
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
#			assert(line[:9] == "Name\tChrA")
			continue
		strs = line.strip().split("\t")
		deletions.append( DEL(strs[0]+"_0", strs[1], int(strs[2])-1, int(strs[3])) )
		# the first translocation
		chr1 = strs[4]
		startpos1 = 0
		endpos1 = int(strs[5]) - 1
		chr2 = strs[1]
		startpos2 = int(strs[2]) - 1
		endpos2 = int(strs[3])
		if (chr1 == chr2) and (endpos2 < endpos1):
			startpos1 = endpos2
		translocations.append( TRA(strs[0]+"_1", chr1, startpos1, endpos1, False, chr2, startpos2, endpos2, False) )
		# the second translocation
		chr1 = strs[1]
		startpos1 = int(strs[2]) - 1
		endpos1 = int(strs[3])
		chr2 = strs[4]
		startpos2 = int(strs[5])
		endpos2 = ChrLengths[chr2]
		if (chr1 == chr2) and (startpos1 > startpos2):
			endpos2 = startpos1
		translocations.append( TRA(strs[0]+"_2", chr1, startpos1, endpos1, False, chr2, startpos2, endpos2, False) )
	fp.close()
	return deletions, translocations


def ReadTranslocations(filename, ChrLengths):
	translocations = []
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
#			assert(line[:9] == "Name\tChrA")
			continue
		strs = line.strip().split("\t")
		#assert(strs[9] == "TRUE")
		id = strs[0]
		# a translocation in RSVsim record indicate 2 sequences
		# 5' of one chromosome is exchanged with 5' of another chromosome
		if int(strs[2]) == 1 and int(strs[6]) == 1:
			# the first sequence: beginning of the second + end of the first chromosome in record
			chr1 = strs[5]
			startpos1 = int(strs[6]) - 1
			endpos1 = int(strs[7])
			invert1 = False
			chr2 = strs[1]
			startpos2 = int(strs[2])
			endpos2 = ChrLengths[chr2]
			invert2 = False
			translocations.append( TRA(id+"_1", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
			# the second sequence: beginning of the first + end of the second chromosome in record
			chr1 = strs[1]
			startpos1 = int(strs[2]) - 1
			endpos1 = int(strs[3])
			invert1 = False
			chr2 = strs[5]
			startpos2 = int(strs[7])
			endpos2 = ChrLengths[chr2]
			invert2 = False
			translocations.append( TRA(id+"_2", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
		# 5' of one chromosome is exchanged with 3; of another chromosome in a reverse complement way
		elif int(strs[2]) == 1 and int(strs[6]) != 1:
			# the first sequence: reverse complement end of the second + end of the first
			chr1 = strs[5]
			startpos1 = int(strs[6]) - 1
			endpos1 = int(strs[7])
			assert(endpos1 == ChrLengths[chr1])
			invert1 = True
			chr2 = strs[1]
			startpos2 = int(strs[2])
			endpos2 = ChrLengths[chr2]
			invert2 = False
			translocations.append( TRA(id+"_1", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
			# the second sequence: beginning of the second + reverse complement of the beginning of the first
			chr1 = strs[5]
			startpos1 = 0
			endpos1 = int(strs[6])
			invert1 = False
			chr2 = strs[1]
			startpos2 = 0
			endpos2 = int(strs[2])
			invert2 = True
			translocations.append( TRA(id+"_2", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
		# 3' of one chromosome is exchanged with 5' of another chromosome in a reverse complement way
		elif int(strs[2]) != 1 and int(strs[6]) == 1:
			# the first seequence: the beginning of the first + reverse complement of the beginning of the second
			chr1 = strs[1]
			startpos1 = 0
			endpos1 = int(strs[2])
			invert1 = False
			chr2 = strs[5]
			startpos2 = 0
			endpos2 = int(strs[7])
			invert2 = True
			translocations.append( TRA(id+"_1", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
			# the second sequence: reverse complement of the end of the first + the end of the second
			chr1 = strs[1]
			startpos1 = int(strs[2]) - 1
			endpos1 = int(strs[3])
			assert(endpos1 == ChrLengths[chr1])
			invert1 = True
			chr2 = strs[5]
			startpos2 = int(strs[7])
			endpos2 = ChrLengths[chr2]
			invert2 = False
			translocations.append( TRA(id+"_2", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
		# 3' of one chromosome is exchanged with 3; of another chromosome
		elif int(strs[2]) != 1 and int(strs[6]) != 1:
			# the first sequence: the beginning of the first + the end of the second
			chr1 = strs[1]
			startpos1 = 0
			endpos1 = int(strs[2])
			invert1 = False
			chr2 = strs[5]
			startpos2 = int(strs[6]) - 1
			endpos2 = int(strs[7])
			assert(endpos2 == ChrLengths[chr2])
			invert2 = False
			translocations.append( TRA(id+"_1", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
			# the second sequence: the beginning of the second + the end of the first
			chr1 = strs[5]
			startpos1 = 0
			endpos1 = int(strs[6])
			invert1 = False
			chr2 = strs[1]
			startpos2 = int(strs[2]) - 1
			endpos2 = int(strs[3])
			assert(endpos2 == ChrLengths[chr2])
			invert2 = False
			translocations.append( TRA(id+"_2", chr1, startpos1, endpos1, invert1, chr2, startpos2, endpos2, invert2) )
	fp.close()
	return translocations


def UpdateAnnotation_del(transcripts, deletions, tlocator, genome):
	# overwrite on the original transcript names
	for d in deletions:
		tids_b1 = [tname for tname in tlocator.BinarySearchPosition(d.chr, d.startpos) if tname in transcripts]
		tids_b2 = [tname for tname in tlocator.BinarySearchPosition(d.chr, d.endpos) if tname in transcripts]
		for tname in set(tids_b1 + tids_b2):
			# the deletion has the same effect as alternative splicing or truncation, so novel adjacency is not recorded.
			# but rename with deletion, and avoid to use inversion or translocation to further modify the transcripts.
			t = transcripts[tname]
			assert((t.StartPos < d.startpos and t.EndPos > d.startpos) or (t.StartPos < d.endpos and t.EndPos > d.endpos))
			newt = copy.deepcopy(t)
			newt.TransID = t.TransID + ".del"
			newt.Exons = []
			for e in t.Exons:
				# exon and deleted region does not overlap
				if not (e.chr == d.chr and max(e.startpos, d.startpos) < min(e.endpos, d.endpos)):
					newt.Exons.append(e)
				# else, the exon is partially overlap
				elif e.startpos < d.startpos:
					newt.Exons.append( ITVL(e.chr, e.startpos, d.startpos, e.strand) )
				# another type of partial overlap
				elif e.endpos > d.endpos:
					newt.Exons.append( ITVL(e.chr, d.endpos, e.endpos, e.strand) )
			# update transcripts
			newt.StartPos = min([e.startpos for e in newt.Exons])
			newt.EndPos = max([e.endpos for e in newt.Exons])
			transcripts[newt.TransID] = newt
			transcripts.pop(tname, None)
	return transcripts


def UpdateAnnotation_inv(transcripts, inversions, tlocator, genome):
	alteration_indicator = {tname:False for tname in transcripts.keys()}
	noveladjacencies = []
	corresponding_trans = []
	for v in inversions:
		tids_start = [tname for tname in tlocator.BinarySearchPosition(v.chr, v.startpos) if tname in transcripts]
		tids_end = [tname for tname in tlocator.BinarySearchPosition(v.chr, v.endpos) if tname in transcripts]
		assert( np.all([tname in transcripts for tname in tids_start+tids_end]) )
		# check whether a single gene is spanning the whole inversion region
		tids_coverwhole = list(set(tids_start)&set(tids_end))
		tids_start = list(set(tids_start) - set(tids_coverwhole))
		tids_end = list(set(tids_end) - set(tids_coverwhole))
		# for whole spanning transcripts
		for tname in tids_coverwhole:
			t = transcripts[tname]
			assert(t.Chr == v.chr)
			assert(t.StartPos < v.startpos and t.EndPos > v.endpos)
			newt = copy.deepcopy(t)
			newt.TransID = tname+".inv"
			# separate the exons to the left, in the middle, and to the right of the inversion
			left_exons = []
			middle_exons = []
			right_exons = []
			for e in t.Exons:
				# if e overlap with the inversion region
				if max(e.startpos, v.startpos) < min(e.endpos, v.endpos):
					middle_exons.append( ITVL(e.chr, max(e.startpos, v.startpos), min(e.endpos, v.endpos), not e.strand) )
				# if e overlap with the left side of the inversion region
				if e.startpos < v.startpos:
					left_exons.append( ITVL(e.chr, e.startpos, min(e.endpos, v.startpos), e.strand) )
				# if e overlap with the right side of the inversion region
				if e.endpos > v.endpos:
					right_exons.append( ITVL(e.chr, max(e.startpos, v.endpos), e.endpos, e.strand) )
			# sort the exons to the genomic coordinate
			assert(len(left_exons) > 0 and len(right_exons) > 0)
			# if the inversion is totally with intron, the transcript is assumed to be unchanged
			if len(middle_exons) == 0:
				continue
			# update exons in the new transcript
			middle_exons = middle_exons[::-1]
			if left_exons[0].strand:
				newt.Exons = left_exons + middle_exons + right_exons
			else:
				newt.Exons = right_exons + middle_exons + left_exons
			# update transcripts
			transcripts[newt.TransID] = newt
			alteration_indicator[tname] = True
			# update novel adjacency
			if left_exons[0].strand:
				noveladjacencies.append( ADJACENCY(v.id, left_exons[-1].chr, left_exons[-1].startpos, left_exons[-1].endpos, False, middle_exons[0].chr, middle_exons[0].startpos, middle_exons[0].endpos, False) )
				noveladjacencies.append( ADJACENCY(v.id, middle_exons[-1].chr, middle_exons[-1].startpos, middle_exons[-1].endpos, True, right_exons[0].chr, right_exons[0].startpos, right_exons[0].endpos, True) )
			else:
				noveladjacencies.append( ADJACENCY(v.id, right_exons[-1].chr, right_exons[-1].startpos, right_exons[-1].endpos, True, middle_exons[0].chr, middle_exons[0].startpos, middle_exons[0].endpos, True) )
				noveladjacencies.append( ADJACENCY(v.id, middle_exons[-1].chr, middle_exons[-1].startpos, middle_exons[-1].endpos, False, left_exons[0].chr, left_exons[0].startpos, left_exons[0].endpos, False) )
			corresponding_trans += [newt.TransID, newt.TransID]
		assert( len(corresponding_trans) == len(noveladjacencies) )
		# for a pair of transcripts occupying separate ends: combine two transcripts if they are on different strand
		for tname1 in tids_start:
			for tname2 in tids_end:
				t1 = transcripts[tname1]
				t2 = transcripts[tname2]
				assert(t1.StartPos < v.startpos and t1.EndPos > v.startpos and t2.StartPos < v.endpos and t2.EndPos > v.endpos)
				# check whether they are on different strand
				if t1.Strand == t2.Strand:
					continue
				# collect the exons of t1 / t2 that are outside inversion and inside inversion
				out_exons_1 = []
				in_exons_1 = []
				out_exons_2 = []
				in_exons_2 = []
				for e in t1.Exons:
					if not (e.chr == v.chr and max(e.startpos, v.startpos) < min(e.endpos, v.endpos)):
						out_exons_1.append( copy.deepcopy(e) )
					elif e.chr == v.chr and e.startpos >= v.startpos and e.endpos <= v.endpos:
						in_exons_1.append( copy.deepcopy(e) )
					else:
						# assert that the exon must properly overlap with the inversion
						assert( (e.startpos < v.startpos and e.endpos < v.endpos) )
						out_exons_1.append( ITVL(e.chr, e.startpos, v.startpos, e.strand) )
						in_exons_1.append( ITVL(e.chr, v.startpos, e.endpos, e.strand) )
				for e in t2.Exons:
					if not (e.chr == v.chr and max(e.startpos, v.startpos) < min(e.endpos, v.endpos)):
						out_exons_2.append( copy.deepcopy(e) )
					elif e.chr == v.chr and e.startpos >= v.startpos and e.endpos <= v.endpos:
						in_exons_2.append( copy.deepcopy(e) )
					else:
						assert( (e.endpos > v.endpos and e.startpos > v.startpos) )
						out_exons_2.append( ITVL(e.chr, v.endpos, e.endpos, e.strand) )
						in_exons_2.append( ITVL(e.chr, e.startpos, v.endpos, e.strand) )
				assert( len(out_exons_1) > 0 and len(in_exons_1) > 0 and len(out_exons_2) > 0 and len(in_exons_2) > 0 )
				assert( np.all([e.strand == out_exons_1[0].strand for e in out_exons_1]) )
				assert( np.all([e.strand == in_exons_1[0].strand for e in in_exons_1]) )
				assert( np.all([e.strand == out_exons_2[0].strand for e in out_exons_2]) )
				assert( np.all([e.strand == in_exons_2[0].strand for e in in_exons_2]) )
				# the left breakpoint
				assert( np.all([e.strand != out_exons_1[0].strand for e in in_exons_2]) )
				newt = Transcript_t(tname1+"|"+tname2+".inv", t1.GeneID+"|"+t2.GeneID, v.chr, out_exons_1[0].strand, out_exons_1[0].startpos, in_exons_2[0].endpos)
				if out_exons_1[0].strand:
					newt.Exons = out_exons_1 + in_exons_2
				else:
					newt.Exons = in_exons_2 + out_exons_1
				transcripts[newt.TransID] = newt
				# update novel adjacency
				if out_exons_1[0].strand:
					noveladjacencies.append( ADJACENCY(v.id, out_exons_1[-1].chr, out_exons_1[-1].startpos, out_exons_1[-1].endpos, False, in_exons_2[0].chr, in_exons_2[0].startpos, in_exons_2[0].endpos, False) )
				else:
					noveladjacencies.append( ADJACENCY(v.id, in_exons_2[-1].chr, in_exons_2[-1].startpos, in_exons_2[-1].endpos, False, out_exons_1[0].chr, out_exons_1[0].startpos, out_exons_1[0].endpos, False) )
				corresponding_trans.append(newt.TransID)
				# the right breakpoint
				newt = Transcript_t(tname1+".inv|"+tname2, t1.GeneID+"|"+t2.GeneID, v.chr, in_exons_1[0].strand, in_exons_1[-1].startpos, out_exons_2[-1].endpos)
				if out_exons_2[0].strand:
					newt.Exons = in_exons_1 + out_exons_2
				else:
					newt.Exons = out_exons_2 + in_exons_1
				transcripts[newt.TransID] = newt
				# update novel adjacency
				if out_exons_2[0].strand:
					noveladjacencies.append( ADJACENCY(v.id, in_exons_1[-1].chr, in_exons_1[-1].startpos, in_exons_1[-1].endpos, True, out_exons_2[0].chr, out_exons_2[0].startpos, out_exons_2[0].endpos, True) )
				else:
					noveladjacencies.append( ADJACENCY(v.id, out_exons_2[-1].chr, out_exons_2[-1].startpos, out_exons_2[-1].endpos, True, in_exons_1[0].chr, in_exons_1[0].startpos, in_exons_1[0].endpos, True) )
				corresponding_trans.append(newt.TransID)
				# update alteration indicator
				alteration_indicator[tname1] = True
				alteration_indicator[tname2] = True
			assert( len(corresponding_trans) == len(noveladjacencies) )
			# if all transcripts are on the same strand, do the search for stop codon to simulate new sequence
			for tname in tids_start + tids_end:
				if alteration_indicator[tname]:
					continue
				t = transcripts[tname]
				assert((t.StartPos < v.startpos and t.EndPos > v.startpos) or (t.StartPos < v.endpos and t.EndPos > v.endpos))
				newt = copy.deepcopy(t)
				newt.TransID = tname + ".inv"
				newt.Exons = []
				# check whether the TSS is within inversion or outside
				tss = t.Exons[0].startpos if t.Strand else t.Exons[0].endpos
				assert( tss != v.startpos and tss != v.endpos)
				if tss < v.startpos: # the transcript must be on forward strand
					assert(t.Strand)
					for e in t.Exons:
						if e.startpos < v.startpos:
							newt.Exons.append( ITVL(e.chr, e.startpos, min(e.endpos,v.startpos), e.strand) )
					assert(len(newt.Exons) > 0)
					# the next line are commented out, that is, not extend the last unchanged exon to the inversion breakpoint
					# newt.Exons[-1] = ITVL(newt.Exons[-1].chr, newt.Exons[-1].startpos, v.startpos, newt.Exons[-1].strand)
					# extend the last exon over the inverted region to find the first stop codon. TSS is forward strand, within-inversion stop codon is reverse strand
					position_stopdodon = FindStopCodon( ReverseComplement(genome[v.chr][(v.startpos):(v.endpos)]) )
					if position_stopdodon == -1:
						newt.Exons.append( ITVL(v.chr, v.startpos, v.endpos, not t.Strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, False, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, False) )
						corresponding_trans.append(newt.TransID)
						# extend to the right of the inversion
						position_stopdodon = FindStopCodon( genome[v.chr][(v.endpos):] )
						assert(position_stopdodon > 0)
						newt.Exons.append( ITVL(v.chr, v.endpos, v.endpos + position_stopdodon + 3, t.Strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, True, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, True) )
						corresponding_trans.append(newt.TransID)
					else:
						newt.Exons.append( ITVL(v.chr, v.endpos - position_stopdodon - 3, v.endpos, not t.Strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, False, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, False) )
						corresponding_trans.append(newt.TransID)
				elif tss > v.startpos and tss < v.endpos:
					for e in t.Exons:
						# select the region within the inversion
						if max(e.startpos, v.startpos) < min(e.endpos, v.endpos):
							newt.Exons.append( ITVL(e.chr, max(e.startpos, v.startpos), min(e.endpos, v.endpos), e.strand) )
					assert(len(newt.Exons) > 0)
					# the next 4 lines are commented out, that is, not extend the last unchanged exon to the inversion breakpoint
					# if newt.Exons[0].strand:
					# 	newt.Exons[-1] = ITVL(newt.Exons[-1].chr, newt.Exons[-1].startpos, v.endpos, newt.Exons[-1].strand)
					# else:
					# 	newt.Exons[-1] = ITVL(newt.Exons[-1].chr, v.startpos, newt.Exons[-1].endpos, newt.Exons[-1].strand)
					# extend to find stop codon
					if newt.Exons[0].strand: # after inversion, TSS becomes reverse strand and extend to the sequence before v.startpos
						position_stopdodon = FindStopCodon( ReverseComplement(genome[v.chr][:(v.startpos)]) )
						if position_stopdodon >= 0:
							newt.Exons.append( ITVL(v.chr, v.startpos - position_stopdodon - 3, v.startpos, not newt.Exons[0].strand) )
						else:
							newt.Exons.append( ITVL(v.chr, 0, v.startpos, not newt.Exons[0].strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, False, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, False) )
						corresponding_trans.append(newt.TransID)
					else: # after the inversion, TSS become forward strand, and extend to the sequence after v.endpos
						position_stopdodon = FindStopCodon( genome[v.chr][(v.endpos):] )
						if position_stopdodon >= 0:
							newt.Exons.append( ITVL(v.chr, v.endpos, v.endpos + position_stopdodon + 3, not newt.Exons[0].strand) )
						else:
							newt.Exons.append( ITVL(v.chr, v.endpos, len(genome[v.chr]), not newt.Exons[0].strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, True, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, True) )
						corresponding_trans.append(newt.TransID)
				elif tss > v.endpos:
					assert(not t.Strand)
					for e in t.Exons:
						if e.endpos > v.endpos:
							newt.Exons.append( ITVL(e.chr, max(v.endpos, e.startpos), e.endpos, e.strand) )
					assert(len(newt.Exons) > 0)
					# the next line are commented out, that is, not extend the last unchanged exon to the inversion breakpoint
					# newt.Exons[-1] = ITVL(newt.Exons[-1].chr, v.endpos, newt.Exons[-1].endpos, newt.Exons[-1].strand)
					# extend the last exon over the inverted region to find the first stop codon. TSS is reverse strand, within-inversion stop codon is forward strand
					position_stopdodon = FindStopCodon( genome[v.chr][(v.startpos):(v.endpos)] )
					if position_stopdodon == -1:
						newt.Exons.append( ITVL(v.chr, v.startpos, v.endpos, not t.Strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, True, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, True) )
						corresponding_trans.append(newt.TransID)
						position_stopdodon = FindStopCodon( ReverseComplement(genome[v.chr][:(v.startpos)]) )
						assert(position_stopdodon >= 0)
						newt.Exons.append( ITVL(v.chr, v.startpos - position_stopdodon - 3, v.startpos, t.Strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, False, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, False) )
						corresponding_trans.append(newt.TransID)
					else:
						newt.Exons.append( ITVL(v.chr, v.startpos, v.startpos + position_stopdodon + 3, not t.Strand) )
						# update novel adjacency
						noveladjacencies.append( ADJACENCY(v.id, v.chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, True, v.chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, True) )
						corresponding_trans.append(newt.TransID)
				# update transcript
				transcripts[newt.TransID] = newt
				# update alteration indicator
				alteration_indicator[tname] = True
		assert( len(corresponding_trans) == len(noveladjacencies) )
	# pop up the altered transcripts
	for tname in [tname for tname,v in alteration_indicator.items() if v]:
		transcripts.pop(tname, None)
	return transcripts, noveladjacencies, corresponding_trans


def UpdateAnnotation_tra(transcripts, translocations, tlocator, genome):
	alteration_indicator = {tname:False for tname in transcripts.keys()}
	noveladjacencies = []
	corresponding_trans = []
	for a in translocations:
		if not a.invert1:
			tids_b1 = [tname for tname in tlocator.BinarySearchPosition(a.chr1, a.endpos1) if tname in transcripts]
		else:
			tids_b1 = [tname for tname in tlocator.BinarySearchPosition(a.chr1, a.startpos1) if tname in transcripts]
		if not a.invert2:
			if a.chr2 == "chr5" and a.startpos2 == 181435605:
			    print(a)
			tids_b2 = [tname for tname in tlocator.BinarySearchPosition(a.chr2, a.startpos2) if tname in transcripts]
		else:
			tids_b2 = [tname for tname in tlocator.BinarySearchPosition(a.chr2, a.endpos2) if tname in transcripts]
		assert( np.all([tname in transcripts for tname in tids_b1+tids_b2]) )
		# separate the transcripts spanning whole vs occupying single breakpoint
		tids_coverwhole = set(tids_b1) & set(tids_b2)
		tids_b1 = list(set(tids_b1) - set(tids_coverwhole))
		tids_b2 = list(set(tids_b2) - set(tids_coverwhole))
		for tname in tids_coverwhole:
			# the translocation is within one single chromosome, must be one of the insertion event
			# re-order the exons based on the translocation
			t = transcripts[tname]
			assert(t.Chr == a.chr1 and t.Chr == a.chr2)
			assert(not a.invert1 and not a.invert2)
			if not (t.StartPos < a.endpos1 and t.EndPos > a.endpos1 and t.StartPos < a.startpos2 and t.EndPos > a.startpos2):
				print(tname)
				print(a)
				print(t.Exons)
			assert(t.StartPos < a.endpos1 and t.EndPos > a.endpos1 and t.StartPos < a.startpos2 and t.EndPos > a.startpos2)
			newt = copy.deepcopy(t)
			newt.TransID = tname+".tra"
			# get the moved region and insertion point
			if a.id[-2:] == "_1":
				assert( not a.invert1 )
				moved_region = (a.startpos2, a.endpos2)
				insertion_point = a.endpos1
			else:
				assert( not a.invert2 )
				moved_region = (a.startpos1, a.endpos1)
				insertion_point = a.startpos2
			# get the exons inside the moved region, and outside + to the left / right of the insertion point
			inside_exons = []
			outside_left = []
			outside_right = []
			for e in t.Exons:
				# add the overlap with the moved region
				if max(e.startpos, moved_region[0]) < min(e.endpos, moved_region[1]):
					inside_exons.append( ITVL(e.chr, max(e.startpos, moved_region[0]), min(e.endpos, moved_region[1]), e.strand) )
				# if there are region outside the moved region, add them to either outside_left or outside_right
				if e.startpos < moved_region[0]:
					subexon = ITVL(e.chr, e.startpos, min(e.endpos, moved_region[0]), e.strand)
					if subexon.startpos < insertion_point:
						outside_left.append( ITVL(subexon.chr, subexon.startpos, min(subexon.endpos, insertion_point), subexon.strand) )
					if subexon.endpos > insertion_point:
						outside_right.append( ITVL(subexon.chr, max(subexon.startpos, insertion_point), subexon.endpos, subexon.strand) )
				# if there are region outside the moved region, add them to either outside_left or outside_right
				if e.endpos > moved_region[1]:
					subexon = ITVL(e.chr, max(e.startpos, moved_region[1]), e.endpos, e.strand)
					if subexon.startpos < insertion_point:
						outside_left.append( ITVL(subexon.chr, subexon.startpos, min(subexon.endpos, insertion_point), subexon.strand) )
					if subexon.endpos > insertion_point:
						outside_right.append( ITVL(subexon.chr, max(subexon.startpos, insertion_point), subexon.endpos, subexon.strand) )
			assert(len(outside_left) > 0 or len(outside_right) > 0)
			# if the moved region is totally within the intron, nothing is changed
			if len(inside_exons) == 0:
				continue
			# sort the exons
			if a.invert1 != a.invert2:
				inside_exons = [ITVL(e.chr, e.startpos, e.endpos, not e.strand) for e in inside_exons[::-1]]
			# merge the exons
			if (len(outside_left) > 0 and outside_left[0].strand) or (len(outside_right) > 0 and outside_right[0].strand):
				newt.Exons = outside_left + inside_exons + outside_right
			else:
				newt.Exons = outside_right + inside_exons + outside_left
			# update transcripts
			transcripts[newt.TransID] = newt
			alteration_indicator[tname] = True
			# update novel adjacency
			if (len(outside_left) > 0 and outside_left[0].strand) or (len(outside_right) > 0 and outside_right[0].strand):
				if len(outside_left) > 0:
					noveladjacencies.append( ADJACENCY(a.id, outside_left[-1].chr, outside_left[-1].startpos, outside_left[-1].endpos, not outside_left[-1].strand, inside_exons[0].chr, inside_exons[0].startpos, inside_exons[0].endpos, inside_exons[0].strand) )
					corresponding_trans += [newt.TransID]
				if len(outside_right) > 0:
					noveladjacencies.append( ADJACENCY(a.id, inside_exons[-1].chr, inside_exons[-1].startpos, inside_exons[-1].endpos, not inside_exons[-1].strand, outside_right[0].chr, outside_right[0].startpos, outside_right[0].endpos, outside_right[0].strand) )
					corresponding_trans += [newt.TransID]
			else:
				if len(outside_right) > 0:
					noveladjacencies.append( ADJACENCY(a.id, outside_right[-1].chr, outside_right[-1].startpos, outside_right[-1].endpos, not outside_right[-1].strand, inside_exons[0].chr, inside_exons[0].startpos, inside_exons[0].endpos, inside_exons[0].strand) )
					corresponding_trans += [newt.TransID]
				if len(outside_left) > 0:
					noveladjacencies.append( ADJACENCY(a.id, inside_exons[-1].chr, inside_exons[-1].startpos, inside_exons[-1].endpos, not inside_exons[-1].strand, outside_left[0].chr, outside_left[0].startpos, outside_left[0].endpos, outside_left[0].strand) )
					corresponding_trans += [newt.TransID]
		assert( len(corresponding_trans) == len(noveladjacencies) )
		# for a pair of transcripts occupying separate ends: combine two transcripts if their strand are the same after a.invert1 and a.invert2
		for tname1 in tids_b1:
			for tname2 in tids_b2:
				t1 = transcripts[tname1]
				t2 = transcripts[tname2]
				if (a.invert1 == a.invert2 and t1.Strand != t2.Strand) or (a.invert1 != a.invert2 and t1.Strand == t2.Strand):
					continue
				part1_exons = []
				part2_exons = []
				for e in t1.Exons:
					if max(e.startpos, a.startpos1) < min(e.endpos, a.endpos1):
						part1_exons.append( ITVL(e.chr, max(e.startpos, a.startpos1), min(e.endpos, a.endpos1), e.strand) )
				for e in t2.Exons:
					if max(e.startpos, a.startpos2) < min(e.endpos, a.endpos2):
						part2_exons.append( ITVL(e.chr, max(e.startpos, a.startpos2), min(e.endpos, a.endpos2), e.strand) )
				# In the insertion case, the inserted region can be too short to be within an intron. Ignore the case when one part of exons is nothing
				if len(part1_exons) == 0 or len(part2_exons) == 0:
					continue
				# sort and merge
				# part1_exons.sort(key = lambda x:x.startpos)
				# if a.invert1:
				# 	part1_exons = part1_exons[::-1]
				# strand_correction = part1_exons[0].strand
				# if a.invert1 != a.invert2:
				# 	strand_correction = not strand_correction
				# part2_exons = [ITVL(e.chr, e.startpos, e.endpos, strand_correction) for e in part2_exons]
				# part2_exons.sort(key = lambda x:x.startpos)
				# if a.invert2:
				# 	part2_exons = part2_exons[::-1]
				strand_correction = part1_exons[0].strand if (not a.invert1) else (not part1_exons[0].strand)
				if a.invert1 != part1_exons[0].strand:
					merged_exons = part1_exons + part2_exons
				else:
					merged_exons = part2_exons + part1_exons
				# update transcripts
				newt = Transcript_t(tname1+"|"+tname2+"."+a.id, t1.GeneID+"|"+t2.GeneID, t1.Chr, strand_correction, -1, -1)
				newt.Exons = merged_exons
				transcripts[newt.TransID] = newt
				# add to novel adjacency
				if a.invert1 != part1_exons[0].strand:
					noveladjacencies.append( ADJACENCY(a.id, part1_exons[-1].chr, part1_exons[-1].startpos, part1_exons[-1].endpos, not part1_exons[-1].strand, part2_exons[0].chr, part2_exons[0].startpos, part2_exons[0].endpos, part2_exons[0].strand) )
				else:
					noveladjacencies.append( ADJACENCY(a.id, part1_exons[0].chr, part1_exons[0].startpos, part1_exons[0].endpos, part1_exons[0].strand, part2_exons[-1].chr, part2_exons[-1].startpos, part2_exons[-1].endpos, not part2_exons[-1].strand) )
				corresponding_trans.append(newt.TransID)
				# update alteration indicator
				alteration_indicator[tname1] = True
				alteration_indicator[tname2] = True
		assert( len(corresponding_trans) == len(noveladjacencies) )
		# if there is transcript that cannot form consistent strand with another, extend it over the breakpoint and find the first stop codon
		for tname in tids_b1 + tids_b2:
			if alteration_indicator[tname]:
				continue
			t = transcripts[tname]
			newt = copy.deepcopy(t)
			newt.TransID = tname + ".tra"
			newt.Exons = []
			# find the tss to see whether it is within either of the connected region
			tss = t.Exons[0].startpos if t.Strand else t.Exons[0].endpos
			# strand and invert1 and invert2 need to satisfy some relation, then it can be extended to the translocated region
			if t.Chr == a.chr1 and a.startpos1 < tss and a.endpos1 > tss and ((not a.invert1) == t.Strand):
				for e in t.Exons:
					if e.chr == a.chr1 and max(e.startpos, a.startpos1) < min(e.endpos, a.endpos1):
						newt.Exons.append( ITVL(e.chr, max(e.startpos, a.startpos1), min(e.endpos, a.endpos1), e.strand) )
				assert(len(newt.Exons) > 0)
				# the next 4 lines are commented out, that is, not extend the last unchanged exon to the translocation breakpoint
				# if t.Strand:
				# 	newt.Exons[-1] = ITVL(newt.Exons[-1].chr, newt.Exons[-1].startpos, a.endpos1, newt.Exons[-1].strand)
				# else:
				# 	newt.Exons[-1] = ITVL(newt.Exons[-1].chr, a.startpos1, newt.Exons[-1].endpos, newt.Exons[-1].strand)
				position_stopdodon = FindStopCodon(genome[a.chr2][(a.startpos2):(a.endpos2)]) if not a.invert2 else FindStopCodon( ReverseComplement(genome[a.chr2][(a.startpos2):(a.endpos2)]) )
				if position_stopdodon == -1:
					if a.invert1 == a.invert2:
						newt.Exons.append( ITVL(a.chr2, a.startpos2, a.endpos2, newt.Exons[0].strand) )
					else:
						newt.Exons.append( ITVL(a.chr2, a.startpos2, a.endpos2, not newt.Exons[0].strand) )
				else:
					if a.invert1 and a.invert2:
						newt.Exons.append( ITVL(a.chr2, a.endpos2 - position_stopdodon - 3, a.endpos2, newt.Exons[0].strand) )
					elif (not a.invert1) and a.invert2:
						newt.Exons.append( ITVL(a.chr2, a.endpos2 - position_stopdodon - 3, a.endpos2, not newt.Exons[0].strand) )
					elif a.invert1 and (not a.invert2):
						newt.Exons.append( ITVL(a.chr2, a.startpos2, a.startpos2 + position_stopdodon + 3, not newt.Exons[0].strand) )
					elif (not a.invert1) and (not a.invert2):
						newt.Exons.append( ITVL(a.chr2, a.startpos2, a.startpos2 + position_stopdodon + 3, newt.Exons[0].strand) )
				assert(newt.Exons[-2].chr == a.chr1 and max(newt.Exons[-2].startpos, a.startpos1) < min(newt.Exons[-2].endpos, a.endpos1) )
				assert(newt.Exons[-1].chr == a.chr2 and max(newt.Exons[-1].startpos, a.startpos2) < min(newt.Exons[-1].endpos, a.endpos2))
				# update alteration_indicator
				alteration_indicator[tname] = True
				# update transcript
				transcripts[newt.TransID] = newt
				# update adjacency
				noveladjacencies.append( ADJACENCY(a.id, newt.Exons[-2].chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, not newt.Exons[-2].strand, newt.Exons[-1].chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, newt.Exons[-1].strand) )
				corresponding_trans.append(newt.TransID)
			elif t.Chr == a.chr2 and a.startpos2 < tss and a.endpos2 > tss and (a.invert2 == t.Strand):
				for e in t.Exons:
					if e.chr == a.chr2 and max(e.startpos, a.startpos2) < min(e.endpos, a.endpos2):
						newt.Exons.append( ITVL(e.chr, max(e.startpos, a.startpos2), min(e.endpos, a.endpos2), e.strand) )
				assert(len(newt.Exons) > 0)
				# the next 4 lines are commented out, that is, not extend the last unchanged exon to the inversion breakpoint
				# if t.Strand:
				# 	newt.Exons[-1] = ITVL(newt.Exons[-1].chr, newt.Exons[-1].startpos, a.endpos2, newt.Exons[-1].strand)
				# else:
				# 	newt.Exons[-1] = ITVL(newt.Exons[-1].chr, a.startpos2, newt.Exons[-1].endpos, newt.Exons[-1].strand)
				position_stopdodon = FindStopCodon(genome[a.chr1][(a.startpos1):(a.endpos1)]) if a.invert1 else FindStopCodon( ReverseComplement(genome[a.chr1][(a.startpos1):(a.endpos1)]) )
				if position_stopdodon == -1:
					if a.invert1 == a.invert2:
						newt.Exons.append( ITVL(a.chr1, a.startpos1, a.endpos1, newt.Exons[0].strand) )
					else:
						newt.Exons.append( ITVL(a.chr1, a.startpos1, a.endpos1, not newt.Exons[0].strand) )
				else:
					if a.invert1 and a.invert2:
						newt.Exons.append( ITVL(a.chr1, a.startpos1, a.startpos1 + position_stopdodon + 3, newt.Exons[0].strand) )
					elif (not a.invert1) and a.invert2:
						newt.Exons.append( ITVL(a.chr1, a.endpos1 - position_stopdodon - 3, a.endpos1, not newt.Exons[0].strand) )
					elif a.invert1 and (not a.invert2):
						newt.Exons.append( ITVL(a.chr1, a.startpos1, a.startpos1 + position_stopdodon + 3, not newt.Exons[0].strand) )
					elif (not a.invert1) and (not a.invert2):
						newt.Exons.append( ITVL(a.chr1, a.endpos1 - position_stopdodon - 3, a.endpos1, newt.Exons[0].strand) )
				assert(newt.Exons[-2].chr == a.chr2 and max(newt.Exons[-2].startpos, a.startpos2) < min(newt.Exons[-2].endpos, a.endpos2) )
				assert(newt.Exons[-1].chr == a.chr1 and max(newt.Exons[-1].startpos, a.startpos1) < min(newt.Exons[-1].endpos, a.endpos1) )
				# update alteration_indicator
				alteration_indicator[tname] = True
				# update transcript
				transcripts[newt.TransID] = newt
				# update adjacency
				noveladjacencies.append( ADJACENCY(a.id, newt.Exons[-2].chr, newt.Exons[-2].startpos, newt.Exons[-2].endpos, not newt.Exons[-2].strand, newt.Exons[-1].chr, newt.Exons[-1].startpos, newt.Exons[-1].endpos, newt.Exons[-1].strand) )
				corresponding_trans.append(newt.TransID)
		assert( len(corresponding_trans) == len(noveladjacencies) )
	# pop up the altered transcripts
	for tname in [tname for tname,v in alteration_indicator.items() if v]:
		transcripts.pop(tname, None)
	return transcripts, noveladjacencies, corresponding_trans


Nucleotide={'A':'T', 'C':'G', 'G':'C', 'T':'A', 'R':'Y', 'Y':'R', 'S':'W', 'W':'S', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N', '.':'.', '-':'-'}


def ReverseComplement(seq):
	rcseq = "".join([Nucleotide[x] for x in seq])
	return rcseq[::-1]


def ReadGenome(fafile):
	genome={}
	fp=open(fafile,'r')
	line=fp.readline().strip()
	tmpseq=''
	tmpname=''
	while line!='':
		if line[0]=='>':
			if len(tmpname)!=0:
				genome[tmpname]=tmpseq
			tmpseq=''
			strs=line.split()
			tmpname=strs[0][1:]
		else:
			tmpseq+=line
		line=fp.readline().strip()
	genome[tmpname]=tmpseq
	fp.close()
	return genome


def GetChrLengths(genome):
	return {k:len(v) for k,v in genome.items()}


def WriteGTF(outputfile, transcripts):
	fp = open(outputfile, 'w')
	for tname,t in transcripts.items():
		strandstr = "+"
		if not t.Strand:
			strandstr = "-"
		fp.write("{}\tVARIATION\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\";\n".format(t.Chr, t.StartPos, t.EndPos, strandstr, t.GeneID, t.TransID))
		for e in t.Exons:
			strandstr = "+"
			if not e.strand:
				strandstr = "-"
			fp.write("{}\tVARIATION\texon\t{}\t{}\t.\t{}\tgene_id \"{}\"; transcript_id \"{}\";\n".format(e.chr, e.startpos, e.endpos, strandstr, t.GeneID, t.TransID))
	fp.close()


def WriteTranscriptSequences(outputfile, transcripts, genome):
	fp = open(outputfile, 'w')
	for tname, t in transcripts.items():
		seq = ""
		for e in t.Exons:
			if e.strand:
				seq += genome[e.chr][e.startpos:e.endpos]
			else:
				seq += ReverseComplement(genome[e.chr][e.startpos:e.endpos])
		fp.write(">" + tname +"\n")
		count = 0
		while count < len(seq):
			fp.write(seq[count:min(count+70,len(seq))] + "\n")
			count += 70
	fp.close()


def WriteNovelAdjacency(outputfile, noveladjacencies_inv, noveladjacencies_tra):
	fp = open(outputfile, 'w')
	fp.write("# chr1\tstartpos1\tendpos1\tchr2\tstartpos2\tendpos2\tvar_id\tscore\tstrand1\tstrand2\n")
	for a in noveladjacencies_inv:
		strand1 = "+"
		if a.isleft1:
			strand1 = "-"
		strand2 = "+"
		if a.isleft2:
			strand2 = "-"
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\n".format(a.chr1, a.startpos1, a.endpos1, a.chr2, a.startpos2, a.endpos2, a.var_id, strand1, strand2))
	for a in noveladjacencies_tra:
		strand1 = "+"
		if a.isleft1:
			strand1 = "-"
		strand2 = "+"
		if a.isleft2:
			strand2 = "-"
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\n".format(a.chr1, a.startpos1, a.endpos1, a.chr2, a.startpos2, a.endpos2, a.var_id, strand1, strand2))
	fp.close()


def WriteNovelAdjacency_maptrans(outputfile, map_noveladj_trans_inv, map_noveladj_trans_tra):
	fp = open(outputfile, 'w')
	fp.write("# chr1\tstartpos1\tendpos1\tchr2\tstartpos2\tendpos2\tvar_id\tscore\tstrand1\tstrand2\ttranscript_ids\n")
	for a,tnames in map_noveladj_trans_inv.items():
		strand1 = "+"
		if a.isleft1:
			strand1 = "-"
		strand2 = "+"
		if a.isleft2:
			strand2 = "-"
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\n".format(a.chr1, a.startpos1, a.endpos1, a.chr2, a.startpos2, a.endpos2, a.var_id, strand1, strand2, ",".join(tnames)))
	for a,tnames in map_noveladj_trans_tra.items():
		strand1 = "+"
		if a.isleft1:
			strand1 = "-"
		strand2 = "+"
		if a.isleft2:
			strand2 = "-"
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\n".format(a.chr1, a.startpos1, a.endpos1, a.chr2, a.startpos2, a.endpos2, a.var_id, strand1, strand2, ",".join(tnames)))
	fp.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python RearrangedTranscripts.py <rsvsim_path> <genome_file> <gtf_file> <output_prefix>")
	else:
		rsvsim_path = sys.argv[1]
		genome_file = sys.argv[2]
		gtf_file = sys.argv[3]
		output_prefix = sys.argv[4]

		genome = ReadGenome(genome_file)
		transcripts = ReadGTF(gtf_file)
		#ChrLengths = GetChrLengths(genome)
		new_chrom_sizes = pd.read_csv("./chrom_sizes.txt", sep="\t", names=['chr','size'], header=None)
		ChrLengths = dict(new_chrom_sizes.to_dict(orient ="split")['data'])
		deletions = ReadDeletions(rsvsim_path + "/deletions.csv")
		inversions = ReadInversions(rsvsim_path + "/inversions.csv")
		translocations = ReadTranslocations(rsvsim_path + "/translocations.csv", ChrLengths)
		tmpdeletions, tmptranslocations = ReadTranslocations_insertion(rsvsim_path + "/insertions.csv", ChrLengths)
		deletions += tmpdeletions
		translocations += tmptranslocations

		tlocator = TranscriptLocator(transcripts)
		# transcripts = UpdateAnnotation_del(transcripts, deletions, tlocator, genome)
		transcripts, noveladjacencies_inv, corresponding_trans_inv = UpdateAnnotation_inv(transcripts, inversions, tlocator, genome)
		transcripts, noveladjacencies_tra, corresponding_trans_tra = UpdateAnnotation_tra(transcripts, translocations, tlocator, genome)
		
		# make unique of the novel adjacencies
		map_noveladj_trans_inv = {}
		map_noveladj_trans_tra = {}
		assert( len(corresponding_trans_inv) == len(noveladjacencies_inv) )
		assert( len(corresponding_trans_tra) == len(noveladjacencies_tra) )
		for i in range(len(noveladjacencies_inv)):
			if noveladjacencies_inv[i] in map_noveladj_trans_inv:
				map_noveladj_trans_inv[ noveladjacencies_inv[i] ].append( corresponding_trans_inv[i] )
			else:
				map_noveladj_trans_inv[ noveladjacencies_inv[i] ] = [ corresponding_trans_inv[i] ]
		for i in range(len(noveladjacencies_tra)):
			if noveladjacencies_tra[i] in map_noveladj_trans_tra:
				map_noveladj_trans_tra[ noveladjacencies_tra[i] ].append( corresponding_trans_tra[i] )
			else:
				map_noveladj_trans_tra[ noveladjacencies_tra[i] ] = [ corresponding_trans_tra[i] ]
		map_noveladj_trans_inv = {k:list(set(map_noveladj_trans_inv[k])) for k in map_noveladj_trans_inv.keys()}
		map_noveladj_trans_tra = {k:list(set(map_noveladj_trans_tra[k])) for k in map_noveladj_trans_tra.keys()}

		WriteGTF(output_prefix + "_annotation.gtf", transcripts)
		WriteTranscriptSequences(output_prefix+"_transcripts.fa", transcripts, genome)
		# WriteNovelAdjacency(output_prefix+"_adjacency.tsv", noveladjacencies_inv, noveladjacencies_tra)
		WriteNovelAdjacency_maptrans(output_prefix+"_adjacency.tsv", map_noveladj_trans_inv, map_noveladj_trans_tra)
