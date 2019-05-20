#!/bin/python

import numpy as np
import collections

ITVL = collections.namedtuple("ITVL", ["chr", "startpos", "endpos", "strand"])

class Transcript_t(object):
	def __init__(self, _TransID, _GeneID, _Chr, _Strand, _StartPos, _EndPos):
		self.TransID=_TransID
		self.GeneID=_GeneID
		self.Chr = _Chr
		self.Strand = _Strand
		self.StartPos = _StartPos
		self.EndPos = _EndPos
		self.Exons = [] # each exon is an ITVL object
	def __eq__(self, other):
		if isinstance(other, Transcript_t):
			return (self.Chr==other.Chr and self.Strand==other.Strand and len(self.Exons)==len(other.Exons) and \
				np.all([self.Exons[i]==other.Exons[i] for i in range(len(self.Exons))]))
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, Transcript_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos<other.StartPos
			else:
				return self.EndPos<other.EndPos
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, Transcript_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos>other.StartPos
			else:
				return self.EndPos>other.EndPos
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result


def GetFeature(line, key):
	s=line.index(key)
	t=line.index(";", s+1)
	return line[(s+len(key)+2):(t-1)]


def ReadGTF(gtffile):
	Transcripts={}
	strand=""
	fp=open(gtffile, 'r')
	tmptransname=""
	tmptranscript=None
	extraExons = []
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[2]=="transcript":
			if tmptransname!="" and not (tmptranscript is None):
				Transcripts[tmptransname]=tmptranscript
			tmptransname=GetFeature(line, "transcript_id")
			tmpgeneid=GetFeature(line, "gene_id")
			tmptranscript=Transcript_t(tmptransname, tmpgeneid, strs[0], (strs[6]=="+"), int(strs[3])-1, int(strs[4]))
		elif strs[2]=="exon":
			thistransid=GetFeature(line, "transcript_id")
			if thistransid == tmptransname and not (tmptranscript is None):
				tmptranscript.Exons.append( ITVL(strs[0], int(strs[3])-1, int(strs[4]), (strs[6]=="+")) )
			else:
				extraExons.append([thistransid, GetFeature(line, "gene_id"), ITVL(strs[0], int(strs[3])-1, int(strs[4]), (strs[6]=="+"))])
	if tmptransname!="" and not (tmptranscript is None):
		Transcripts[tmptransname]=tmptranscript
	for e in extraExons:
		if e[0] in Transcripts:
			Transcripts[e[0]].Exons.append( e[2] )
		else:
			tmptranscript = Transcript_t(e[0], e[1], e[2].chr, e[2].strand, e[2].startpos, e[2].endpos)
			tmptranscript.Exons.append( e[2] )
			Transcripts[tmptranscript.TransID] = tmptranscript
	for t in Transcripts.keys():
		Transcripts[t].StartPos = min([e.startpos for e in Transcripts[t].Exons])
		Transcripts[t].EndPos = max([e.endpos for e in Transcripts[t].Exons])
		Transcripts[t].Exons.sort(key=lambda x:x.startpos)
		if not Transcripts[t].Strand:
			Transcripts[t].Exons = Transcripts[t].Exons[::-1]
	fp.close()
	return Transcripts


def Map_Gene_Trans(Transcripts):
	GeneTransMap={}
	TransGeneMap={}
	for v in Transcripts.values():
		TransGeneMap[v.TransID]=v.GeneID
		if v.GeneID in GeneTransMap:
			GeneTransMap[v.GeneID].append(v.TransID)
		else:
			GeneTransMap[v.GeneID]=[v.TransID]
	for g,v in GeneTransMap.items():
		sortedv = sorted(v)
		GeneTransMap[g] = sortedv
	return [GeneTransMap, TransGeneMap]


def GetTransLength(Transcripts):
	TransLength = {t:np.sum([e.endpos-e.startpos for e in v.Exons]) for t,v in Transcripts.items()}
	return TransLength


def ReadTranscriptFasta(filename, namesplitter = " "):
	TransSequence = {}
	fp = open(filename, 'r')
	name = ""
	seq = ""
	for line in fp:
		if line[0] == '>':
			if len(name) != 0:
				TransSequence[name] = seq
			name = line.strip().split(namesplitter)[0][1:]
			seq = ""
		else:
			seq += line.strip()
	if len(name) != "":
		TransSequence[name] = seq
	fp.close()
	return TransSequence


class TranscriptLocator(object):
	def __init__(self, transcripts):
		self.trans_ranges = [(tname, t.Chr, t.StartPos, t.EndPos) for tname,t in transcripts.items()] # tuple of <trans ID, trans chr, trans start, trans end>
		self.trans_ranges.sort(key = lambda x:(x[1], x[2], x[3]))

	def BinarySearchPosition(self, chr, pos, surroundings = 100):
		lo = 0
		hi = len(self.trans_ranges)
		while lo < hi:
			mid = int( (lo + hi) / 2)
			if self.trans_ranges[mid][1] < chr or (self.trans_ranges[mid][1] == chr and self.trans_ranges[mid][3] <= pos):
				lo = mid + 1
			elif self.trans_ranges[mid][1] == chr and self.trans_ranges[mid][2] <= pos and self.trans_ranges[mid][3] > pos:
				lo = mid
				hi = mid
				break
			elif self.trans_ranges[mid][1] > chr or (self.trans_ranges[mid][1] == chr and self.trans_ranges[mid][2] > pos):
				hi = mid - 1
			else:
				print("error: chr = {}, pos = {}, lo = {}, mid = {}, hi = {}".format(chr, pos, self.trans_ranges[lo], self.trans_ranges[mid], self.trans_ranges[hi]))
		# check the surroundings of lo and hi to include all genes that cover pos
		tids = []
		i = lo
		while i >= 0 and i >= lo - surroundings:
			if self.trans_ranges[i][1] == chr and self.trans_ranges[i][2] <= pos and self.trans_ranges[i][3] > pos:
				tids.append(self.trans_ranges[i][0])
			i -= 1
		i = hi
		while i < len(self.trans_ranges) and i <= hi + 50:
			if self.trans_ranges[i][1] == chr and self.trans_ranges[i][2] <= pos and self.trans_ranges[i][3] > pos:
				tids.append(self.trans_ranges[i][0])
			i += 1
		return list(set(tids))
