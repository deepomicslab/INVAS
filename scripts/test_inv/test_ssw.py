import sys
import os.path as op
import argparse as ap
import ctypes as ct
import timeit as ti
import gzip
import math
import ssw_lib

sLibPath="./"
nMatch=2
nMismatch=2
nOpen=3
nExt=1
sMatrix=""
nThr=0
bBest=False

def to_int(seq, lEle, dEle2Int):
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num

def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, nMaskLen)

    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)

def buildPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = ''
    sA = ''
    sR = ''
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += q[nQOff : nQOff+n]
            sA += ''.join(['|' if q[nQOff+j] == r[nROff+j] else '*' for j in range(n)])
            sR += r[nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += q[nQOff : nQOff+n]
            sA += ' ' * n
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sA += ' ' * n
            sR += r[nROff : nROff+n]
            nROff += n
    return sCigar, sQ, sA, sR



# init DNA score matrix
dEle2Int = {}
dInt2Ele = {}
lEle = ['A', 'C', 'G', 'T', 'N']
dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'T', 'c':'G', 'g':'C', 't':'A'} 
for i,ele in enumerate(lEle):
    dEle2Int[ele] = i
    dEle2Int[ele.lower()] = i
    dInt2Ele[i] = ele
nEleNum = len(lEle)
lScore = [0 for i in range(nEleNum**2)]
for i in range(nEleNum-1):
    for j in range(nEleNum-1):
        if lEle[i] == lEle[j]:
            lScore[i*nEleNum+j] = nMatch
        else:
            lScore[i*nEleNum+j] = -nMismatch

# translate score matrix to ctypes
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore

ssw = ssw_lib.CSsw(sLibPath)

# query info
sQId,sQSeq,sQQual = "","CGGGTTCATGTGGGAGGACCTACTTACGTGGTCACGGTACTCG",""
# build query profile
qNum = to_int(sQSeq, lEle, dEle2Int)
qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)

# build rc query profile

sQRcSeq = ''.join([dRc[x] for x in sQSeq[::-1]])
qRcNum = to_int(sQRcSeq, lEle, dEle2Int)
qRcProfile = ssw.ssw_init(qRcNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)

# set mask len
nMaskLen = len(sQSeq) // 2

# target info
sRId,sRSeq,_ = "","GCCGGGTTGGGGTGTTGGGCCCTGGAGGGTGCACAGACTCTCCTCTCGGCCCGGACCCCCAGGCCCAAGTACACCCTCCTGGATGAATGCACCAGTGCCATGAGCATCGACGTGGAAGGCAAGATCTTCCAGGCGGCCAAGGACGCAGGCATTGCCCTGCTCTCCATCACCCACCGGCCCTCCCTGTGGTAGGTGCCCTGTCTCCCTTCCTGGGGTGAGTGGGAGTGGCTGCCTGAGGGGAGGAGGTGGC",""
rNum = to_int(sRSeq, lEle, dEle2Int)

# format of res: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
nFlag = 0
res = align_one(ssw, qProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)

resPrint = res
strand = 0
sCigar, sQ, sA, sR = buildPath(sQSeq, sRSeq, res[4], res[2], res[8])
print(resPrint)
print('target_name: {}\nquery_name: {}\noptimal_alignment_score: {}\t'.format(sRId, sQId, resPrint[0]))
if resPrint[1] > 0:
    print('suboptimal_alignment_score: {}\t'.format(resPrint[1])),
if resPrint[2] + 1:
    print('target_begin: {}\t'.format(resPrint[2] + 1)),
print('target_end: {}\t'.format(resPrint[3] + 1)),
if resPrint[4] + 1:
    print('query_begin: {}\t'.format(resPrint[4] + 1))
print('query_end: {}\n'.format(resPrint[5] + 1))
print("cigar:",sCigar)
print("strand")
print('target_begin: {}'.format(resPrint[2] + 1))
print('target_end: {}'.format(resPrint[3] + 1))