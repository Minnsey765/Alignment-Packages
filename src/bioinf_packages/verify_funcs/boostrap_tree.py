# randomly sample a given gene into smaller constituant parts and align it

from Bio import SeqIO
import random
from Bio import Entrez

#seq = sequence
#n = subsequence length
#N = number of subsequences sampled
#output1 = output directory of fastas
#output2 = output directory of nexus

def boostrap_tree(seq: str, n: int, N: int, output1: str, output2: str):
    #resample sequence for N times
    for i in range(N):
        start = random.randint(0, len(seq) - n)
        sample = seq[start:start+n]
        