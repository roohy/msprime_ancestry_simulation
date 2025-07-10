#This code was originally written and provided to me by Dr. Christopher Gignoux
# first, calculate LD scores for each SNP, and dictionaries of SNPs "in LD"

import numpy as np 
import sys

class LDScore:
    def __init__(self) -> None:
        self.ld_scores = {}
        self.in_ld = {}

    def add_ld_score(self, snp_id, r2):
        if snp_id not in self.ld_scores:
            self.ld_scores[snp_id] = 1 + r2
        else:
            self.ld_score[snp_id]  += r2 
        if snp_id not in self.in_ld:
            self.in_ld[snp_id] = set()
    def add_ld_pair(self,snp1,snp2):
        self.in_ld[snp1].add(snp2)
        self.in_ld[snp2].add(snp1)
    def 

def main():
    input_addr = sys.argv[1]
    min_r2 = float(sys.argv[2])
    output_addr = sys.argv[3]
    ld_score_keeper = LDScore()

    in_ld = {}
    ld_scores = {}
    with open(input_addr,'r') as input_file:
        header = input_file.readline()
        for line in input_file:
            data = line.strip().split()
            r2 = float(line[-1])
            snp1 = line[2]
            snp2 = line[5]
            ld_score_keeper.add_ld_score(snp1,r2)
            ld_score_keeper.add_ld_score(snp2,r2)
            if r2 >= min_r2:
                ld_score_keeper.add_ld_pair(snp1,snp2)
        

#min_r2 = 0.5

#infile = 'plink.ld'

# first go through the file, calculate ld scores and what snps are in pairs

in_ld = {}
ld_scores = {}

for i, line in enumerate(open('test_files/'+infile)):
	if i > 0:
		line = line.strip().split()
		r2 = float(line[-1])
		snp1 = line[2]
		snp2 = line[5]
		for snp in (snp1,snp2):
			if snp in ld_scores:
				ld_scores[snp] += r2
			else:
				ld_scores[snp] = 1+r2
			if snp not in in_ld:
				in_ld[snp] = set()
		if r2>=min_r2:
			in_ld[snp1].add(snp2)
			in_ld[snp2].add(snp1)

print('found %s total snps' % (len(ld_scores)))

# greedy tag snp selection

rsids = [i for i in ld_scores.keys()]
ld_scoresc = [i for i in ld_scores.values()]

in_tags = set()

outfile = open(outfilename,'w')
print('writing output to %s' % (outfile))

print('running greedy portion now')
ranked_rsids = array(rsids)[argsort(ld_scoresc)[::-1]]
ranked_ldscores = array(ld_scoresc)[argsort(ld_scoresc)[::-1]]
# set up the greedy algorithm
in_tag_ld = set()
for i,rsid in enumerate(ranked_rsids):
	if rsid in in_tag_ld:
		continue
	else:
		in_tags.add(rsid)
		outfile.write('%s\t%s\n' % (rsid,ranked_ldscores[i]))
		[in_tag_ld.add(j) for j in in_ld[rsid]]

print('found %s tags' % (len(in_tags)))
outfile.close()
