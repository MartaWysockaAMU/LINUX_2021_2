import re
codons = {}
stop = [] # TAA, TAG, TGA
with open("codon_table.txt") as file1:
    for line in file1:
        nt = line.strip().split()[0]
        aa = line.strip().split()[1]
        if aa == "*":
            stop.append(nt)
            codons[nt] = aa
        else:
            codons[nt] = aa

file2 = open("rosalind_orf.txt").readlines()
seq = ""
for line in file2:
    if not line.startswith(">"):
        seq += line.strip()

starts = [orf.start() for orf in re.finditer('ATG', seq)]
ORFs = {}
for start in starts:
    ORFs[start] = ["", False]
    for i in range(start,len(seq),3):
        codon = ""
        for j in range(i,i+3):
            codon += seq[j]
        if codon in stop:
            ORFs[start][1] = True
            break
        else:
            ORFs[start][0] += codons[codon]
reverse = {"A":"T","T":"A","G":"C","C":"G"}
seq2 = ""
for s in seq[::-1]:
    seq2 += reverse[s]

starts = [orf.start() for orf in re.finditer('ATG', seq2)]
#print(starts)
for start in starts:
    ORFs["-"+str(start)] = ["", False]
    for i in range(start,len(seq2),3):
        codon = ""
        for j in range(i,i+3):
            codon += seq2[j]
        #print(codon, codons[codon])
        if codon in stop:
            ORFs["-"+str(start)][1] = True
            break
        else:
            ORFs["-"+str(start)][0] += codons[codon]

aminos = []
for start,aa in ORFs.items():
    if aa[1]:
        if aa[0] not in aminos:
            aminos.append(aa[0])

print("\n".join(aminos))

