#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 22:45:26 2020

@author: diyasaha
"""
#Data files too big to upload to GitHub
input_file = open('slim_chr2_seq.fasta','r')
name = []
seq = []

content = input_file.readlines() #excluding the reference
content = [x.strip() for x in content]
input_file.close()

i = 1
for line in content:
    if i % 2 == 0 :
        seq.append(line)
    else :
        name.append(line)       
    i += 1

    
file = open('slim_chr2_SNPS.vcf','r')
Chr = []
rsid = []
slim_pos = []
genome_pos = []
ref = []
alt = []
geneinfo = []
for x in file.readlines()[1:]:
    x = x.split("\t")
    Chr.append(x[0])
    rsid.append(x[1])
    slim_pos.append(x[2])
    genome_pos.append(x[3])
    ref.append(x[4])
    alt.append(x[5])
    geneinfo.append(x[6])
file.close()

#Part 1, Q3
array = []
frequencies = []
for l in range(len(alt)):
    counter = 0
    for k in range(1,len(seq)):
        if seq[k][int(slim_pos[l])-1] == alt[l]:
            counter += 1
    array.append(counter)

for x in array: 
    frequencies.append(x/90) #Frequency = how many times variant occurs/all the seq

#Part 1, Q4
min_freq_index = [i for i, j in enumerate(frequencies) if j == min(frequencies)]
max_freq_index = [i for i, j in enumerate(frequencies) if j == max(frequencies)]

min_frequencies = []
max_frequencies = []
for b in min_freq_index: 
    min_frequencies.append(rsid[b])
for b in max_freq_index:
    max_frequencies.append(rsid[b])

print('SNPs with lowest occurence:', min_frequencies)
print('SNPs with highest occurence:', max_frequencies)

#Part 1, Q5
stdout = open('slim_chr2_SNPs_withfrequencies.vcf','w')

print('Chromosome_Numb','\t SNP_Name','\t SNP_slim_pos','\t SNP_ref','\t SNP_alt','\t SNP_freq','\t Gene_info', file = stdout)

for j in range(len(frequencies)):
    print(Chr[j],'\t',rsid[j],'\t',slim_pos[j],'\t',ref[j],'\t',alt[j],'\t',frequencies[j],'\t',geneinfo[j],file = stdout)
stdout.close()

#Part 1 Optional 
def mean(freq = list):
    return sum(freq)/len(freq)
print('mean across all SNP frequencies in the population=', mean(frequencies))

#Part 2 
output1 = open('lymphoma_variants_small.txt','w')
print('Individual_ID','\t','variantID','\t','GeneDiseaseInfo', file = output1)
#index 79 - 82 
for h in range(79,83):
    print(name[1:][h],rsid[h],geneinfo[h],sep = '\t', file = output1)
output1.close()

'''
Part 2 Q3

SNP: rs3948464
What kind of SNP: A to G change 
Function of gene it is in and its association with the disease: major role in the outcome of tuberculosis infection. 27 SNPs 
in the SP110 gene in 219 Gambian familities and identified 3 that were associated with tuberculosis, including C/T in exon 
11 (rs3948464). 
SNP frequency: 0.9889
'''

#Part 3

dict_SNP = {}
gene_ref = seq[0]

for position in range(len(gene_ref)):
    for individual_seq in seq: 
        if individual_seq[position] != gene_ref[position]:
            mykey = str(position+1) + ";" + individual_seq[position] + ';' + gene_ref[position]        
            if mykey in dict_SNP:
                dict_SNP[mykey] = dict_SNP[mykey] + 1
            else:
                dict_SNP[mykey] = 1

output2 = open('novel_variants.txt','w')

keys_lst = dict_SNP.keys()
print('SNP postion', '\t', 'Ref SNP', '\t', 'Alt SNP', '\t', 'Alt SNP Freq', file = output2)
    
for key,value in dict_SNP.items():
    pos_in_seq, alt_SNP, ref_SNP = key.split(';')
    SNP_freq = value/len(seq)
    print(pos_in_seq,'\t',alt_SNP,'\t',ref_SNP,'\t',SNP_freq, file = output2)


output2.close()


