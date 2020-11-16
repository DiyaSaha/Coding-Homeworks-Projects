#!/usr/bin/env python3
# -*- coding: utf-8 -*-

greg_file = open('GregMendel_SNPs.txt','r')
lilly_file = open('LillyMendel_SNPs.txt','r')

g_id_list = []
g_chr = []
g_pos = []
g_SNP1 = []
g_SNP2 = []
for line in greg_file.readlines():
    line = line.strip()
    line.rstrip()
    g_id_list.append(line.split('chr')[0])
    g_SNP1.append(line.split('chr')[1].split('(')[1].split(',')[0])
    g_SNP2.append(line.split('chr')[1].split('(')[1].split(',')[1].split(')')[0])
    g_chr.append((line.split('chr')[1]).split('-')[0])
    g_pos.append(((line.split('chr')[1]).split('-')[1]).split('(')[0])

l_id_list = []
l_chr = []
l_pos = []
l_SNP1 = []
l_SNP2 = []
for line in lilly_file.readlines():
    line = line.strip()
    line.rstrip()
    l_id_list.append(line.split('chr')[0])
    l_SNP1.append(line.split('chr')[1].split('(')[1].split(',')[0])
    l_SNP2.append(line.split('chr')[1].split('(')[1].split(',')[1].split(')')[0])
    l_chr.append((line.split('chr')[1]).split('-')[0])
    l_pos.append(((line.split('chr')[1]).split('-')[1]).split('(')[0])
    
greg_file.close()
lilly_file.close()

   
#Part b 
def read_SNP_file(filename): 
    id_list = []
    chrom = []
    pos = []
    SNP = []
    file = open(str(filename),'r')
    for line in file.readlines():
        line = line.strip()
        line.rstrip()
        id_list.append(line.split('chr')[0])
        SNP1 = line.split('chr')[1].split('(')[1].split(',')[0]
        SNP2 = line.split('chr')[1].split('(')[1].split(',')[1].split(')')[0]
        SNP.append(SNP1+SNP2)
        chrom.append((line.split('chr')[1]).split('-')[0])
        pos.append(((line.split('chr')[1]).split('-')[1]).split('(')[0])
    dict_SNP = {'ID':id_list,'SNP': SNP,'chr':chrom, 'position':pos}
    return dict_SNP

Greg_dict = read_SNP_file('GregMendel_SNPs.txt')
Lilly_dict = read_SNP_file('LillyMendel_SNPs.txt')
 
#------------------------------------Part c------------------------------------ 
global_start = -10
global_end = -10
global_max =  0
local_start = -5
local_end = -5
local_max = 0
for index in range(len(Greg_dict['ID'])):
    mychr = Greg_dict['chr'][index]
    if mychr == '10':
        mypos = Greg_dict['position'][index]
        Greg_SNP = Greg_dict['SNP'][index]
        Lilly_SNP = Lilly_dict['SNP'][index]
        if Lilly_SNP == Greg_SNP: 
            if local_max == 0:
                local_start = Greg_dict['position'][index]
                local_max += 1
            else:
                local_max += 1
        else:
            local_end = Greg_dict['position'][index - 1]
            if local_max > global_max: #Current shared region has more SNPs than the maximum region so far
                global_max = local_max
                global_start = local_start
                global_end = local_end
            #reset/reinitialize your local variables
            local_max = 0
            local_start = 0
            local_end = 0

print('max shared SNPs region:',global_max)
print(global_start,global_end)
print('-----------------------------')
#--------------------------------Part d--------------------------------------- 
SNP_id = []
descrip = []
genotype = []
SNPdef_file = open('SNP_definitions.txt','r')

for line in SNPdef_file.readlines()[1:]:
    line = line.strip()
    line.rstrip()
    line = line.split('\t')
    SNP_id.append(line[0])
    genotype.append(line[1])
    descrip.append(line[2])
SNPdef_dict = {}
for i in range(len(descrip)):
    mykey = SNP_id[i], genotype[i]
    SNPdef_dict[mykey] = descrip[i]
        
SNPdef_file.close()

#---------------------------------Part e---------------------------------------
     
def disease(dictionary):
    disease_list = []
    for i in range(len(dictionary['chr'])):
        if (dictionary['chr'][i] == '9') and (22070000 < int(dictionary['position'][i]) < 22106000):
            disease_SNP = dictionary['SNP'][i]
            disease_ID = dictionary['ID'][i]
            for keys in SNPdef_dict.keys():
                if (disease_ID,disease_SNP) == keys:
                    disease_list.append(SNPdef_dict[keys])
    return print(*disease_list,sep = '\n')

print('Gregs risk for heart attack')
disease(Greg_dict)
print('-----------------------------')
print('Lillys risk for heart attack')
disease(Lilly_dict)

#----------------------------------Part f----------------------------------------

'''
mySNP = rs17817449(G,G)
This is not the SNP that shows the strongest association with body weight but is one of 
the co-inherited SNPs in the FTO gene region. This is on
chromosome 16 at position 53779455 

'''
#check if it exists in Lilly's and Greg's dictionaries, and if it does check if it 
#exists in the SNP descriptions 
    
def check_status(dictionary,mySNP):
    if str(mySNP) in dictionary['ID']:
        key = (str(mySNP),dictionary['SNP'][dictionary['ID'].index(str(mySNP))])
        pos = dictionary['position'][dictionary['ID'].index(str(mySNP))]
        if key in SNPdef_dict.keys():
            return (SNPdef_dict[key],pos)
print('-----------------------------')
print('Gregs status for mySNP and postion:',check_status(Greg_dict,'rs17817449'))
print('Lillys status for mySNP and position:',check_status(Lilly_dict,'rs17817449'))
