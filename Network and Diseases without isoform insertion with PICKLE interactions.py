import csv
import numpy as np
import pandas as pd
tiny = np.finfo(np.float32).tiny
import time


# Set starting time:
start_time = time.time();


# Importing Human proteome and Pickle interaction files:
# Converting Pickle binary interaction format to uniprot 'interacts with' format.
# Adding (replacing) interactions column of the Pickle database.
# Writing Human proteome data with Pickle interactions.

pickle_ints = pd.read_csv('UniProtNormalizedTabular-default.txt', delimiter='\t', usecols=['InteractorA', 'InteractorB'])
Human_proteome = pd.read_csv('uniprot-proteome%3AUP000005640.tab', delimiter='\t')

grp1 = pickle_ints.groupby('InteractorA')['InteractorB'].apply(list).apply(lambda x : '; '.join(str(elem) for elem in x)).reset_index(name='Interacts with').rename(columns={'InteractorA':'Entry'})
grp2 = pickle_ints.groupby('InteractorB')['InteractorA'].apply(list).apply(lambda x : '; '.join(str(elem) for elem in x)).reset_index(name='Interacts with').rename(columns={'InteractorB':'Entry'})

ints_in_uniprot_format = grp1.append(grp2).groupby('Entry')['Interacts with'].apply(list).apply(lambda x : '; '.join(str(elem) for elem in x)).reset_index(name='Interacts with')
ints_in_uniprot_format['Interacts with'] = ints_in_uniprot_format['Interacts with'].apply(lambda x : x.split('; ')).apply(lambda x : list(set(list(x)))).apply(lambda x : '; '.join(str(elem) for elem in x))
Human_proteome = Human_proteome.merge(ints_in_uniprot_format, how='left', on='Entry')
Human_proteome.to_csv('Human proteome with Pickle interactions.csv', index=False)


# Importing .csv file and creating dictionaries and lists:
# indx is a dictionary that will be used to find a given protein entry and return its index.
# ints is the list of interaction lists (list of lists).
# seqs is the list for sequence strings.


indx = {} #entry indexes dictionary
ints = [] #list of interactions
n_ints =[] #list of number of interactions
seqs = [] #list of sequence strings
Entry= [] #list of Entries
prot_name = [] #list of protein names
orpha_uniprot = [] #list of orphanet numbers
involvement = [] #list of involvement in diseases
status = [] #list of Sprot/trembl status
gene_name = [] #list of Gene names (primary)

with open('Human proteome with Pickle interactions.csv', newline = '') as csvfile:
    n = 0 #number of row
    for row in csv.DictReader(csvfile):
        indx[row['Entry']] = n
        #create list of lists for interactions, renaming "Itself", triming coforms and making empty lists as needed:
        ints.append([entry for entry in row['Interacts with'].split('; ')] if row['Interacts with'] is not '' else [])
        n_ints.append(row['Interacts with'].count('; ') + 1 if row['Interacts with'] is not '' else 0)
        seqs.append(row['Sequence'])
        Entry.append(row['Entry'])
        prot_name.append(row['Protein names'])
        orpha_uniprot.append(row['Cross-reference (Orphanet)'])
        involvement.append(row['Involvement in disease'])
        status.append(row['Status'])
        gene_name.append(row['Gene names  (primary )'])
        n += 1


# Here we count all of alphabet letters' repititions in sequence strings and form the rep_ent matrix, which its rows are representative of protein entries and its columns are representative of A-Z letters respectively.


rep_ent = np.zeros((n, 22), dtype = np.float32) #matrix of number of repitition of amino acids in entry sequences initialized with zero
for i in range(n):
    rep_ent[(i,  0)] = seqs[i].count('A')

    rep_ent[(i,  1)] = seqs[i].count('C')
    rep_ent[(i,  2)] = seqs[i].count('D') + seqs[i].count('B') / 2
    rep_ent[(i,  3)] = seqs[i].count('E') + seqs[i].count('Z') / 2
    rep_ent[(i,  4)] = seqs[i].count('F')
    rep_ent[(i,  5)] = seqs[i].count('G')
    rep_ent[(i,  6)] = seqs[i].count('H')
    rep_ent[(i,  7)] = seqs[i].count('I') + seqs[i].count('J') / 2

    rep_ent[(i, 8)] = seqs[i].count('K')
    rep_ent[(i, 9)] = seqs[i].count('L') + seqs[i].count('J') / 2
    rep_ent[(i, 10)] = seqs[i].count('M')
    rep_ent[(i, 11)] = seqs[i].count('N') + seqs[i].count('B') / 2
    rep_ent[(i, 12)] = seqs[i].count('O')
    rep_ent[(i, 13)] = seqs[i].count('P')
    rep_ent[(i, 14)] = seqs[i].count('Q') + seqs[i].count('Z') / 2
    rep_ent[(i, 15)] = seqs[i].count('R')
    rep_ent[(i, 16)] = seqs[i].count('S')
    rep_ent[(i, 17)] = seqs[i].count('T')
    rep_ent[(i, 18)] = seqs[i].count('U')
    rep_ent[(i, 19)] = seqs[i].count('V')
    rep_ent[(i, 20)] = seqs[i].count('W')

    rep_ent[(i, 21)] = seqs[i].count('Y')


# Calculating length of entry sequences, calculating probability of each letter per entry and then calculationg pit. The entry pit result is stored in pit_ent which is a (n * 1) matrix; it gives pit number for each entry.


len_ent = rep_ent @ np.ones((22, 1), dtype = np.float32) #length of entry sequences
P = rep_ent / (len_ent) #matrix of probability of each amino acid in each entry
pit_ent = -np.sum(P * np.log2(P + tiny), axis = 1, keepdims = True) #matrix for pit of each entry (used tiny to avoid log0)


# Calculating sum of pit * length for network members (except the hub). It is stored in pits_net which is a (n * 1) matrix; it gives pit * length number for each network entry.
# And besides that we form the network of reactions and create a matrix for sum of alphabet repititions in a network (columns and rows are like rep_ent).
# We do these two in one place to save time of dictionary look-ups.


y_summation = np.zeros((n, 1), dtype = np.float32) #matrix of sum of pit*length of network members initialized with zero
rep_net = np.copy(rep_ent) #matrix of number of repitition of amino acids in interacting networks initialized with rep_ent (to include network hubs)
y_concat = np.zeros((n,22), dtype = np.float32)
for i in range(n):
    for entry in ints[i]:
        try:
            ind = indx[entry]
            y_summation[(i, 0)] += pit_ent[(ind, 0)] * len_ent[(ind, 0)]
            rep_net[(i), 0:22] += rep_ent[(ind), 0:22]
            y_concat[(i), 0:22] += rep_ent[(ind), 0:22]
        except KeyError:
            None



# Calculating length of each network, calculating probability of each letter per network entry and then calculationg pit. The network pit result is stored in pit_net which is a (n * 1) matrix; it gives pit number for each network entry.


len_net = rep_net @ np.ones((22, 1), dtype = np.float32) #length of networks
len_y_concat = y_concat @ np.ones((22,1), dtype = np.float32) # length of ys
P = rep_net / (len_net) #matrix of probability of each amino acid in each network
P_y_concat = y_concat / (len_y_concat + tiny) #matrix of probability of each amino acid in ys
pit_net = -np.sum(P * np.log2(P + tiny), axis = 1, keepdims = True) #matrix for pit of each network (used tiny to avoid log0)
pit_y_concat = -np.sum(P_y_concat * np.log2(P_y_concat + tiny), axis = 1, keepdims = True) #matrix for pit of ys (used tiny to avoid log0)


# Calcualting network efficiency:


Total_net = len_net * pit_net
Total_y_concat = pit_y_concat * (len_net - len_ent)
B = pit_net * len_net - y_summation

Age = pd.read_csv('main_HUMAN.csv', usecols=['Entry', 'modeAge'])

orphadata_unprocessed = pd.read_csv('en_product9_prev.csv', encoding = "ISO-8859-1",
                                    usecols=['OrphaNumber', 'ExpertLink', 'Name', 'Name6', 'Source', 'Name12', 'Name20', 'Name24',
                                             'Name28']).rename(columns={'OrphaNumber': 'OrphaNumber (source: Orphanet)',
                                                                        'Name': 'Disease name (source: Orphanet)', 'Name6': 'Disease type',
                                                                        'Name12': 'Occurance type', 'Name20': 'Occurance value',
                                                                        'Name24': 'Geo distrib', 'Name28': 'Validity'})
orphadata_unprocessed['OrphaNumber (source: Orphanet)'] = orphadata_unprocessed['OrphaNumber (source: Orphanet)'].astype(str)

orphadata = pd.read_csv('en_product9_prev.csv', encoding = "ISO-8859-1",
                        usecols=['OrphaNumber', 'ExpertLink', 'Name', 'Name6', 'Source', 'Name12', 'Name20', 'Name24',
                                 'Name28']).rename(columns={'OrphaNumber': 'OrphaNumber (source: Orphanet)',
                                                            'Name': 'Disease name (source: Orphanet)', 'Name6': 'Disease type',
                                                            'Name12': 'Occurance type', 'Name20': 'Occurance value',
                                                            'Name24': 'Geo distrib', 'Name28': 'Validity'})
orphadata = orphadata[orphadata['Geo distrib'] == 'Worldwide'] # filtering out all locally distributed diseases -- Only keeps 'Worldwide'
orphadata = orphadata[orphadata['Occurance value'] != 'Unknown'] # filtering out 'Unknown' occurences
orphadata = orphadata[orphadata['Occurance value'] != 'Not yet documented'] # filtering out missing occurence values
orphadata = orphadata[orphadata['Occurance type'] != 'Lifetime Prevalence'] # filtering out 'Lifetime Prevalence'
orphadata = orphadata[orphadata['Occurance type'] != 'Cases/families'] # filtering out 'Cases/families'
orphadata = orphadata.dropna(subset=['Occurance value'])
orphadata.drop_duplicates('OrphaNumber (source: Orphanet)', keep='first', inplace=True)
orphadata['OrphaNumber (source: Orphanet)'] = orphadata['OrphaNumber (source: Orphanet)'].astype(str)

network = pd.DataFrame({'Status': status,
                        'Entry': Entry,
                        'Gene name': gene_name,
                        'Protein name': prot_name,
                        'Protein length': len_ent[:, 0],
                        'Protein entropy (in pits)': pit_ent[:, 0],
                        'Number of interactions': n_ints,
                        'Network length': len_net[:, 0],
                        'Network entropy (in pits)': pit_net[:, 0],
                        'Total network entropy (in pits)': Total_net[:, 0],
                        'Scalar summation of interactors entropy (in pits)': y_summation[:, 0],
                        'Estimation of the joint entropy of interactors (in pits)': pit_y_concat[:, 0],
                        'Total joint entropy estimation (in pits)': Total_y_concat[:, 0],
                        'Mutual information estimation (in pits)': B[:, 0],
                        'Orpha no. (source: UniProt)': orpha_uniprot,
                        'Involvement in disease (source: UniProt)': involvement})
network = network.merge(Age, how='left', on='Entry').rename(columns={'modeAge': 'Gene age'})
out1 = network.copy()
out1.to_csv('Disease probability.csv', index = False)

network['Orpha no. (source: UniProt)'] = network['Orpha no. (source: UniProt)'].apply(lambda x : x.split(';')[:-1]).apply(list)
network = network.explode('Orpha no. (source: UniProt)') # requires pandas version 0.25 or later
out2 = network.merge(orphadata_unprocessed, how='left', left_on = 'Orpha no. (source: UniProt)', right_on = 'OrphaNumber (source: Orphanet)')
out2 = out2.reindex(columns=['Orpha no. (source: UniProt)',
                             'Involvement in disease (source: UniProt)', 
                             'OrphaNumber (source: Orphanet)',
                             'Disease name (source: Orphanet)',
                             'Disease type', 
                             'Occurance value',
                             'Occurance type',
                             'Geo distrib',
                             'Validity',
                             'Source',
                             'Status',
                             'Entry',
                             'Gene name',
                             'Protein name',
                             'Protein length',
                             'Protein entropy (in pits)',
                             'Number of interactions',
                             'Network length',
                             'Network entropy (in pits)',
                             'Total network entropy (in pits)',
                             'Scalar summation of interactors entropy (in pits)',
                             'Estimation of the joint entropy of interactors (in pits)',
                             'Total joint entropy estimation (in pits)', 
                             'Mutual information estimation (in pits)',
                             'Gene age',
                             'ExpertLink'])
out2['Orpha no. (source: UniProt)'] = out2['Orpha no. (source: UniProt)'].astype(float)
out2 = out2.sort_values(by=['Orpha no. (source: UniProt)'])
out2.to_csv('Disease analysis of all proteins.csv', index = False)

out3 = network.merge(orphadata, how='left', left_on = 'Orpha no. (source: UniProt)', right_on = 'OrphaNumber (source: Orphanet)')
out3 = out3.dropna(subset=['OrphaNumber (source: Orphanet)'])
out3 = out3.reindex(columns=['OrphaNumber (source: Orphanet)',
                             'Disease name (source: Orphanet)',
                             'Involvement in disease (source: UniProt)',
                             'Disease type',
                             'Occurance value',
                             'Occurance type',
                             'Geo distrib',
                             'Validity',
                             'Source',
                             'Status',
                             'Entry',
                             'Gene name',
                             'Protein name',
                             'Protein length',
                             'Protein entropy (in pits)',
                             'Number of interactions',
                             'Network length',
                             'Network entropy (in pits)',
                             'Total network entropy (in pits)',
                             'Scalar summation of interactors entropy (in pits)',
                             'Estimation of the joint entropy of interactors (in pits)',
                             'Total joint entropy estimation (in pits)',
                             'Mutual information estimation (in pits)',
                             'Gene age',
                             'ExpertLink'])
out3.to_csv('Disease occurances associations.csv', index=False)

out4 = network.merge(orphadata, how='left', left_on = 'Orpha no. (source: UniProt)', right_on = 'OrphaNumber (source: Orphanet)')
out4['Orpha no. (source: UniProt)'] = out4['Orpha no. (source: UniProt)'].fillna(0)
out4['OrphaNumber (source: Orphanet)'] = out4['OrphaNumber (source: Orphanet)'].fillna(0)
out4['Occurance value'] = out4.apply(lambda row : row['Occurance value'] if row['Orpha no. (source: UniProt)'] == row['OrphaNumber (source: Orphanet)'] else tiny, axis=1)
out4['Occurance value'] = out4['Occurance value'].fillna(0)
out4['Occurance value'] = out4['Occurance value'].replace({'<1 / 1 000 000' : 0.000001, '1-9 / 100 000': 0.00005, '1-9 / 1 000 000': 0.000005,
                                                           '1-5 / 10 000': 0.0003, '>1 / 1000': 0.001, '6-9 / 10 000': 0.00075})

Tot_occur = out4[out4['Occurance value'] != 0].groupby('Entry')['Occurance value'].apply(sum).reset_index(name='Total occurance value')
Tot_occur['Total occurance value'] = Tot_occur['Total occurance value'].apply(lambda x : 'NA' if x < 0.000000001 else x)
out4 = out1.merge(Tot_occur, how='left', on='Entry')
out4['Total occurance value'] = out4['Total occurance value'].fillna(0)
out4 = out4.reindex(columns=['Status',
                             'Entry',
                             'Gene name',
                             'Protein name',
                             'Protein length',
                             'Protein entropy (in pits)',
                             'Number of interactions',
                             'Network length',
                             'Network entropy (in pits)',
                             'Total network entropy (in pits)',
                             'Scalar summation of interactors entropy (in pits)',
                             'Estimation of the joint entropy of interactors (in pits)',
                             'Total joint entropy estimation (in pits)',
                             'Mutual information estimation (in pits)',
                             'Orpha no. (source: UniProt)',
                             'Total occurance value',
                             'Involvement in disease (source: UniProt)',
                             'Gene age'])
out4.to_csv('Total occurances per disease.csv', index=False)


# Calculating total time elapsed. It should be around half a minute!


end_time = time.time()
print('Elapsed time:', end_time - start_time, 'seconds')

