# estimates synonymous, non-synonymous, and intergenic  SNPs
# from a list of SNPs, a gene annotation file, and reference genome
# outputs bar plot and text file list of SNP details

# import required libraries
import matplotlib.pyplot as plt
import csv
import numpy as np

csv.field_size_limit(131072000000000000)

#dictionary of codons for each amino acid
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I',
    'ATG':'M', 'ACA':'T', 'ACC':'T',
    'ACG':'T', 'ACT':'T', 'AAC':'N',
    'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R',
    'AGG':'R', 'CTA':'L', 'CTC':'L',
    'CTG':'L', 'CTT':'L', 'CCA':'P',
    'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q',
    'CAG':'Q', 'CGA':'R', 'CGC':'R',
    'CGG':'R', 'CGT':'R', 'GTA':'V',
    'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A',
    'GCT':'A', 'GAC':'D', 'GAT':'D',
    'GAA':'E', 'GAG':'E', 'GGA':'G',
    'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S',
    'TCT':'S', 'TTC':'F', 'TTT':'F',
    'TTA':'L', 'TTG':'L', 'TAC':'Y',
    'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*',
    'TGG':'W'}

# translation function
def translate(seq):
    protein = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        protein += codon_table[codon]
    return(protein)

# reverse complement function
def reverse_complement(seq):
    comp = ""
    result = ""
    for i in seq:
        if i == 'A':
            comp = comp + 'T'
        elif i == 'C':
            comp = comp + 'G'
        elif i == 'T':
            comp = comp + 'A'
        elif i == 'G':
            comp = comp + 'C'
    for i in range(len(comp) - 1, -1, -1):
        result += comp[i]
    return(result)

# open genes file
genes_file = open('genes.txt')
genes_data = csv.reader(genes_file, delimiter="\t")

# open genome fasta
genome_file = open('genome.fasta')
genome_data = csv.reader(genome_file, delimiter="\t")

# open SNP file
snp_file = open('SNPs.txt')
snp_data = csv.reader(snp_file, delimiter="\t")

# create output table
output_file = open("output_table", "w+")
output_file.write('contig\tsnpPosition\tType\tAnnotation\tgeneID\t'
                  'referenceCodon\tmodifiedCodon\treferenceAA\tmodifiedAA\n')

# convert iterable to list
snp_list = list(snp_data)
genes_list = list(genes_data)
genome_list = list(genome_data)

# progress counter
i=0

# iterate through each row in the snp file, [1:-1] - skip header and ignore the newline at the end of the file
for row in snp_list[1:-1]:
    # show progress to user
    i += 1
    print('Currently processing SNP #' + str(i))
    print(row)



    # create variable for keeping track if SNP is intergenic
    intergenic_counter = 0

    # analyze SNP
    for line in genes_list:
        # search for intergenic SNPs
        if row[0] == line[0] \
                and (int(line[2]) <= int(row[1]) <= int(line[3])
                     or int(line[2]) >= int(row[1]) >= int(line[3])):

            # write contig number to output file
            output_file.write(row[0])

            # write SNP position within contig to output file
            output_file.write('\t' + row[1])

            #print(line[2] + ' ' + row[1] + ' ' + line[3])
            intergenic_counter += 1

            # determine if SNP is synonymous or non-synonymous substitution
            if row[0] == 'Contig_1':
                
                # reverse transcribe if start > end position
                if int(line[3]) > int(line[2]):
                    gene_sequence = genome_list[1][0][int(line[2]):int(line[3]) + 1]
                    gene_sequence = reverse_complement(gene_sequence)

                    gene_end = int(line[3])
                    snp_pos = int(row[1])
                    snp_pos_adjusted = snp_pos - gene_end
                    snp_nt = reverse_complement(str(row[2]))

                    # finds which codon the SNP exists in
                    codon_number = int((snp_pos_adjusted) / 3) + 1

                    # creates a duplicate of the respective gene sequence with the respective SNP replaced
                    gene_snp = gene_sequence[:snp_pos_adjusted] + snp_nt + gene_sequence[snp_pos_adjusted + 1:]

                    # gets the reference and modified codon sequences
                    codon_reference = gene_sequence[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]
                    codon_modified = gene_snp[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

                    print(gene_sequence)
                    print(snp_pos)
                    print('reversed')
                    print(codon_reference)
                    print(codon_modified)

                    if translate(codon_reference) == translate(codon_modified):
                        print('Synonymous mutation')
                        substitution_type = 'Synonymous mutation'

                    else:
                        print('Non-synonymous mutation')
                        substitution_type = 'Non-Synonymous mutation'

                else:
                    gene_sequence = genome_list[1][0][int(line[3]):int(line[2]) + 1]
                    gene_sequence = gene_sequence[::-1]
                    gene_start = int(line[2])
                    snp_pos = int(row[1])
                    snp_pos_adjusted = snp_pos - gene_start
                    snp_nt = str(row[2])

                    codon_number = int((snp_pos_adjusted) / 3) + 1

                    gene_snp = gene_sequence[:snp_pos_adjusted] + snp_nt + gene_sequence[snp_pos_adjusted + 1:]

                    codon_reference = gene_sequence[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]
                    codon_modified = gene_snp[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

                    print(gene_sequence)
                    print(snp_pos)
                    print(codon_reference)
                    print(codon_modified)

                    if translate(codon_reference) == translate(codon_modified):
                        print('Synonymous mutation')
                        substitution_type = 'Synonymous mutation'
                    else:
                        print('Non-synonymous mutation')
                        substitution_type = 'Non-Synonymous mutation'

            elif row[0] == 'Contig_2':
                print(row[1])
                print("hello")
                print(line[2])
                print(line[3])
                print("goodbye")
                
                if int(line[3]) > int(line[2]):
                    gene_sequence = genome_list[3][0][int(line[2]):int(line[3]) + 1]
                    gene_sequence = reverse_complement(gene_sequence)

                    gene_end = int(line[3])
                    snp_pos = int(row[1])
                    snp_pos_adjusted = snp_pos - gene_end
                    np_nt = reverse_complement(str(row[2]))

                    codon_number = int((snp_pos_adjusted) / 3) + 1

                    gene_snp = gene_sequence[:snp_pos_adjusted] + snp_nt + gene_sequence[snp_pos_adjusted + 1:]

                    codon_reference = gene_sequence[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]
                    codon_modified = gene_snp[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

                    print(gene_sequence)
                    print(snp_pos)
                    print('reversed')
                    print(codon_reference)
                    print(codon_modified)

                    if translate(codon_reference) == translate(codon_modified):
                        print('Synonymous mutation')
                        substitution_type = 'Synonymous mutation'

                    else:
                        print('Non-synonymous mutation')
                        substitution_type = 'Non-Synonymous mutation'

                else:
                    gene_sequence = genome_list[3][0][int(line[3]):int(line[2]) + 1]
                    gene_sequence = gene_sequence[::-1]
                    gene_start = int(line[2])
                    snp_pos = int(row[1])
                    snp_pos_adjusted = snp_pos - gene_start
                    snp_nt = str(row[2])

                    codon_number = int((snp_pos_adjusted) / 3) + 1

                    gene_snp = gene_sequence[:snp_pos_adjusted] + snp_nt + gene_sequence[snp_pos_adjusted + 1:]

                    codon_reference = gene_sequence[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]
                    codon_modified = gene_snp[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

                    print(gene_sequence)
                    print(snp_pos)
                    print(codon_reference)
                    print(codon_modified)

                    if translate(codon_reference) == translate(codon_modified):
                        print('Synonymous mutation')
                        substitution_type = 'Synonymous mutation'

                    else:
                        print('Non-synonymous mutation')
                        substitution_type = 'Non-Synonymous mutation'

            elif row[0] == 'Contig_3':
                
                if int(line[3]) > int(line[2]):
                    gene_sequence = genome_list[5][0][int(line[2]):int(line[3]) + 1]
                    gene_sequence = reverse_complement(gene_sequence)

                    gene_end = int(line[3])
                    snp_pos = int(row[1])
                    snp_pos_adjusted = snp_pos - gene_end
                    snp_nt = reverse_complement(str(row[2]))

                    codon_number = int((snp_pos_adjusted) / 3) + 1

                    gene_snp = gene_sequence[:snp_pos_adjusted] + snp_nt + gene_sequence[snp_pos_adjusted + 1:]

                    codon_reference = gene_sequence[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]
                    codon_modified = gene_snp[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

                    print(gene_sequence)
                    print(snp_pos)
                    print('reversed')
                    print(codon_reference)
                    print(codon_modified)

                    if translate(codon_reference) == translate(codon_modified):
                        print('Synonymous mutation')
                        substitution_type = 'Synonymous mutation'

                    else:
                        print('Non-synonymous mutation')
                        substitution_type = 'Non-Synonymous mutation'

                else:
                    gene_sequence = genome_list[5][0][int(line[3]):int(line[2]) + 1]
                    gene_sequence = gene_sequence[::-1]
                    gene_start = int(line[2])
                    snp_pos = int(row[1])
                    snp_pos_adjusted = snp_pos - gene_start
                    snp_nt = str(row[2])

                    codon_number = int((snp_pos_adjusted) / 3) + 1

                    gene_snp = gene_sequence[:snp_pos_adjusted] + snp_nt + gene_sequence[snp_pos_adjusted + 1:]

                    codon_reference = gene_sequence[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]
                    codon_modified = gene_snp[(codon_number - 1) * 3:(codon_number - 1) * 3 + 3]

                    print(gene_sequence)
                    print(snp_pos)
                    print(codon_reference)
                    print(codon_modified)

                    if translate(codon_reference) == translate(codon_modified):
                        print('Synonymous mutation')
                        substitution_type = 'Synonymous mutation'

                    else:
                        print('Non-synonymous mutation')
                        substitution_type = 'Non-Synonymous mutation'

            # write substitution type to output file
            output_file.write('\t' + str(substitution_type))

            # write gene annotation to output file
            output_file.write('\t' + str(line[4]))

            # write gene ID to output file
            output_file.write('\t' + str(line[1]))

            # write reference codon to output file
            output_file.write('\t' + str(codon_reference))

            # write modified codon to output file
            output_file.write('\t' + str(codon_modified))

            # write reference amino acid to output file
            output_file.write('\t' + str(translate(codon_reference)))

            # write modified amino acid to output file
            output_file.write('\t' + str(translate(codon_modified)))

            output_file.write('\n')

    # write intergenic to SNP type in output file
    if intergenic_counter == 0:
            # write contig number to output file
        output_file.write(row[0])

        # write SNP position within contig to output file
        output_file.write('\t' + row[1])
        output_file.write('\t' + 'Intergenic' + '\t-' + '\t-' + '\t-' + '\t-' + '\t-' + '\t-' + "\n")
        print('Intergenic')

# open and read input text file for creating bar plot
file = open("output_table")
data = csv.reader(file, delimiter="\t")

# set variables for containing data from output_table
contig_1_intergenic = 0
contig_1_synonymous = 0
contig_1_non_synonymous = 0
contig_2_intergenic = 0
contig_2_synonymous = 0
contig_2_non_synonymous = 0
contig_3_intergenic = 0
contig_3_synonymous = 0
contig_3_non_synonymous = 0

# skips header in text file
next(data)

# append data from text file to lists
for row in data:
    if row[0] == 'Contig_1' and row[2] == 'Intergenic':
        contig_1_intergenic+=1
    elif row[0] == 'Contig_1' and row[2] == 'Synonymous mutation':
        contig_1_synonymous+= 1
    elif row[0] == 'Contig_1' and row[2] == 'Non-Synonymous mutation':
        contig_1_non_synonymous+= 1
    elif row[0] == 'Contig_2' and row[2] == 'Intergenic':
        contig_2_intergenic+=1
    elif row[0] == 'Contig_2' and row[2] == 'Synonymous mutation':
        contig_2_synonymous+= 1
    elif row[0] == 'Contig_2' and row[2] == 'Non-Synonymous mutation':
        contig_2_non_synonymous+= 1
    elif row[0] == 'Contig_3' and row[2] == 'Intergenic':
        contig_3_intergenic+=1
    elif row[0] == 'Contig_3' and row[2] == 'Synonymous mutation':
        contig_3_synonymous+= 1
    elif row[0] == 'Contig_3' and row[2] == 'Non-Synonymous mutation':
        contig_3_non_synonymous+= 1


print(contig_1_intergenic)
print(contig_1_synonymous)
print(contig_1_non_synonymous)
print(contig_2_intergenic)
print(contig_2_synonymous)
print(contig_2_non_synonymous)
print(contig_3_intergenic)
print(contig_3_synonymous)
print(contig_3_non_synonymous)

print(contig_1_non_synonymous + contig_1_synonymous + contig_1_intergenic + contig_2_non_synonymous + contig_2_synonymous + contig_3_non_synonymous + contig_3_synonymous + contig_3_intergenic)

N = 3
ind = np.arange(N)  # the x locations for the groups
width = 0.27       # the width of the bars

fig = plt.figure()
ax = fig.add_subplot(111)

yvals = [contig_1_intergenic, contig_2_intergenic, contig_3_intergenic]
rects1 = ax.bar(ind, yvals, width, color='r')
zvals = [contig_1_synonymous,contig_2_synonymous,contig_3_synonymous]
rects2 = ax.bar(ind+width, zvals, width, color='g')
kvals = [contig_1_non_synonymous,contig_2_non_synonymous,contig_3_non_synonymous]
rects3 = ax.bar(ind+width*2, kvals, width, color='b')

ax.set_ylabel('Number of SNPs')
ax.set_xticks(ind+width)
ax.set_xticklabels( ('Contig 1', 'Contig 2', 'Contig 3') )
ax.legend( (rects1[0], rects2[0], rects3[0]), ('Intergenic', 'Synonymous', 'Non-Synonymous') )

def autolabel(rects):
    for rect in rects:
        h = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                ha='center', va='bottom')

#autolabel(rects1)
#autolabel(rects2)
#autolabel(rects3)

#plt.show()

# export final figure
plt.savefig("final-assignment-figure-mg", format="pdf")

print('\n' + '\x1b[6;30;42m' + 'Script complete. See output_table and final-assignment-figure-mg for output.' + '\x1b[0m')