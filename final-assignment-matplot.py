import numpy as np
import matplotlib.pyplot as plt
import csv

# open and read input text file
file = open("output_table")
data = csv.reader(file, delimiter="\t")

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
    if row[0] == 'Contig_2' and row[2] == 'Intergenic':
        contig_2_intergenic+=1
    elif row[0] == 'Contig_2' and row[2] == 'Synonymous mutation':
        contig_2_synonymous+= 1
    elif row[0] == 'Contig_2' and row[2] == 'Non-Synonymous mutation':
        contig_2_non_synonymous+= 1
    if row[0] == 'Contig_3' and row[2] == 'Intergenic':
        contig_3_intergenic+=1
    elif row[0] == 'Contig_3' and row[2] == 'Synonymous mutation':
        contig_3_synonymous+= 1
    elif row[0] == 'Contig_3' and row[2] == 'Non-Synonymous mutation':
        contig_3_non_synonymous+= 1

'''print(contig_1_intergenic)
print(contig_1_synonymous)
print(contig_1_non_synonymous)
print(contig_2_intergenic)
print(contig_2_synonymous)
print(contig_2_non_synonymous)
print(contig_3_intergenic)
print(contig_3_synonymous)
print(contig_3_non_synonymous)'''

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

ax.set_ylabel('Number of SNPS')
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

plt.show()

# export final figure
#plt.savefig("final-assignment-figure-mg", format="pdf")