import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# Obtain FC values from the cells

import D_GL_VI_148_lila_density_neuraminidase, F_GL_VI_182_pbmc_illumina



plt.rcParams.update({'font.size': 7, 'font.family': 'arial', 'pdf.fonttype': 42})

density = D_GL_VI_148_lila_density_neuraminidase.test
pbmcs = F_GL_VI_182_pbmc_illumina.test



df = pd.concat([density, pbmcs])
print(df)
lectin = 'diCBM40-[73]'
order = ['CMAS_n', 'CMAS', 'WT_n', 'B-cells', 'NK', 'Monocytes', 'T-cells',
           'WT']
data = df.loc[df['lectin'] ==lectin].set_index('cells').loc[order].reset_index()
# x, y = data['cells'], data['fold']

#sns.set(rc={'figure.figsize': (2, 3.5)}, style='ticks')
fig, ax = plt.subplots(figsize=(4,2))
r1, r2, r3 = data[::3],data[1::3],data[2::3]

plot = sns.barplot(x='lectin', y='fold', hue='cells',
                    data=data,
                   palette=['#ffffff','#f2ebf4','#e6d7ea','#d9c3df','#cdb0d4',
                            '#c09dca','#b38abf','#a677b5','#9965aa','#8c52a0'],

                   edgecolor='black',
                   ci=None, errcolor='black', capsize=0.01, errwidth=1.5,
                   dodge=True

                   )
sns.stripplot(x='lectin', y='fold', hue='cells',
                         jitter=.2, dodge=True, data=r1,
                         palette=['#000000'
                                  ], edgecolor='black', linewidth=1, size=4)

sns.stripplot(x='lectin', y='fold', hue='cells',

                         jitter=0.2, dodge=True, data=r2,
                         palette=['#808080'
                                  ], edgecolor='black', linewidth=1, size=4)

sns.stripplot(x='lectin', y='fold', hue='cells',
                         jitter=0.2, dodge=True, data=r3,
                         palette=['#ffffff'
                                  ], edgecolor='black', linewidth=1, size=4)

plt.yscale('log')
plt.legend('')
plt.ylim(0.5, 1000)
labels = [0.5, 1, 10, 100, 1000]
plt.yticks(labels, labels)
plt.ylabel("Fold Change (FC)")
plt.savefig(f'{lectin}-{order[-2]}horizontally.pdf', transparent=True)

plt.show()





# Create the plot


"""
This is the xy graph
plt.figure(figsize=(15, 6))
plt.plot(x, y, marker='o', linestyle='-', color='b', label='fold')
plt.xlabel('Cell Types')
plt.ylabel('Fold Change (FC)')
plt.yscale('log')
plt.legend('')
plt.xticks(rotation=45)  # Rotate the labels for better visibility
labels = [1, 10, 100, 1000]
ticks = labels
plt.yticks(ticks, labels)

plt.savefig('cbm50.pdf', transparent=True)
plt.show()


"""

