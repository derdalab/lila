import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({'font.size': 7, 'font.family': 'arial', 'pdf.fonttype': 42})

pd.set_option('display.max_rows', None)


file_name = 'GL_VII_108_illumina_CHST1_and_CHST2_master_list.csv'
file_path = os.path.join(os.getcwd(), file_name)
dictionary = 'GF.xlsx'
number_sdbs = pd.read_excel(os.path.join(os.getcwd(), dictionary)).shape[0]
replicates = 3
targets = 3
figure = 'Fig3H.xlsx'

### PART 1

# 1 - Creates a master list of all SDBs and associated lectins
if not os.path.exists(file_path):

    master_list = pd.read_excel(os.path.join(os.getcwd(), dictionary)).iloc[:,[0,2,5]].set_index('SDB')

    # 2 - Opens the folder containing all sequencing files
    folder = os.path.join(os.getcwd(), 'Sequencing Files', file_name[:-16])
    files = os.listdir(folder)
    for file in files:
        # 3 - Creates panda dataframe from illumina file
        data = pd.read_csv(f'{folder}/{file}', sep=' ', header=9)
        primers = list(data.columns.values)[6:]

        # 4 - Parses through each primer combination (one replicate)
        for reads in primers:
            if file == '20230825-1728LJooBO-GL.txt' and reads not in ['R6F12_RN1RP1', 'R6F13_RN1RP2', 'R6F14_RN1RP3']:
                pass
            else:
                total_reads = int(data[reads][0])

                # 5 - Combine reads with same mindex and delete all mindex != 0
                y = 0
                for x in data['mindex']:
                    if x != 0:
                        data.loc[data.index == x, reads] += data[reads][y]
                    y += 1
                mindex_0 = data[data['mindex'] == 0]


                # 6 - Creates a list of desirable SDBs and corresponding reads in ppm
                df = []
                for x in master_list['Sequence']:

                    for y in mindex_0['Nuc']:
                        if x.upper() == y:
                            raw = (mindex_0.loc[mindex_0['Nuc'] == y][reads].tolist()[0])
                            ppm = ((mindex_0.loc[mindex_0['Nuc'] == y][reads].tolist()[0])/total_reads)*1000000
                            sequence = master_list.loc[master_list['Sequence'] == x].index[0]

                            df.append([raw, ppm, sequence])
                df = pd.DataFrame(df, columns=[reads[:-7]+'_Raw', reads[:-7]+'_CPM', 'SDB']).set_index('SDB')

                # 7 - Export data to master list


                master_list= pd.concat([master_list, df], axis=1)

                print(f'File {file}, read {reads} completed')
                print(master_list)

    # Separate column of raw reads from CPM and puts values in the correct order (i.e. input comes first followed
    # by the corresponding targets
    raw_cpm = [0,1] + list(range(2, master_list.shape[1], 2)) + list(range(3, master_list.shape[1], 2))
    master_list.iloc[:,raw_cpm].to_csv(file_path)

### PART 2


# 1 - Import data from step 1, fill Nan with 0 and convert to numpy array
data = pd.read_csv(file_path).fillna(0).to_numpy()

# 4 - Obtain Fold changes
input, wt, chst1, chst2 = np.mean(data[:,15:18 ], axis=1).reshape(number_sdbs,1), data[:, 18:21],data[:, 21:24], data[:, 24:27]
fc = np.concatenate((np.divide(wt, input, out=np.zeros_like(wt), where=input!=0).reshape(number_sdbs*replicates, 1),
                     np.divide(chst1, input, out=np.zeros_like(chst1), where=input!=0).reshape(number_sdbs*replicates, 1),
                     np.divide(chst2, input, out=np.zeros_like(chst2), where=input!=0).reshape(number_sdbs*replicates,1),
                     input.reshape(number_sdbs,1)))


cells = np.array(['Wt']*number_sdbs*replicates +
                 ['CHST1']*number_sdbs*replicates + ['CHST2']*number_sdbs*replicates
                 + ['Input']*number_sdbs).\
    reshape((number_sdbs*replicates*(targets+1)-(number_sdbs*2), 1))

lectin = np.array([lectin for lectin in data[:,2] for _ in range(replicates)]*targets +
                  [lectin for lectin in data[:,2]]).\
    reshape((number_sdbs*replicates*(targets+1))-(number_sdbs*2), 1)

df = pd.DataFrame(np.concatenate((cells, lectin, fc), axis=1), columns=['cells', 'lectin', 'fold'])
df_cells = df.loc[(df['cells'] != 'Input')]
df_input = df.loc[(df['cells'] == 'Input')]


def mbp_norm(df_cell):
    df2_cell = df_cell.copy()
    mbp_values = df2_cell.loc[df2_cell['lectin'].isin(['MBP'])]
    interval = int(mbp_values.shape[0]/targets)
    means = [mbp_values[0:interval]['fold'].mean(), mbp_values[interval:interval*2]['fold'].mean(),
             mbp_values[interval*2:interval*3]['fold'].mean()]

    new_fc = []
    n = 0
    counter = 0
    for x in df_cell['fold']:

        if counter < (number_sdbs*replicates-1):
            new_fc.append(x/means[n])
            counter += 1
        else:
            new_fc.append(x / means[n])
            n += 1
            counter = 0

    df_cell['fold'] = new_fc

    return df_cell


def quadratic_scale(df_cell):
    new_fc = []
    for x in df_cell['fold']:
        new_fc.append(x**(1/2))

    df_cell['fold'] = new_fc

    return df_cell

test = mbp_norm(df_cells)


export = pd.read_csv(file_path).fillna(0)
fc1, fc2, fc3, ix, nfc1, nfc2, nfc3  = df_cells[::3], df_cells[1::3], df_cells[2::3], df_input, \
                                 test[::3], test[1::3], test[2::3]

everything = [fc1, fc2, fc3, ix, nfc1, nfc2, nfc3]

for replicate in everything:
    for cell_type in set(replicate['cells']):
        desired_column = replicate.loc[replicate['cells'] == cell_type]['fold']
        desired_column = desired_column.rename(cell_type).reset_index(drop=True)
        export = pd.concat([export,desired_column], axis=1)

export.to_excel(figure, freeze_panes=(1,4))




if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(5.5,1.5))
    r1, r2, r3= test[::3], test[1::3], test[2::3]

    #sns.set(rc={'figure.figsize': (3.5, 3.5)}, style='ticks')
    order=['CHST1', 'CHST2', 'Wt']
    hue_order=['(Siglec-7)2','Siglec-7', 'Siglec-7R',
                                  'diCBM40-[290]', 'diCBM40-[145]', 'diCBM40-[73]', 'diCBM40-[26]', 'diCBM40-[3]', 'MBP',
                                   'Blank','Sc' ]
    sns.barplot(x='cells', y='fold', hue='lectin',
                order=order,
                       hue_order=hue_order,
                       data=test,
                       palette=[ '#44AD56','#71ad44', '#d6e2cd',
                                '#a06236', '#c67942', '#dc9969', '#e4af8a', '#f1d7c6','#506da5',
                                 '#ffffff','#c2c4c6'
                                ],
                       edgecolor='black',
                       ci=None, errcolor='black', capsize=0.01, errwidth=1.5

                       )


    sns.stripplot(x='cells', y='fold', hue='lectin', order=order,
                  hue_order=hue_order,
                         jitter=0.2, dodge=True, data=r1,
                         palette=['#000000'
                                  ], edgecolor='black', linewidth=1, size=3)

    sns.stripplot(x='cells', y='fold', hue='lectin', order=order,
                  hue_order=hue_order,
                         jitter=0.2, dodge=True, data=r2,
                         palette=['#ABABAB'
                                  ], edgecolor='black', linewidth=1, size=3)

    sns.stripplot(x='cells', y='fold', hue='lectin', order=order,
                  hue_order=hue_order,
                         jitter=0.2, dodge=True, data=r3,
                         palette=['#ffffff'
                                  ], edgecolor='black', linewidth=1, size=3)



    plt.yscale('log')
    plt.xlabel("Target")
    plt.ylabel("Fold Change (FC)")
    plt.ylim(0.05,10000)
    labels = [0.05, 0.1, 1, 10, 100, 1000, 10000]
    plt.yticks(labels, labels)
    plt.legend('')


    plt.savefig('GL-VII-108 chsts.pdf', transparent=True)

    plt.show()
