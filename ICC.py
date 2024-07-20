#%% https://pingouin-stats.org/generated/pingouin.intraclass_corr.html
# # --------------------ICC calculate ------------------------
import pingouin as pg
import pandas as pd
import numpy as np
import os
OUTDATED_IGNORE = 1
np.seterr(invalid='ignore')

cc = 38
for k in (1, 2):
    dataDir = r'C:\Users\Tracy\Desktop\ICC'
    readerA = '/ReaderA-T1-1.csv'
    if k == 1:
        readerB = '/ReaderA-T1-2.csv'
    if k == 2:
        readerB = '/ReaderB-T1-1.csv'

    data_ = pd.read_csv(dataDir + readerA)
    data_1 = pd.read_csv(dataDir + readerA)
    data_2 = pd.read_csv(dataDir + readerB)
    print('data_1\n', data_1[1:3], '\n')
    print('data_2\n', data_2[1:3], '\n\n')

    data_ = data_.iloc[:, cc:]
    data_1 = data_1.iloc[:, cc:]
    data_2 = data_2.iloc[:, cc:]
    print(data_[1:3])

    # %%
    list = data_.columns.values.tolist()
    data_feature = pd.DataFrame({"feature": list})
    data_2_ = pd.DataFrame()

    # %%
    data_1.insert(0, "reader", np.ones(data_1.shape[0]))
    data_2.insert(0, "reader", np.ones(data_2.shape[0]) * 2)

    data_1.insert(0, "target", range(data_1.shape[0]))
    data_2.insert(0, "target", range(data_2.shape[0]))
    data = pd.concat([data_1, data_2])

    for feature in list:
        icc = pg.intraclass_corr(data=data, targets="target", raters="reader", ratings="%s" % feature)
        print('icc', icc[1:2])
        icc_2 = pd.DataFrame(icc[1:2])  # icc(2,1)

        data_2_ = pd.concat([data_2_, icc_2])
        data_2_.reset_index(drop=True, inplace=True)

    data_2 = pd.concat([data_feature, data_2_], axis=1)
    index_75 = (data_2['ICC'] > 0.75)
    data_icc75 = data_[data_.columns[index_75]]

    # %%  excel
    if k == 1:
        data_2.to_csv(dataDir + '/Intra_ICC_2_value.csv', encoding='utf_8_sig', index=False)
        data_icc75.to_csv(dataDir + '/Intra_ICC_2_0.75_feature.csv', encoding='utf_8_sig', index=False)
        print('Intra_ICC finished')
    else:
        data_2.to_csv(dataDir + '/Inter_ICC_2_value.csv', encoding='utf_8_sig', index=False)
        data_icc75.to_csv(dataDir + '/Inter_ICC_2_0.75_feature.csv', encoding='utf_8_sig', index=False)
        print('Inter_ICC finished')

