#https://pyradiomics.readthedocs.io/en/latest/features.html #
# ------------- Radiomics features extraction -------------------
from radiomics import featureextractor
import pandas as pd
import os
import glob
import SimpleITK as sitk
import traceback

os.chdir(r'C:\Users\Tracy\Desktop\Feature_extract')
num_label = 1
film_list = os.listdir()
film_list = sorted(film_list)
cc = 0
film = pd.DataFrame(film_list[::2])
# ------------------- 1. image mask test ---------------
for ii in range(cc * 2, len(film_list), 2):
    filename = film_list[ii]
    K0 = 0
    K = 3
    print('filename', ii, ':', str(filename[K0:K]))
    start = str(filename[K0:K])
    image_List = glob.glob(start + '_image.nii.gz')
    print('image_List_test', image_List)
    if any(image_List)==False:
       print('the image error', input())

    mask_List = glob.glob(start + '_Merge.nii')
    print('mask_List_test', mask_List)
    if any(mask_List)==False:
       print('the mask error', input())

# ---------------- 2. feature extraction ------------------
df = pd.DataFrame()
Error_list = pd.DataFrame({'error file': [], 'error reason': []})
for ii in range(cc * 2, len(film_list), 2):
    filename = film_list[ii]
    print('filename', ii, ':', str(filename[K0:K]))
    start = str(filename[K0:K])
    image_List = glob.glob(start + '_image.nii.gz')
    print('image_List', image_List)
    mask_List = glob.glob(start + '_Merge.nii')
    print('mask_List', mask_List)
    settings = {'binWidth': 25, 'resampledPixelSpacing': [1, 1, 1],
                'interpolator': sitk.sitkBSpline, 'normalize': True, 'normalizeScale': 100}
    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.enableImageTypes(Original={}, Square={}, Wavelet={})
    for imageName, maskName in zip(image_List, mask_List):
        try:
            featureVector = extractor.execute(imageName, maskName, label=num_label)
            df_add = pd.DataFrame.from_dict(featureVector.values()).T
            df_add.columns = featureVector.keys()
            df = df.append(pd.concat([pd.DataFrame({'p_ID': filename}, index=[0]), df_add], axis=1))
            print('process：', ii/2+1, '/', len(film_list) / 2, '=', round((ii + 2) / len(film_list), 3), '\n\n')
        except:
            print(traceback.format_exc(limit=None, chain=True))
            Error_list.loc[ii, :] = [filename, traceback.format_exc()]
            print('error patient id:', filename)

df.to_csv('Radiomic results.csv', encoding='gbk', index=False)
Error_list.to_csv('Error_list.csv', encoding='gbk', index=False)
print('Feture Extration is Finish')
