# ---------- N4 bias field correction and resampling ----------
import nibabel as nib
import numpy as np
import os
import glob
import SimpleITK as sitk
import pandas as pd
import traceback

dataset = r'C:\Users\Tracy\Desktop\image_mask'
image_S = 'image'
mask_S = 'ROI'
image_suffix = 'image.nii.gz'
mask_suffix = 'Merge.nii'

border = 120 # How many pixels remain around the tumor

scans = ['CT', 'MR'][1]
w_width = 400; w_center = 40 # Window width and window position

def calculate_voi_layers(filepath): # Calculate the maximum number of VOI layers
    img = nib.load(filepath)
    data = img.get_fdata()
    num_layers = data.shape[2]
    num_voi_layers = len([layer for layer in range(num_layers) if data[:,:,layer].any()])-1
    return num_voi_layers

def calculate_output_size(up_down_layers): # Set resampling voxel size
    output_size = [224,224,2*up_down_layers+1]
    return output_size

def crop_image(input_image_path, output_image_path, roi_image_path, output_roi_path, up_down_layers):
    input_image = sitk.ReadImage(roi_image_path)
    roi_arr = sitk.GetArrayFromImage(input_image)
    input_image = sitk.ReadImage(input_image_path)
    input_arr = sitk.GetArrayFromImage(input_image)
    if np.sum(input_arr.shape) != np.sum(roi_arr.shape):
        print('--------------------', input_image_path, ' Image ROI mismatch-----------------'); Stop

    if np.max(roi_arr) == 0:
        print('--------------------', input_image_path, ' ROI blank -------------------------'); Stop
    max_slices = np.argmax(np.sum(roi_arr, axis=(1, 2)))
    K = 1; Upsample = 1
    if '_1_' in input_image_path:
        Upsample = 0
    for ii in [0]:
        if K <= Upsample:
            slices = max_slices + list(range(-up_down_layers, up_down_layers + 1)) + ii
            try:
                cropped_arr = roi_arr[slices, :, :]; K = K + 1
            except:
                if max(slices) > min(roi_arr.shape)-up_down_layers:
                    slices = list(range(roi_arr.shape[0] - up_down_layers * 2 - 1, roi_arr.shape[0]))
                    cropped_arr = roi_arr[slices, :, :]
                else:
                    slices = list(range(0, up_down_layers * 2 + 1))
                    cropped_arr = roi_arr[slices, :, :]

            cropped_roi = sitk.GetImageFromArray(cropped_arr)
            cropped_roi.SetSpacing(input_image.GetSpacing())
            cropped_roi.SetOrigin((input_image.GetOrigin()[0], input_image.GetOrigin()[1],
                                   input_image.GetOrigin()[2] + slices[0] * input_image.GetSpacing()[2]))
            cropped_roi.SetDirection(input_image.GetDirection())
            sitk.WriteImage(cropped_roi, output_roi_path[0].split('Crop_' + roi_image_path)[0] + 'Crop_Max' + str(ii) + '_' +roi_image_path.split('\\')[-1])
            cropped_arr = input_arr[slices, :, :]
            cropped_image = sitk.GetImageFromArray(cropped_arr)
            cropped_image.SetSpacing(input_image.GetSpacing())
            cropped_image.SetOrigin((input_image.GetOrigin()[0], input_image.GetOrigin()[1],
                                     input_image.GetOrigin()[2] + slices[0] * input_image.GetSpacing()[2]))
            cropped_image.SetDirection(input_image.GetDirection())
            sitk.WriteImage(cropped_image, output_image_path[0].split('Crop_' + input_image_path)[0] + 'Crop_Max' + str(ii) + '_' +input_image_path.split('\\')[-1])

def resample_image(image_path, mask_path, output_size, image_interpolator=sitk.sitkLinear,
                   mask_interpolator=sitk.sitkNearestNeighbor):
    image = sitk.ReadImage(image_path)
    mask = sitk.ReadImage(mask_path)
    original_size = image.GetSize()
    original_spacing = image.GetSpacing()
    new_spacing = [original_spacing[i] * (original_size[i] / output_size[i]) for i in range(len(output_size))]
    resampler = sitk.ResampleImageFilter()
    resampler.SetSize(output_size)
    resampler.SetOutputOrigin(image.GetOrigin())
    resampler.SetOutputDirection(image.GetDirection())
    resampler.SetInterpolator(image_interpolator)
    resampler.SetOutputSpacing(new_spacing)
    resampled_image = resampler.Execute(image)
    resampler.SetInterpolator(mask_interpolator)
    resampled_mask = resampler.Execute(mask)
    return resampled_image, resampled_mask

def crop_image_with_mask(image_path, mask_path, output_image_path, output_roi_path, border):
    image = nib.load(image_path)
    mask = nib.load(mask_path)
    image_data = image.get_fdata()
    mask_data = mask.get_fdata()
    nonzero_coords = np.nonzero(mask_data)
    min_x, max_x = min(nonzero_coords[0]), max(nonzero_coords[0])
    min_y, max_y = min(nonzero_coords[1]), max(nonzero_coords[1])
    center = np.mean(nonzero_coords, axis=1)
    Max = max(max_x-min_x, max_y-min_y)/2
    cropped_image = image_data[(int(center[0]-Max-border)):(int(center[0]+Max+border)),
                               (int(center[1]-Max-border)):(int(center[1]+Max+border)), ]
    cropped_mask  = mask_data[(int(center[0]-Max-border)):(int(center[0]+Max+border)),
                              (int(center[1]-Max-border)):(int(center[1]+Max+border)), ]
    cropped_image_nifti = nib.Nifti1Image(cropped_image, image.affine, image.header)
    cropped_mask_nifti = nib.Nifti1Image(cropped_mask, mask.affine, mask.header)
    nib.save(cropped_image_nifti, str(output_image_path[0]))
    nib.save(cropped_mask_nifti, str(output_roi_path[0]))


def adjustMethod1(data_resampled, w_width, w_center):
    val_min = w_center - (w_width / 2)
    val_max = w_center + (w_width / 2)

    data_resampled[data_resampled < val_min] = val_min
    data_resampled[data_resampled > val_max] = val_max

    return data_resampled
os.chdir(dataset)
parent_path = os.path.dirname(dataset)
new_folder_name = "Crop_"+dataset.split('\\')[-1]
save_path = os.path.join(parent_path, new_folder_name)
output_path = os.path.join(parent_path, new_folder_name)
if not os.path.exists(save_path):
    os.makedirs(save_path)
else:
    print(f"{save_path} already exist")

Image_List = sorted(glob.glob('*' + image_suffix))
Mask_List = sorted(glob.glob('*' + mask_suffix))

if len(Image_List)-len(Mask_List) == 0:
    Error_list = pd.DataFrame({'file': [], 'error reason': []})
    for jj in range(0, len(Image_List)):
        image_path = Image_List[jj]; mask_path = Mask_List[jj]
        try:
            os.chdir(dataset)
            image_List = [Image_List[jj]]; mask_List = [Mask_List[jj]]
            for input_image_path, roi_image_path in zip(image_List, mask_List):
                print(jj+1, '/', len(Mask_List), ' ImageName =', input_image_path, 'MaskName =', roi_image_path, '\n')
                output_image_path = [output_path + '\\Crop_' + input_image_path]
                output_roi_path = [output_path + '\\Crop_' + roi_image_path]
                crop_image(input_image_path, output_image_path, roi_image_path, output_roi_path, calculate_voi_layers(roi_image_path))

            os.chdir(save_path)
            if scans == 'CT':
                for filename in os.listdir(save_path):
                    if "mage" in filename:
                        imagePath = os.path.join(save_path, filename)
                        input_image = sitk.ReadImage(imagePath)
                        roi_arr = sitk.GetArrayFromImage(input_image)
                        origin = input_image.GetOrigin(); direction = input_image.GetDirection(); space = input_image.GetSpacing()
                        roi_arr = adjustMethod1(roi_arr, w_width, w_center)
                        savedImg = sitk.GetImageFromArray(roi_arr)
                        savedImg.SetOrigin(origin)
                        savedImg.SetDirection(direction)
                        savedImg.SetSpacing(space)
                        sitk.WriteImage(savedImg, os.path.join(save_path, imagePath.split('\\')[-1]))

            image_List = glob.glob('*'+str(input_image_path))
            mask_List = glob.glob('*'+str(roi_image_path))
            for image_path, mask_path in zip(image_List, mask_List):
                output_image_path = [output_path + '\\' + image_path]
                output_roi_path = [output_path + '\\' + mask_path]
                crop_image_with_mask(image_path, mask_path, output_image_path, output_roi_path, border)

            if scans == 'MR':
                people = 0
                for i in os.listdir(save_path):
                    if "mage" in i.lower():
                        people = people + 1
                        imagePath = os.path.join(save_path, i)
                        input_image = sitk.ReadImage(imagePath)
                        mask_image = sitk.OtsuThreshold(input_image, 0, 1, 200)
                        input_image = sitk.Cast(input_image, sitk.sitkFloat32)
                        corrector = sitk.N4BiasFieldCorrectionImageFilter()
                        output_image = corrector.Execute(input_image, mask_image)
                        output_image = sitk.Cast(output_image, sitk.sitkInt16)
                        sitk.WriteImage(output_image, imagePath)

            for image_path, mask_path in zip(image_List, mask_List):
                resampled_image, resampled_mask = resample_image(image_path, mask_path, calculate_output_size(calculate_voi_layers(mask_path)))
                sitk.WriteImage(resampled_image, image_path)
                sitk.WriteImage(resampled_mask, mask_path)

        except:
            os.chdir(save_path)
            if os.path.exists(image_path):
                os.remove(image_path)
            if os.path.exists(mask_path):
                os.remove(mask_path)
            print(traceback.format_exc(limit=None, chain=True))
            Error_list.loc[jj, :] = [Image_List[jj], traceback.format_exc()]
            print('error patient:', Image_List[jj])
    os.chdir(dataset.split('\\'+dataset.split('\\')[-1])[0])
    Error_list.to_csv(dataset.split('\\')[-1]+'Error_list.csv', encoding='gbk', index=False)
    print('Feture Extration is Finish')
