import SimpleITK as sitk
import sys

loc_im = sys.argv[1]
loc_mask = sys.argv[2]
loc_saveTo = sys.argv[3]

image = sitk.ReadImage(loc_im)
mask = sitk.ReadImage(loc_mask)

masked_image = sitk.Mask(image, mask)

sitk.WriteImage(masked_image, loc_saveTo)
