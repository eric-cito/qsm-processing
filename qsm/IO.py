from typing import AnyStr, Optional, List
import SimpleITK as sitk
import os


def WriteImageIfPathProvided(img:sitk.Image, path:Optional[AnyStr]) -> bool:
    '''
    Writes the image to file if a path is provided.
    Returns true if the file was written
    '''
    write = path is not None
    if write:
        sitk.WriteImage(img, path, useCompression=path.endswith(".gz"), imageIO="NiftiImageIO")
    return write
        

