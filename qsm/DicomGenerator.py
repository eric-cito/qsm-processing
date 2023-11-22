import SimpleITK as sitk
import tempfile
import json
import time
import os
import warnings
from typing import List, Tuple

class DicomGenerator:
    '''
    Writes the input image as DICOMS 
    '''

    def __init__(self, image:sitk.Image, dir_out:str, loc_dicomJSON:str):
        '''
        :image: The image to convert to DICOM
        :dir_out: The directory path to save to
        :loc_dicomJSON: Location of a JSON created by dcm2niix when the original image, from which this was derived, was converted from DICOM
        '''
        self.source = image
        self.dir_out = dir_out
        self.loc_dicomJSON = loc_dicomJSON

    def Run(self):
        
        tags = self.GetSharedTags()

        self.writer = sitk.ImageFileWriter()
        self.writer.KeepOriginalImageUIDOn()

        os.makedirs(self.dir_out, exist_ok=True)

        print(self.source.GetSize())
        for i in range(self.source.GetDepth()):
            self.WriteSlice(tags, i)

        # eagerly clean up in case it holds file handles or similar
        del self.writer

    def WriteSlice(self, tags, i):
        # Based on code here: https://simpleitk.readthedocs.io/en/v1.2.2/Examples/DicomSeriesReadModifyWrite/Documentation.html
        
        image_slice = self.source[:,:,i]
        
        # Tags shared by the series.
        for tag, value in tags:
            image_slice.SetMetaData(tag, value)

        # Slice specific tags.
        image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d")) # Instance Creation Date
        image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S")) # Instance Creation Time
        image_slice.SetMetaData("0020|0032", '\\'.join(map(str, self.source.TransformIndexToPhysicalPoint((0,0,i))))) # Image Position (Patient)
        image_slice.SetMetaData("0020|0013", str(i)) # Instance Number

        # Write to the output directory and add the extension dcm, to force writing in DICOM format.
        self.writer.SetFileName(os.path.join(self.dir_out, str(i)+'.dcm'))
        self.writer.Execute(image_slice)

    def GetSharedTags(self):
        '''
        Returns tags all dicoms will have
        '''

        inputProperties = self.ReadJSON()
        tags_to_copy = {
                "PatientName": "0010|0010",
                "PatientID": "0010|0020",
                "PatientBirthDate": "0010|0030",
                "StudyInstanceUID":"0020|000D",
                "StudyID":"0020|0010",
                "StudyDate":"0008|0020",
                "StudyTime":"0008|0030",
                "AccessionNumber":"0008|0050",
                "Modality":"0008|0060"
        }

        tags = []
        for key,value in tags_to_copy.items():
            if key in inputProperties:
                tags.append((value, inputProperties[key]))
            else:
                warnings.warn(key + " was not found in the original dicom and so is not in the output dicom")
        

        tags.append(("0008|103e", inputProperties["SeriesDescription"] + " Processed"))  # Series Description
        tags.append(("0008|0008","DERIVED\\SECONDARY")), # Image Type


        tags = tags + self.GetPixelTypeTag()


        # Image Orientation (Patient)
        # Based on code here: https://simpleitk.readthedocs.io/en/v1.2.2/Examples/DicomSeriesReadModifyWrite/Documentation.html
        orientation = self.source.GetDirection()
        # -- Reorder and format to DICOM standard
        orientation = [orientation[0], orientation[3], orientation[6],
                       orientation[1], orientation[4], orientation[7]]
        orientation = '\\'.join(map(str,orientation))
        tags.append(("0020|0037", orientation))# Image Orientation (Patient)

        modification_date = time.strftime("%Y%m%d")
        modification_time = time.strftime("%H%M%S")
        tags = tags + [
                        ("0008|0031",modification_time), # Series Time
                        ("0008|0021",modification_date), # Series Date
                        ("0020|000e", "1.2.826.0.1.3680043.2.1125."+modification_date+".1"+modification_time) # Series Instance UID
                    ]
                    
        return tags


    def GetPixelTypeTag(self) -> List[Tuple[str]]:
        if self.source.GetPixelID() == sitk.sitkInt16:
            return [
                ('0028|0100', '16'), # - bits allocated to 16
                ('0028|0101', '16'),# - bits stored to 16
                ('0028|0102', '15'),# - high bit to 15
                ('0028|0103','1') # - pixel representation to 1
            ]
        else:
            raise Exception("".join("Unsupported pixel type:", self.source.GetPixelIDTypeAsString()))


    def ReadJSON(self):
        with open(self.loc_dicomJSON) as f:
            data = json.load(f)
        return data
    
