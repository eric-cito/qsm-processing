from typing import Union
import os

class FileLocations:

    def __init__(self, dir_top):
        self.dir_top =  dir_top
        self.dir_dicoms_in = os.path.join(dir_top, 'qsm_dicoms/brain_ax_3d_swi_v1/')

        self.dir_out_top =os.path.join(dir_top, 'processed_QSM/')
        self.dir_raw_nii = os.path.join(self.dir_out_top, 'raw/')
        
        
        self.dir_phase_mag = os.path.join(self.dir_out_top, 'phase_mag/')

        self.loc_phase = os.path.join(self.dir_phase_mag, 'Phase.nii.gz')
        self.loc_magnitude = os.path.join(self.dir_phase_mag, 'Magni.nii.gz')


    def GetLoc(self, dir:str, echoNumber:Union[str,int], type:str, suffix:str="nii"):
        
        if int(echoNumber) < 10:
            echoNumber = "0" + str(echoNumber)
        else:
            echoNumber = str(echoNumber)

        loc = dir + "echo" + echoNumber + "_" +type+ "."+ suffix
        return loc