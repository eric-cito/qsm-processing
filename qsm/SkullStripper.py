from typing import List
from .ImageGetter import ImageGetter
from .IO import WriteImageIfPathProvided
import SimpleITK as sitk
import os
import tempfile
import brainles_hd_bet.run

brainles_hd_bet.utils.folder_with_parameter_files = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")


class SkullStripper:
    def __init__(self, mrImgOrLocation):
        self.mrGetter = ImageGetter(mrImgOrLocation)


    def GetOrCalcBrainmask(self, loc_saveTo=None)  -> sitk.Image:
        '''
        Returns a mask of the brain
        '''

        if loc_saveTo is not None and os.path.exists(loc_saveTo):
            print(loc_saveTo, "found, skullstrip skipped")
            return sitk.ReadImage(loc_saveTo)

        loc_in = self.mrGetter.location
        writeAndDelete = loc_in is None
        if writeAndDelete:
            loc_in = tempfile.mktemp(suffix=".nii")
            sitk.WriteImage(self.mrGetter.GetImage(), loc_in)

        
        mask = self.RunHDBet_CPU(loc_in)
        
        WriteImageIfPathProvided(mask, loc_saveTo)

        # Cleanup
        if writeAndDelete:
            os.remove(loc_in)

        return mask
    
    def RunHDBet_CPU(self, input_files:[str, List[str]], mode:str = 'fast', tta=False):
        '''Helper method to call HD bet using the CPU. Returns the mask'''

        if mode != "fast" and mode != "accurate":
            raise Exception("Bad mode " + mode)
        
        loc_out = tempfile.mktemp(suffix="nii.gz")

        
        return brainles_hd_bet.run.run_hd_bet(input_files, 
                                        loc_out, 
                                        mode, 
                                        config_file=None, 
                                        device='cpu', 
                                        postprocess=True, 
                                        do_tta=tta, 
                                        keep_mask=False, 
                                        overwrite=True, 
                                        bet=False)[0]
