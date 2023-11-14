from . import RunHDBet
from .ImageGetter import ImageGetter
from .IO import WriteImageIfPathProvided
import SimpleITK as sitk
import os
import tempfile


class SkullStripper:
    def __init__(self, mrImgOrLocation):
        self.mrGetter = ImageGetter(mrImgOrLocation)


    def CalcBrainmask(self, loc_saveTo=None)  -> sitk.Image:
        '''
        Returns a mask of the brain
        '''

        if loc_saveTo is not None and os.path.exists(loc_saveTo):
            print(loc_saveTo, "found, skullstrip skipped")
            return sitk.ReadImage(loc_saveTo)

        loc_in = self.mrGetter.location
        writeAndDelete = loc_in is None
        if writeAndDelete:
            loc_in = tempfile.mktemp(suffix="nii")
            sitk.WriteImage(self.mrGetter.GetImage(), loc_in)

        
        mask = RunHDBet.RunHDBet_CPU(loc_in)
        
        WriteImageIfPathProvided(mask, loc_saveTo)

        # Cleanup
        if writeAndDelete:
            os.remove(loc_in)

        return mask