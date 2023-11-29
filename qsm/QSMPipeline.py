import os
from typing import AnyStr, Any, Union, List

import SimpleITK as sitk

from .DicomToPhasePipeline import DicomToPhasePipeline
from .UnwrapPhasePipeline import UnwrapPhasePipeline
from .FileLocations import FileLocations
from .SkullStripper import SkullStripper
from .DicomGenerator import DicomGenerator

class BackgroundFieldRemovalAndDipoleInversionPipeline():

    def __init__(self, unwrappedPhase:sitk.Image, brainmask:sitk.Image, locs:FileLocations):
        self.phase = unwrappedPhase
        self.brainmask = brainmask
        self.locs = locs


    def Run(self) -> sitk.Image:
        raise Exception("Not implemented")


class QSMPipeline:

    def __init__(self,locs:FileLocations):
        self.locs = locs

    def Run(self):
        os.makedirs(self.locs.dir_out_top, exist_ok=True)

        (phase, mag, TEs) = DicomToPhasePipeline(self.locs).Run()

        print("Skull stripping")
        brainMask = self.GetOrCalcBrainmask(mag)
        unwrapped = UnwrapPhasePipeline(phase, mag, brainMask).Run()

        dipoleInverted = BackgroundFieldRemovalAndDipoleInversionPipeline(unwrapped, brainMask, self.locs).Run()

        final = dipoleInverted

        DicomGenerator(sitk.Cast(final * 1000, sitk.sitkInt16), self.locs.dir_dicoms_out, self.locs.loc_dicomHeader).Run()

    def GetOrCalcBrainmask(self, mag):
        # The last echo has the least skull visible
        # so probably will do best in general
        # Though the brain looks strange.
        # If having issues, change to the first echo
        # where the skull and brain may look more like
        # the NN is expecting
        lastEcho = self.ExtractLastInSeries(mag)
        brainMask = SkullStripper(lastEcho).GetOrCalcBrainmask(self.locs.loc_brainmask)
        return brainMask

    def ExtractLastInSeries(self, image:sitk.Image):
        '''
        Extracts the last frame in a series of 3D imaages
        '''
        size = list(image.GetSize())
        newSize = [size[0], size[1], size[2], 0]
        index = [0, 0, 0, size[3]-1]
        return sitk.Extract(image, size=newSize, index=index)



        
        

        

