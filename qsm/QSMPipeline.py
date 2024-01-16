import os
from typing import AnyStr, Any, Union, List

import SimpleITK as sitk

from .DicomToPhasePipeline import DicomToPhasePipeline
from .UnwrapPhasePipeline import UnwrapPhasePipeline
from .FileLocations import FileLocations
from PythonUtils.SkullStripper import SkullStripper
from .DicomGenerator import DicomGenerator

class BackgroundFieldRemovalAndDipoleInversionPipeline():

    def __init__(self, unwrappedPhase:sitk.Image, brainmask:sitk.Image, locs:FileLocations):
        self.phase = unwrappedPhase
        self.brainmask = brainmask
        self.locs = locs


    def Run(self) -> sitk.Image:
        raise Exception("Not implemented")


class QSMPipeline:

    def __init__(self, locs:FileLocations, useGPU=False):
        self.locs = locs
        self.useGPU = useGPU

    def Run(self) -> None:
        os.makedirs(self.locs.dir_out_top, exist_ok=True)

        (phase, mag, TEs) = DicomToPhasePipeline(self.locs).Run()

        print("Skull stripping")
        brainMask = self.GetOrCalcBrainmask(mag)
        unwrapped = UnwrapPhasePipeline(phase, mag, brainMask, TEs, self.locs).Run()

        print("Rescaling...")
        unwrapped = unwrapped * 

        print("Running background field removal and dipole inversion")
        dipoleInverted = BackgroundFieldRemovalAndDipoleInversionPipeline(unwrapped, brainMask, self.locs).Run()

        final = dipoleInverted

        print("Generating DICOMs")
        DicomGenerator(sitk.Cast(final * 1000, sitk.sitkInt16), self.locs.dir_dicoms_out, self.locs.loc_dicomHeader).Run()

    def GetOrCalcBrainmask(self, mag) -> sitk.Image:
        # The last echo has the least skull visible
        # and first seems to have the most signal
        # Have had issues using the last echo
        lastEcho = self.ExtractFirstInSeries(mag)
        brainMask = SkullStripper(lastEcho, useGPU=self.useGPU).CalcBrainmask(self.locs.loc_brainmask)
        return brainMask

    def ExtractFirstInSeries(self, image:sitk.Image) -> sitk.Image:
        '''
        Extracts the last frame in a series of 3D imaages
        '''
        size = list(image.GetSize())
        newSize = [size[0], size[1], size[2], 0]
        index = [0, 0, 0, 0] # Last index here for the series. e.g. size[3]-1
        return sitk.Extract(image, size=newSize, index=index)



        
        

        

