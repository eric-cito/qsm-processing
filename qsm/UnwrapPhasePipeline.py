from PythonUtils import AssertImage
from .FileLocations import FileLocations
from .Romeo import RunRomeo_sitk
import SimpleITK as sitk
import os
from typing import List


class UnwrapPhasePipeline:

    def __init__(self, phase:sitk.Image, mag:sitk.Image, brainmask:sitk.Image, TEs:List[float], locs:FileLocations):
        self.phase = phase
        self.magnitude = mag
        self.brainmask = brainmask
        self.TEs = TEs
        self.locs = locs
        AssertImage.AssertAreSameSize(self.phase, self.magnitude)
        AssertImage.AssertAreSameSize(self.phase, self.brainmask)


    def Run(self) -> sitk.Image:
        if os.path.exists(self.locs.phase_unwrapped):
            print("Unwrapped phase found. Generation skipped")
            return sitk.ReadImage(self.locs.phase_unwrapped)

        #self.EnsureOddSliceCount()

        unwrapped =  self.UnwrapPhase()

        return unwrapped


    def EnsureOddSliceCount(self):
        if self.phase.GetSize()[2] % 2 == 0:
            # Even number of slices
            self.phase = sitk.ConstantPad(self.phase, sitk.VectorInt32(0,0,1))
            self.magnitude = sitk.ConstantPad(self.magnitude, sitk.VectorInt32(0,0,1))
            self.brainmask = sitk.ConstantPad(self.brainmask, sitk.VectorInt32(0,0,1))


    def UnwrapPhase(self):
        if os.path.exists(self.locs.phase_unwrapped):
            print(self.locs.phase_unwrapped, "found. Unwrapping not re-performed")
            return sitk.ReadImage(self.locs.phase_unwrapped)
        else:
            return RunRomeo_sitk(self.phase, self.magnitude, self.TEs, self.brainmask, 
                                 None, # we have no use for the corrected wrapped image and it is huge
                                 self.locs.romeo_mask,
                                 self.locs.phase_unwrapped,
                                 self.locs.romeo_b0)