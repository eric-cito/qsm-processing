from typing import Optional
import os
import subprocess
import SimpleITK as sitk
from SimpleITK import Image
from typing import List
from PythonUtils.TempInMemory import TempDirInMemory
from PythonUtils.OSUtils import GetParentDirectory, AllExist, AllExistOrNone
from PythonUtils.IO import MoveGzSafe

def RunRomeo_sitk(phase:Image, 
                  magnitude:Image, 
                  echoTimes:List[float], 
                  mask:Optional[Image],
                  loc_saveTo_correctedPhase:Optional[str],
                  loc_saveTo_mask:Optional[str],
                  loc_saveTo_unwrapped:Optional[str],
                  loc_saveTo_B0:Optional[str]):

    if AllExistOrNone([loc_saveTo_correctedPhase, loc_saveTo_mask, loc_saveTo_unwrapped, loc_saveTo_B0]):
        print("Results found. Romeo unwrapping skipped.")
        return

    # Register and keep only the files we want
    # We work in a temporary directory to avoid file clashes
    with TempDirInMemory() as tempDir:
        dir = tempDir

        loc_phase = os.path.join(dir, "phase.nii")
        loc_magnitude = os.path.join(dir, "mag.nii")
        loc_mask = os.path.join(dir, "mask.nii")

        sitk.WriteImage(phase, loc_phase)
        sitk.WriteImage(magnitude, loc_magnitude)
        sitk.WriteImage(mask, loc_mask)

        RunRomeo(loc_phase, loc_magnitude, echoTimes, dir, loc_mask=loc_mask)

        # Generates the following
        print(os.listdir(dir))
        MoveGzSafe(os.path.join(dir, 'corrected_phase.nii'), loc_saveTo_correctedPhase, True)
        if loc_mask is None:
            MoveGzSafe(os.path.join(dir, 'mask.nii'), loc_saveTo_mask, True)
        MoveGzSafe(os.path.join(dir, 'unwrapped.nii'), loc_saveTo_unwrapped, True)
        MoveGzSafe(os.path.join(dir, 'B0.nii'), loc_saveTo_B0, True)


def RunRomeo(loc_phase:str, magnitude:str, echoTimes:List[float], dir_out:str, loc_mask:Optional[str] = None):
    binaryLocation = os.path.abspath(os.path.join(GetParentDirectory(__file__), "../romeo/romeo/bin/romeo"))
    args = [binaryLocation, "-p", loc_phase, "-m", magnitude, "-B", "-t", str(echoTimes), "-o",  dir_out]
    if loc_mask is not None:
        args.append("--mask")
        args.append(loc_mask)
    subprocess.run(args)