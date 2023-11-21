import os
import subprocess
import SimpleITK as sitk
from SimpleITK import Image
from typing import List
import tempfile

def RunRomeo_sitk(phase:Image, magnitude:Image, echoTimes:List[float], loc_saveTo:str):

    if os.path.exists(loc_saveTo):
        print("Result found. Romeo unwrapping skipped.")
        return

    # Register and keep only the files we want
    # We work in a temporary directory to avoid file clashes
    tempDir = tempfile.mkdtemp()
    try:
        loc_phase = os.path.join(tempDir, "phase.nii")
        loc_magnitude = os.path.join(tempDir, "mag.nii")

        sitk.WriteImage(phase, loc_phase, phase)
        sitk.WriteImage(phase, loc_magnitude, magnitude)


        RunRomeo(loc_phase, loc_magnitude, echoTimes, tempDir)

        # shutil.move(pipeline.loc_mrInCTSpace, loc_saveTo)
    finally:
        # Clean up
        raise Exception("Get file out above")
        shutil.rmtree(tempDir)    

def RunRomeo(loc_phase:str, magnitude:str, echoTimes:List[float], dir_out:str):
    subprocess.Run(["../romeo/romeo/bin/romeo", oc_phase, "-m", magnitude, "-B", "-t", str(echoTimes), "-o",  dir_out])