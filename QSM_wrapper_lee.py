import sys
import os
from typing import List, AnyStr, Any, Union, Tuple
import subprocess
import dcm2niix
from OSUtils import ListDirWithFullPaths
import SimpleITK as sitk
# QSM_processing 
# Written by Jingwen
# Modified by Melanie

# %% add path

# warning('off');
# addpath(genpath('/home/lreid/qsm/qsm/jingwen_code'));
# addpath('/data/morrison/scripts/matlab_snippets/nifti_tools');
# addpath('/home/lreid/qsm');

# %% set up data paths


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


class QSMPipeline:

    def __init__(self,locs:FileLocations):
        self.locs = locs

    def Run(self):
        os.makedirs(self.locs.dir_out_top, exist_ok=True)
        [locs_mag, locs_imag, locs_real] = self.ConvertDicoms()

        noEchoes = len(locs_mag)

        phaseImages = list()

        for iEcho in range(0,noEchoes):
            phaseImages.append(self.CreatePhaseForEcho(locs_imag[iEcho], locs_real[iEcho]))
        

    def ConvertDicoms(self) -> Tuple[str]:

    
        def GetNiftis() -> Tuple[str]:
            locs_mag = ListDirWithFullPaths(self.locs.dir_raw_nii, "_mag.nii")
            locs_imag = ListDirWithFullPaths(self.locs.dir_raw_nii, "_imaginary.nii")
            locs_real = ListDirWithFullPaths(self.locs.dir_raw_nii, "_real.nii")
            return (locs_mag, locs_imag, locs_real)

        def MoveToCorrectLocation(directory, suffix):
            '''
            Renames values < 10 to 01, 02, etc so sorting works
            '''
            for loc in ListDirWithFullPaths(directory, str.join(".",suffix)):
                number = self.ExtractEchoNumberFromFilename(loc)
                type = loc.split("_")[-1].split(".")[0]
                destination = self.locs.GetLoc(directory, number, type, suffix)

                if loc != destination:
                    os.rename(loc, destination)

        [locs_mag, locs_imag, locs_real] = GetNiftis()

        if (len(locs_mag) != 0) and (len(locs_imag) == len(locs_mag)) and (len(locs_real) == len(locs_mag)):
            print("Raw niftis found. Dicom conversion skipped")
            return (locs_mag, locs_imag, locs_real)
        
        dcm2niix.main(["-m", "n", "-f", "echo%e_%z", "-o", self.locs.dir_raw_nii,"-w", "0", self.locs.dir_dicoms_in])

        # Files will be named
        # magnitude: echoX_.nii
        # real: echoX_real.nii
        # imaginary: echoX_imaginary.nii
        # where X is the echo number

        # Fix up file names
        # -- mag files don't have a type listed
        for loc in ListDirWithFullPaths(self.locs.dir_raw_nii, "_.nii"):
            os.rename(loc, loc.replace("_.nii", "_mag.nii"))
        

        for loc in ListDirWithFullPaths(self.locs.dir_raw_nii,"_.json"):
            os.rename(loc, loc.replace(".json", "mag.json"))
        
        MoveToCorrectLocation(self.locs.dir_raw_nii, "nii")
        MoveToCorrectLocation(self.locs.dir_raw_nii, "json")

        return GetNiftis()


    def ExtractEchoNumberFromFilename(self, loc:str) -> int:
        fn = os.path.split(loc)[-1]
        return int(fn[4:].split("_")[0])


    def CreatePhaseForEcho(self, loc_imaginary, loc_real):
        '''
        Gets the phase image if existing, else creates it
        '''

        echoNumber = self.ExtractEchoNumberFromFilename(loc_imaginary)
        locTo = self.locs.GetLoc(self.locs.dir_phase_mag, echoNumber, "phase")

        if os.path.exists(locTo):
            return sitk.ReadImage(locTo)
        else:
            imaginary = sitk.ReadImage(loc_imaginary,sitk.sitkFloat64)
            real = sitk.ReadImage(loc_real,sitk.sitkFloat64)

            asComplex = sitk.RealAndImaginaryToComplex(real,imaginary)

            # Negate every odd slice
            for sliceNo in range(1, asComplex.GetSize()[-1], 2):
                asComplex[:,:,sliceNo] *= -1

            phase = sitk.ComplexToPhase(asComplex)
            sitk.WriteImage(phase, locTo)
            return phase

        
        

        

locs = FileLocations('/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/'); # must end with /s)
pipeline = QSMPipeline(locs)
pipeline.Run()