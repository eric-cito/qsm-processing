import sys
import os
from typing import List, AnyStr, Any, Union, Tuple
import subprocess
import dcm2niix
import SimpleITK as sitk
import json
from .OSUtils import ListDirWithFullPaths, DeleteAllIfFound, GetWithDifferentSuffix
from .FileLocations import FileLocations
from .SkullStripper import SkullStripper

class QSMPipeline:

    def __init__(self,locs:FileLocations):
        self.locs = locs

    def Run(self):
        os.makedirs(self.locs.dir_out_top, exist_ok=True)

        (phase, mag, TEs) = DicomToPhasePipeline(self.locs).Run()

        print("Skull stripping")
        brainMask = self.GetOrCalcBrainmask(mag)


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


class DicomToPhasePipeline:
    '''
    Converts dicoms, calculates magnitude and phase images
    '''

    def __init__(self,locs:FileLocations):
        self.locs = locs

    def Run(self, cleanUpIntermediates=True) -> Tuple[sitk.Image]:

        if (os.path.exists(self.locs.loc_phase) and
            os.path.exists(self.locs.loc_magnitude) and
            os.path.exists(self.locs.loc_TEs)):
            print("Magnitude and phase images found. Skipping generation")
            phase = sitk.ReadImage(self.locs.loc_phase)
            magnitude = sitk.ReadImage(self.locs.loc_magnitude)
            TEs = self.ReadTEFile()
        else:       
            os.makedirs(self.locs.dir_raw_nii, exist_ok=True)
            os.makedirs(self.locs.dir_phase_mag, exist_ok=True)

            [locs_mag, locs_imag, locs_real] = self.ConvertDicoms()

            print("Creating Phase images")
            phase = sitk.JoinSeries(self.CalcPhaseImages(locs_imag, locs_real))
            sitk.WriteImage(phase, self.locs.loc_phase)

            print("Creating Magnitude images")
            magnitude = self.ConcatenateImages(locs_mag)
            sitk.WriteImage(magnitude, self.locs.loc_magnitude)

            TEs = self.GetOrCreateTEFile(locs_imag)

            if cleanUpIntermediates:
                print("Cleaning up")
                DeleteAllIfFound(locs_imag)
                DeleteAllIfFound(locs_mag)
                DeleteAllIfFound(locs_real)

            
        print("Dicom to Phase Pipeline complete")
        return (phase, magnitude, TEs)


    def GetOrCreateTEFile(self, locs_imaginary:List[str]) -> List[float]:
        '''
        Writes TEs to a JSON array file
        '''

        def ReadTE(locJSON):
            with open(locJSON) as f:
                data = json.load(f)
                return data["EchoTime"]
        
        TEs = list()
        for loc in locs_imaginary:
            te =  ReadTE(GetWithDifferentSuffix(loc, 'json'))
            TEs.append(te)
        
        with open(self.locs.loc_TEs, mode='x') as f:
            f.write("[")
            f.write(", ".join([str(te) for te in TEs]))
            f.write("]\n")

        return [float(te) for te in TEs]

    def ReadTEFile(self):
        with open(self.locs.loc_TEs) as f:
            data = json.load(f)
        return data["EchoTime"]

    def CalcPhaseImages(self, locs_imag, locs_real):
        '''
        Calculates phase from imaginary and real.
        NB ITK can cause weirdness where the phase leaves the bounds 
        of -Pi to Pi due to rounding error, or something. This looks
        like striping 
        '''
        phaseImages = list()
        for iEcho in range(0,len(locs_imag)):
            phaseImages.append(self.GetOrCreatePhaseForEcho(locs_imag[iEcho], locs_real[iEcho]))
        return phaseImages

    def ConcatenateImages(self, locs_images):
        allIms = list()
        for loc in locs_images:
            allIms.append(sitk.ReadImage(loc))
        return sitk.JoinSeries(allIms)


    def ConvertDicoms(self) -> Tuple[str]:

    
        def GetNiftis() -> Tuple[str]:
            locs_mag = ListDirWithFullPaths(self.locs.dir_raw_nii, "_mag.nii", errorIfNotFound=False)
            locs_imag = ListDirWithFullPaths(self.locs.dir_raw_nii, "_imaginary.nii", errorIfNotFound=False)
            locs_real = ListDirWithFullPaths(self.locs.dir_raw_nii, "_real.nii", errorIfNotFound=False)
            return (locs_mag, locs_imag, locs_real)

        def MoveToCorrectLocation(directory, suffix):
            '''
            Renames values < 10 to 01, 02, etc so sorting works
            '''
            for loc in ListDirWithFullPaths(directory, "." + suffix):
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


    def GetOrCreatePhaseForEcho(self, loc_imaginary, loc_real):
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
                asComplex[:,:,sliceNo] *= -1.0

            phase = sitk.ComplexToPhase(asComplex)
            sitk.WriteImage(phase, locTo)
            return phase

        
        

        

