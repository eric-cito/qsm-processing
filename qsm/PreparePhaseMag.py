# Runs the dicom to phase pipeline
# supply the dicom directory for qsm, and the output directory
import sys
from .FileLocations import FileLocations
from .DicomToPhasePipeline import DicomToPhasePipeline

locs = FileLocations(sys.argv[1], sys.argv[2])

DicomToPhasePipeline(locs).Run()