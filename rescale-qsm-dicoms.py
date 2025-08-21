#!/usr/bin/env python3
import os, sys, pydicom, numpy as np
from shutil import copyfile

if len(sys.argv) != 3:
    print("Usage: python rescale_qsm_dicoms.py <input_dir> <output_dir>")
    sys.exit(1)

input_dir = sys.argv[1]
output_dir = sys.argv[2]
os.makedirs(output_dir, exist_ok=True)

QSM_MIN = -0.35
QSM_MAX =  0.35
QSM_RANGE = QSM_MAX - QSM_MIN

STORED_MIN = 0
STORED_MAX = 65535
STORED_RANGE = STORED_MAX - STORED_MIN

RESCALE_SLOPE = QSM_RANGE / STORED_RANGE
RESCALE_INTERCEPT = QSM_MIN
WINDOW_CENTER = 32768
WINDOW_WIDTH  = 65535

print(f"Tagging DICOMs from: {input_dir}")
print(f"Writing to: {output_dir}")

for fname in sorted(os.listdir(input_dir)):
    if not fname.lower().endswith(".dcm"):
        continue

    path_in = os.path.join(input_dir, fname)
    path_out = os.path.join(output_dir, fname)

    dcm = pydicom.dcmread(path_in)

    dcm.BitsAllocated = 16
    dcm.BitsStored = 16
    dcm.HighBit = 15
    dcm.PixelRepresentation = getattr(dcm, "PixelRepresentation", 0)

    dcm.RescaleSlope = float(RESCALE_SLOPE)
    dcm.RescaleIntercept = float(RESCALE_INTERCEPT)
    dcm.WindowCenter = float(WINDOW_CENTER)
    dcm.WindowWidth  = float(WINDOW_WIDTH)

    if getattr(dcm, "SpecificCharacterSet", None) == "ISO-IR 100":
        dcm.SpecificCharacterSet = "ISO_IR 100"

    dcm.save_as(path_out)

print("Rescaled DICOMs saved.")
