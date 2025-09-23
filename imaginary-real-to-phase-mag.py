import SimpleITK as sitk
import sys

loc_real = sys.argv[1]
loc_imag = sys.argv[2]
loc_phase = sys.argv[3]

real = sitk.ReadImage(loc_real)
imaginary = sitk.ReadImage(loc_imag)

real = sitk.Cast(real, sitk.sitkFloat64)
imaginary = sitk.Cast(imaginary, sitk.sitkFloat64)


asComplex = sitk.RealAndImaginaryToComplex(real,imaginary)

# Negate every odd slice due to a GE's implementation where
# they get wrapping in the z dimension
for sliceNo in range(1, asComplex.GetSize()[-1], 2):
    asComplex[:,:,sliceNo] *= -1.0

phase = sitk.ComplexToPhase(asComplex)

# GE needs the phase to be negated
phase = phase * -1.0

# NB ITK can cause weirdness where the phase leaves the bounds 
# of -Pi to Pi due to rounding error, or something. This looks
# like striping 
# Subtract 2PI from all values >= PI
phase = sitk.Cast(phase, sitk.sitkFloat64)
phase += sitk.Cast(phase >= 3.14159265358979323846, sitk.sitkFloat64) * (-2 * 3.14159265358979323846)

sitk.WriteImage(phase, loc_phase)
print(f"Phase image written to {loc_phase}")