# copied from https://github.com/kamesy/QSM.jl/tree/main

try
    using QSM;
    using NIfTI;
catch
    import Pkg; 
    Pkg.add("NIfTI");
    Pkg.add("QSM");

    using QSM;
    using NIfTI;
end

# constants
y = 267.52      # gyromagnetic ratio
B0 = 3          # main magnetic field strength
TEs=[0.013, 0.015584, 0.018168, 0.020752, 0.023336, 0.02592, 0.028504, 0.031088, 0.033672, 0.036256]

uphas = niread("/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/processed_QSM/phase_mag/phase_corrected_unwrapped.nii.gz");
mask1 = !=(0).(niread("/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/processed_QSM/phase_mag/romeo_mask.nii.gz").raw);

# Voxel size in mm
# https://github.com/JuliaNeuroscience/NIfTI.jl
vsz = voxel_size(uphas.header) 

# convert units
@views for t in axes(uphas, 4)
    uphas.raw[:,:,:,t] .*= inv(B0 * y * TEs[t])
end

# remove non-harmonic background fields
print("Running VSharp")
fl, mask2 = vsharp(uphas.raw, mask1, (vsz[1], vsz[2], vsz[3])) # voxel size is an array but needs to be a tuple

asNii = NIVolume(uphas.header, uphas.extensions, fl)
niwrite("/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/processed_QSM/background-removed.nii.gz", asNii);

#niwrite("/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/processed_QSM/background-removed.nii.gz", mask2);

print("Done")