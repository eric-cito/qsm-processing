function []=recon_multiband_all_parfor(tnum, dicomfile)
%
% wrapper script around recon_multiband_grid for epi to run on all files on one machine
%
%addpath('/home/jlupo/svn/surbeck/common/multiband_epi/trunk/')

pfile=strcat(tnum, '_multiband_dti_noddi');
ref=strcat(tnum, '_multiband_dti_noddi_.ref.dat');
vrgf=strcat(tnum, '_multiband_dti_noddi_.vrgf.dat');
%[status, dcm]=system(sprintf('ls E*.DCM'));

disp(sprintf('Using dicom file %s', dicomfile))

process_cal(pfile, ref, vrgf);

npass = 98;
nslice = 20;

parfor pass = 1:npass
    for slice = 1:nslice
        recon_multiband_grid(pfile,ref,vrgf, dicomfile, pass, slice);
        disp(slice)
        disp(pass)
    end
end
system(sprintf('mkdir dicom'));
system(sprintf('mv *.dcm dicom/.'));

end

