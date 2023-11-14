% Re/Imag images with slice to slice phase chopping
real1=dicomread('Image_00829.dcm');real1_offset=real1(1,1);real1=real1-real1_offset;
imag1=dicomread('Image_00830.dcm');imag1_offset=imag1(1,1);imag1=imag1-imag1_offset;
complex1=complex(double(real1),double(imag1));phase1=angle(complex1);figure;imagesc(phase1);colormap('gray');colorbar;

real2=dicomread('Image_00841.dcm');real2_offset=real2(1,1);real2=real2-real2_offset;
imag2=dicomread('Image_00842.dcm');imag2_offset=imag2(1,1);imag2=imag2-imag2_offset;
complex2=complex(double(real2),double(imag2));phase2=angle(complex2);figure;imagesc(phase2);colormap('gray');colorbar;

real3=dicomread('Image_00853.dcm');real3_offset=real3(1,1);real3=real3-real3_offset;
imag3=dicomread('Image_00854.dcm');imag3_offset=imag3(1,1);imag3=imag3-imag3_offset;
complex3=complex(double(real3),double(imag3));phase3=angle(complex3);figure;imagesc(phase3);colormap('gray');colorbar;

real4=dicomread('Image_00865.dcm');real4_offset=real4(1,1);real4=real4-real4_offset;
imag4=dicomread('Image_00866.dcm');imag4_offset=imag4(1,1);imag4=imag4-imag4_offset;
complex4=complex(double(real4),double(imag4));phase4=angle(complex4);figure;imagesc(phase4);colormap('gray');colorbar;
diff1=phase1-phase2;figure;imagesc(diff1);colormap('gray');colorbar;
diff2=phase2-phase3;figure;imagesc(diff2);colormap('gray');colorbar;
diff3=phase3-phase4;figure;imagesc(diff3);colormap('gray');colorbar;

% Re/Imag images with slice to slice phase chopping removed
complex1=-1*complex(double(real1),double(imag1));phase1=angle(complex1);figure;imagesc(phase1);colormap('gray');colorbar;
complex2=complex(double(real2),double(imag2));phase2=angle(complex2);figure;imagesc(phase2);colormap('gray');colorbar;
complex3=-1*complex(double(real3),double(imag3));phase3=angle(complex3);figure;imagesc(phase3);colormap('gray');colorbar;
complex4=complex(double(real4),double(imag4));phase4=angle(complex4);figure;imagesc(phase4);colormap('gray');colorbar;
diff1=phase1-phase2;figure;imagesc(diff1);colormap('gray');colorbar;
diff2=phase2-phase3;figure;imagesc(diff2);colormap('gray');colorbar;

diff3=phase3-phase4;figure;imagesc(diff3);colormap('gray');colorbar;