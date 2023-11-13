function res = mtimes(K, seg)
% res = mtimes(K, seg)
%
% seg:  is the segmentation (chi1 labels) image to be run through the foward model. 
%
% K: an object of type KOP which contains all the parameters of the foward model
%
% this code is essentially equivalent to fft_convolve.m

  %Lorentz correction
  bz = K.b0z * seg./(3+K.chi0).*(1+K.chi0);
    
  seg = K.b0z * seg;
  
  s = size(K.kernel);
  
  x = zeros(s);
  
  s2 = size(seg);
  
  x(1:s2(1),1:s2(2),1:s2(3)) = seg(1:s2(1),1:s2(2),1:s2(3));
  
  x_tr = fftn(x);
   
  x_tr_re = real(x_tr);
  x_tr_im = imag(x_tr);
  
  k_tr = fftn(K.kernel);
   
  k_tr_re = real(k_tr);
  k_tr_im = imag(k_tr);
  
  mult_re = x_tr_re.*k_tr_re-x_tr_im.*k_tr_im;
  mult_im = x_tr_re.*k_tr_im+x_tr_im.*k_tr_re;
  
  mult = mult_re + i*mult_im;
  
  conv_result =  ifftn(mult);
  
  conv_result = conv_result(s2(1):s(1),s2(2):s(2),s2(3):s(3));
  
  bz = bz - conv_result;
  
  %scale 
  res = bz .* (K.delta/(1+K.chi0));
  