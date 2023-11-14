
function [ wG ] = threshold_gradient( gm, mask, noise, percentage )
% Generate a binary image that acts as regularization weight, based on the gradient
% of the magnitude data.
%
% Based on the gradient_mask function in MEDI.
%
%
% Parameters:
% gm: gradient image (of the magnitude image, or R2* map, etc), as vector or magnitude.
% mask: binary 3D image that defines the ROI.
% noise: estimaded noise standard deviation in the complex signal
% percentage: fraction of the voxels to be considered as relevant.
%             Please see the following references for recommended values:
% Liu T, Xu W, Spincemaille P, Avestimehr AS and Wang Y. Accuracy of the morphology enabled
% dipole inversion (MEDI) algorithm for quantitative susceptibility mapping in MRI. IEEE Trans Med
% Imaging. 2012 Mar;31(3):816-24.
% Wang S, Chen W, Wang C, Liu T, Wang Y, Pan C, Mu K, Zhu C, Zhang X and Cheng J, Structure
% Prior Effects in Bayesian Approaches of Quantitative Susceptibility Mapping, BioMed Research
% International. 2016;2738231:10 p.
%
% Please note that not thresholding may yield better results:
% Carlos Milovic, Berkin Bilgic, Bo Zhao, Julio Acosta-Cabronero, and Cristian Tejos. 
% Spatially weighted regularization with Magnitude prior knowledge for QSM. 
% ESMRMB 2017, Barcelona, Spain.
%              
% Try instead: wG = max(gm,noise);
%
% Output:
% wG: binary image. 0 denotes a relevant gradient, 1 elsewhere.
%
% Created by Ildar Khalidov in 2010
% Last modified by Carlos Milovic in 2017.03.30

field_noise_level = noise;
denominator = sum(mask(:)==1);
wG = gm;
numerator = sum(wG(:)>field_noise_level);
if  (numerator/denominator)>percentage
    while (numerator/denominator)>percentage
        field_noise_level = field_noise_level*1.05;
        numerator = sum(wG(:)>field_noise_level);
    end
else
    while (numerator/denominator)<percentage
        field_noise_level = field_noise_level*.95;
        numerator = sum(wG(:)>field_noise_level);
    end
end

wG = (wG<=field_noise_level);



end

