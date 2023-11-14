%function [res] = run_KOP_intsources(subj_ind_vec)
%
%
% Purpose: runs code for the subject
% index numbers specified in the subj_ind_vec. The subject index is given
% below.
%
%
%for example [res] = run_fastDipole([1,2]) would run
%on the first 2 subjects, 0013 and 0213
%

%path(path,'~/src/fsl/etc/matlab/');

% The Subject Index - cell array

subj_ind = {

%       Young Subjects
% ROI numbers    Sphase numbers
%
{'1994_2738',	'0013_12192007'};       % 1
{'2087_2800', 	'0213_02212008'};       % 2
{'1993_2776', 	'0012_12192007'};       % 3
{'1978_2657', 	'0156_02062008'};       % 4
{'1805_2754', 	'0037_01042008'};       % 5
{'1795_2742', 	'0152_02052008'};       % 6
{'2061_2784', 	'0162_02072008'};       % 7
{'2094_2810', 	'0066_01112008'};       % 8
{'2165_2777', 	'0157_02062008'};       % 9 
{'2020_2768', 	'0011_12192007'};       % 10
{'2156_2880', 	'0014_12192007'};       % 11


%       Elderly subjects

%ROIs               Sphase

{'1936_2669',	  '0041_01072008'};     % 12
{'1893_2711',	  '0065_01112008'};     % 13
{'2047_2773',	  '0091_01182008'};     % 14
{'2144_2858',	  '0121_01252008'};     % 15
{'1968_2707',	  '0173_02082008'};     % 16
{'2139_2847',	  '0234_02272008'};     % 17
{'2030_2765',	  '0036_01042008'};     % 18 
{'1984_2714',	  '0057_01112008'};     % 19
{'1827_2716',	  '0070_01142008'};     % 20
{'1825_2713',	  '0106_01222008'};     % 21
{'2084_2794',	  '0172_02082008'};     % 22
{'2153_2875',	  '0223_02252008'};     % 23

};

%subj_ind_vec = subj_ind_vec(:);

%for i=1:(size(subj_ind_vec,1))
    
%     ind = subj_ind_vec(i);
%     
%     string_pair = subj_ind{ind};

    subjnum_rois = '';
    subjnum_sph = '';
    
%    eval(['cd ', subjnum_sph])
    
    ext_src_dir = 'results_part1';
       
    dims = [256, 256, 36];
    
    [param chi1_est chi_est chi_est_ppm fmap_est f1 obj_vals noise_thr fmap_lp mask1 mask2 mask3 nfm_final fmap_acq_corr, scales] = KOP_intsources_nfm_corr3(subjnum_sph, subjnum_rois, ext_src_dir, dims); 
    
    mkdir results_part2_chi0
    cd results_part2_chi0
    
    
    %save_avw(chi1_est,'FM_chi1','f',scales);
    save_avw(chi_est,'FM_chi_rel','f',scales);
    save_avw(chi_est_ppm,'FM_chi_ppm_rel','f',scales);
    save_avw(mask1,'FM_mask1','f',scales)
    save_avw(mask2, 'FM_mask2','f',scales)
    save_avw(mask3,'FM_mask3','f',scales)
    save_avw(nfm_final,'nfm_final','f',scales);
    save_avw(fmap_acq_corr,'fmap_acq_corr','f',scales);
    save_avw(fmap_lp,'fmap_lp','f',scales);
    %save results_intsources.mat mean_susc_rel rois noise_thr
    
    %cd ../../../
    
%end

res = 1;
