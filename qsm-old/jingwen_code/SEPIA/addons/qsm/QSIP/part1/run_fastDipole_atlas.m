function run_fastDipole_atlas(subj_ind_vec)
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

path(path,'~/src/fsl/etc/matlab/');


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

subj_ind_vec = subj_ind_vec(:);

for i=1:(size(subj_ind_vec,1))
    
    ind = subj_ind_vec(i);
    
    string_pair = subj_ind{ind};

    subjnum_rois = char(string_pair{1});
    subjnum_sph = char(string_pair{2}) 
    
    eval(['cd ', subjnum_sph])
    
    [nfm_final, nfm_dipole, atlas_thr, msk, x] = fn_fastDipole_1p5T_atlas(subjnum_sph);
    
    !mkdir results_part1
    cd results_part1
    
    scales = [0.9375, 0.9375, 2.5];
    
    save_avw(nfm_final,'FM_nfm_final','f',scales);
    save_avw(nfm_dipole,'nfm_dipole','f',scales);
    save_avw(x,'ext_sources','f',scales);
    save_avw(msk,'FM_msk','f',scales);
    %save_avw(atlas_thr,'FM_atlas_thr','f',scales);
    
    cd ../../
    
end

