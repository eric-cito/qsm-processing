function run_fastDipole_7T_HD_atlas(subj_ind_vec)
%
%
% Purpose: runs code for the subject
% index numbers specified in the subj_ind_vec. The subject index is given
% below.
%
%
%for example run_fastDipole([1,2]) would run
%on the first 2 subjects.
%

%path(path,'~/src/fsl/etc/matlab/');


% The Subject Index - cell array

% Controls = [3, 5:6, 8:11, 13:16, 31:32]

subj_ind = {

%       HD Subjects     //controls = [2,4:5,7:10,12:15,26], patients = [1,3,6,11,16:25];
% Subj numbers    Session numbers
%
{'b3375', 't8225'};     %1   
%{'b3372', 't8211'};     %    - no t1 .idf data
{'b3365', 't8183'};     %2                // Control //
{'b3299', 't7977'};     %3   
{'b3251', 't7801'};     %4                // Control //
{'b3236', 't7748'};     %5                // Control //
{'b3202', 't8015'};     %6
{'b3141', 't7435'};     %7   (slices=40)  // Control //
{'b3140', 't7434'};     %8   (slices=40)  // Control //
{'b3138', 't7433'};     %9  (slices=40)  // Control //
{'b3132', 't7418'};     %10  (slices=40)  // Control //
{'b3129', 't7413'};     %11  (slices=40)
{'b3128', 't7532'};     %12  (slices=40)  // Control //
{'b3094', 't7349'};     %13  (slices=40)  // Control //
{'b3089', 't7334'};     %14  (slices=40)  // Control //
{'b3079', 't7311'};     %15  (slices=40)  // Control // 
%{'b3064', 't7530'};     %  (slices=40) --- (same subj below)   
{'b3064', 't8195'};     %16  
{'b3049', 't7238'};     %17  (slices=40)
{'b3034', 't7194'};     %18  (slices=40)
%{'b3032', 't7187'};     %  (slices=40) --- (same subj as below) 
{'b3032', 't7993'};     %19  
{'b3028', 't7171'};     %20  (slices=40)
{'b3025', 't7158'};     %21  (slices=40)
{'b3023', 't7151'};     %22   -------------- magnitude image highly artifactual 
{'b3016', 't7126'};     %23  (slices=40)
%{'b3008', 't7088'};     %  (slices=40) --- (same subj as below)  
{'b3008', 't8003'};     %24
%{'b2996', 't7060'};     %  (slices=40) --- (same subj as below) 
{'b2996', 't8232'};     %25
{'b3400', 't8351'};     %26                // Control // 
%{'b3401', 't8350'};     %  (slices=46) ----- cmb phase image is unusable //Control //
};

subj_dir='/scratch/brain1/brain_work2/yicheng/3T_tof_swi/clare_qsm/code_for_ucsf_june2016';
%subj_dir='/data/chess1/cpoynton/MRdata_analysis/MGH/temp/hd/hd';

% Input filenames

fname_atlas = 'ROR_atlas'; % the atlas
fname_mask = 'ROR_phase_bet_mask_ero_sph3_man'; % brain mask
fname_mag = 'ROR_mag_cmb'; % magnitude image
fname_uphase = 'ROR_phase_cmb_kh_nan'; % unwrapped phase image

% Set parameters
params.atlas_thr = 0.96;        % Atlas threshold
params.atlas_flag = 1;          % set to 1 to use atlas, 0 to not use it
params.bo = 7.0;                % B0 field in Tesla
params.te = 16e-3;              % echo time in seconds
params.chi0 = 0.4;              % unitless, in ppm (do not change)
params.delta = -9.5;            % unitless, in ppm (do not change) 
params.num_iter = 300;          % number of CG iterations

%residual_all = zeros(size(subj_ind_vec,1),1);
residual_all = zeros(params.num_iter, size(subj_ind_vec,1));

subj_ind_vec = subj_ind_vec(:);

for i=1:(size(subj_ind_vec,1))
    
    ind = subj_ind_vec(i);
    
    string_pair = subj_ind{ind};

    subjnum_subj = char(string_pair{1})
    subjnum_session = char(string_pair{2}); 
    
    eval(['cd ', subjnum_subj, '/', subjnum_session])
    
    [nfm_final, nfm_dipole, atlas_thr, msk, x, scales, num_iter, residual] = fn_fastDipole_7T_HD_atlas(fname_atlas, fname_mask, fname_mag, fname_uphase, params);
    
    %residual_all(i) = residual;
    residual_all(:,i) = residual;

    %!mkdir KH_results_part1_300iters_phasebetmaskeroman_atlas_residual
    %cd KH_results_part1_300iters_phasebetmaskeroman_atlas_residual
    !mkdir results_part1
    cd results_part1
     
    save_avw(nfm_final,'FM_nfm_final','f',scales);
    save_avw(nfm_dipole,'nfm_dipole','f',scales);
    save_avw(x,'ext_sources','f',scales);
    save_avw(msk,'FM_msk','f',scales);
    save params.mat params fname_atlas fname_mag fname_uphase fname_mask
    save residual.mat residual
    disp('saved residual')

    
    eval(['cd ', subj_dir]);
end

save residual_all.mat residual_all
