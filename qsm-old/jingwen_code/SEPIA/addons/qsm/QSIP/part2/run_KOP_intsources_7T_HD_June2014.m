function [res] = run_KOP_intsources_7T_HD_June2014(subj_ind_vec)
%
%
% Purpose: runs code for the subject
% index numbers specified in the subj_ind_vec. The subject index is given
% below.
%
%
%for example [res] = run_KOP_intsources_7T_HD_June2014([1,2]) would run
%on the first 2 subjects, b3375 and b3365
%

%path(path,'~/src/fsl/etc/matlab/');

fsldir=getenv('FSLDIR');
fsldirmpath=sprintf('%s/etc/matlab/',fsldir);
path(path,fsldirmpath);
clear fsldir fsldirmpath


subj_ind = {

%       HD Subjects        //controls = [2,4:5,7:10,12:15,26], patients = [1,3,6,11,16:25];
% Subj numbers    Session numbers
%
{'b3375', 't8225'};     %1   
%{'b3372', 't8211'};     %               - no t1 .idf data
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
{'b0000', 't0000'};     %27 --- ART13618 volunteer, has no meaning
};


subj_dir='/scratch/brain1/brain_work2/yicheng/3T_tof_swi/clare_qsm/code_for_ucsf_june2016';
%subj_dir='/data/chess1/cpoynton/MRdata_analysis/MGH/temp/hd'

subj_ind_vec = subj_ind_vec(:);

objvals_allsubj = cell(size(subj_ind_vec,1), 1);

for i=1:(size(subj_ind_vec,1))
    
    ind = subj_ind_vec(i);
    
    string_pair = subj_ind{ind};

    subjnum_subj = char(string_pair{1})
    subjnum_session = char(string_pair{2}); 
    
    eval(['cd ', subjnum_subj, '/', subjnum_session])
    
    %ext_src_dir = 'KH_results_part1_300iters_phasebetmaskeroman_atlas_residual';
    ext_src_dir = 'results_part1';
    
       
    dims = [256, 256, 46]; %for low res "ROR" data
    
    [param chi1_est chi_est chi_est_ppm fmap_est f1 obj_vals noise_thr fmap_lp mask1 mask2 mask3 nfm_final fmap_acq_corr, scales] = KOP_intsources_nfm_corr3_7T_HD_June2014(subjnum_subj, subjnum_session, ext_src_dir, dims); 
    
    %!mkdir KH_results_part2_300iters_phasebetmaskeroman_atlas_residual_june2014
    
    %!mv susc* KH_results_part2_300iters_phasebetmaskeroman_atlas_residual_june2014
    
    %cd KH_results_part2_300iters_phasebetmaskeroman_atlas_residual_june2014
    
    !mkdir results_part2
    !mv susc* results_part2
    cd results_part2
    
    
    objvals_allsubj{i} = obj_vals;
    
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
    save obj_vals.mat obj_vals

    eval(['cd ', subj_dir]);
end

%save June2014_objvals_allsubj.mat objvals_allsubj

res = 1;
