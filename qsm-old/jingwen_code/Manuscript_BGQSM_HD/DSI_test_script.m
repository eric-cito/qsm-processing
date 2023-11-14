clc; clear;

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];

%% Loop through subjects and recon DTI files

for ii = 1:length(subjList)
    
    dataPath = sprintf('/working/lupolab/jingwen/004_HD_DSI/data/%s/%s', ...
        subjList{ii}, examList{ii});
    
    % do a little clean up
    cd(dataPath);
    mkdir('Archive');
    system('mv *.int2 Archive/');
    system('mv *.idf Archive');
    
    % convert dtitk to fib
    dtiFiles = dir(dataPath);
    dtiFiles(~contains({dtiFiles.name}, 'TOPUP_FA.nii.gz')) = [];
    if ~isempty(dtiFiles)
        dtiBase = strsplit(dtiFiles.name,'_TOPUP_FA.nii.gz');
        dtiBase = dtiBase{1};
        fprintf(' - DTI base : %s \n', dtiBase);
    end
    
    DSI_nii2fib([dataPath '/' dtiBase '_TOPUP_dtitk_aff_diffeo'], ...
        [dataPath '/' dtiBase '_dtitk']);
    
end

%% convert dtitk to fib - template

DSI_nii2fib('/working/lupolab/jingwen/004_HD_DSI/HC_mean_dtitk', ...
    '/working/lupolab/jingwen/004_HD_DSI/HC_mean_dtitk');