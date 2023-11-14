clear; clc;
warning('off');

%% add path

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList.xlsx');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status];
AgeList = [T.age];

HCind = strcmp(statusList,'HC');
subjList = subjList(HCind);
examList = examList(HCind);
AgeList = AgeList(HCind);
statusList = statusList(HCind);

DCMout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/DICOM';
mkdir(DCMout_root);

QSMfile_list = {'QSM_iLSQR_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho'};

subjFolders = dir('/data/7T_hunt');

%% randomly select 10 subjects

rng(1);

ind = randsample(length(subjList),10);
subjSel = subjList(ind);
examSel = examList(ind);
AgeSel = AgeList(ind);

figure;
histogram(AgeList,[20:10:70]); hold on;
histogram(AgeSel,[20:10:70]);

%% randomly assign the case number and record the correspondance

rng(1);

Order = randperm(100)';

subj = reshape(repmat(subjSel,[1 10])',100,1);
QSMfile = repmat(QSMfile_list',[10 1]);

T = table(subj,QSMfile,Order);

save([DCMout_root '/../VisualList.mat'], 'T');

%% Loop through subjects

for ii = 1:length(subjSel)
    
    %     if ~strcmp(statusList{ii},'HC') && ~strcmp(statusList{ii},'EM') && ~strcmp(statusList{ii},'PM')
    %         continue
    %     end
    
    fprintf('%s is %s \n', subjSel{ii}, statusList{ii});
    
    examPath = ['/data/7T_hunt/' subjSel{ii} '/' examSel{ii}];
    
    if exist([examPath '/swan_qsm/HDBET_allQSM/Phase.nii.gz'], 'file') == 2
        
        for qq = 1:length(QSMfile_list)
            
            ind = (ii-1)*length(QSMfile_list) + qq;
            examID = [subjSel{ii} '_' examSel{ii} '_' num2str(qq)];
            
            if exist(sprintf('%s/%s/IM%04i.dcm',DCMout_root,examID,1),'file') == 2
                % continue;
            end
            
            % load nifti file
            QSMfile_root = [examPath '/swan_qsm/HDBET_allQSM/'];
            QSMpath = [QSMfile_root '/' QSMfile_list{qq} '.nii.gz'];
            nii = load_untouch_nii(QSMpath);
            img = double(nii.img);
            
            if qq == 9
                img = img/0.5684;
            end
            
            img = img + 0.15;
            img(img < 0) = 0; img(img > 0.3) = 0.3;
            img = flip(flip(rot90(int16(img/0.3*4096))),3);
            
            mkdir([DCMout_root '/' examID '_ax']);
            
            uid = dicomuid;
            meta.Modality = 'MR';
            meta.SeriesNumber = ii*100+qq;
            meta.SeriesInstanceUID = uid;
            meta.SliceThickness = nii.hdr.dime.pixdim(4);
            meta.SpacingBetweenSlices = meta.SliceThickness;
            meta.PixelSpacing = nii.hdr.dime.pixdim(2:3);
            meta.ImageOrientationPatient = [1 0 0 0 1 0];
            meta.PatientName.FamilyName = sprintf('Case%03i', Order(ind));
            meta.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
            meta.StudyDescription = [examID '_ax'];
            meta.SeriesDescription = sprintf('%i', Order(ind));
            
            z_base = nii.hdr.hist.qoffset_z;
            
            fprintf('%s \n',QSMfile_list{qq});
            for ss = 1:size(img,3)
                meta.InstanceNumber = ss;
                meta.SliceLocation = z_base + (ss-1)*meta.SliceThickness;
                dcm_filename = sprintf('%s/%s/IM%04i.dcm',DCMout_root,[examID '_ax'],ss);
                dicomwrite(img(:,:,ss), dcm_filename, meta,...
                    'CreateMode', 'copy');
            end
            
        end
        
    end
    
end

%% streaking artifact cases

examID = 'Streaking_score2';

if exist(sprintf('%s/%s/IM%04i.dcm',DCMout_root,examID,1),'file') == 2
    % continue;
end

% load nifti file
% Streaking 0
% QSMpath = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/HDQSM/QSM_HDQSM_l10E-4p20.nii.gz';
% Streaking 1
% QSMpath = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/HDQSM/QSM_HDQSM_l10E-4p80.nii.gz';
% Streaking 2
% QSMpath = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/HDQSM/QSM_HDQSM_l10E-5p50.nii.gz';
% Streaking 3
% QSMpath = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/MEDI/QSM_MEDI_l10E5p00.nii.gz';
% OverRegularization
% QSMpath = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/HDQSM/QSM_HDQSM_l10E-3p50.nii.gz';
% signal variation
% QSMpath = '/data/7T_hunt/temp_062921_keep/for_archiving/swan_qsm/HDBET_allQSM/QSM_SSTGV_meanEcho.nii.gz';
% noise
% QSMpath = '/data/7T_hunt/b4480/t12992/swan_qsm/HDBET_allQSM/QSM_FANSI_nonlinearTV_meanEcho.nii.gz';
nii = load_untouch_nii(QSMpath);
img = double(nii.img);

img = img + 0.15;
img(img < 0) = 0; img(img > 0.3) = 0.3;
img = flip(flip(rot90(int16(img/0.3*4096))),3);

mkdir([DCMout_root '/' examID]);

uid = dicomuid;
meta.Modality = 'MR';
meta.SeriesNumber = 2001;
meta.SeriesInstanceUID = uid;
meta.SliceThickness = nii.hdr.dime.pixdim(4);
meta.SpacingBetweenSlices = meta.SliceThickness;
meta.PixelSpacing = nii.hdr.dime.pixdim(2:3);
meta.ImageOrientationPatient = [1 0 0 0 1 0];
meta.PatientName.FamilyName = examID;
meta.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
meta.StudyDescription = examID;
meta.SeriesDescription = sprintf('%i', Order(ind));

z_base = nii.hdr.hist.qoffset_z;

fprintf('%s \n',QSMfile_list{qq});
for ss = 1:size(img,3)
    meta.InstanceNumber = ss;
    meta.SliceLocation = z_base + (ss-1)*meta.SliceThickness;
    dcm_filename = sprintf('%s/%s/IM%04i.dcm',DCMout_root,examID,ss);
    dicomwrite(img(:,:,ss), dcm_filename, meta,...
        'CreateMode', 'copy');
end