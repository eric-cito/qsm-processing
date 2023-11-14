clear; clc;
warning('off');

%% add path

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');

%% read in list

DCMout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/DICOM_revision_120';
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
    'QSM_QSMnet_meanEcho' ...
    'QSM_xQSM2_meanEcho' ...
    'QSM_iQSM2_meanEcho'};

load([DCMout_root '/../VisualList.mat']);

%% load subject list and save data path

subjSel = unique(T.subj);
examSel = cell(size(subjSel));

SubjList = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList.xlsx');
subjList = [SubjList.b_num];
examList = [SubjList.t_num];

for ii = 1:length(subjSel)
    examSel{ii} = examList{strcmp(subjList, subjSel{ii})};
end

subjSel{10} = 'temp_110620_keep';

%% randomly assign the case number and record the correspondance

rng(777);

subj = reshape(repmat(subjSel,[1 12])',120,1);
QSMfile = repmat(QSMfile_list(1:12)',[10 1]);

% randomly order the 120 cases
rng(120);
Order = randperm(120)';

% add corresponding examlist
exam = cell(size(subj));
for ii = 1:length(subj)
    exam{ii} = examSel{strcmp(subjSel, subj{ii})};
end

T = table(subj,exam,QSMfile,Order);

save([DCMout_root '/../VisualList_revision_120.mat'], 'T');

%% Loop through subjects

for ii = 1:length(subj)
    
    examPath = ['/data/7T_hunt/' subj{ii} '/' exam{ii}];
    examID = [subj{ii} '_' num2str(Order(ii))];
    
    if exist([examPath '/swan_qsm/HDBET_allQSM/' QSMfile{ii} '.nii.gz'], 'file') == 0
        
        fprintf('Cannot find %s file %s with index %i \n', subj{ii}, QSMfile{ii}, Order(ii));
        
    elseif exist(sprintf('%s/%s/IM%04i.dcm',DCMout_root,[examID '_ax'],1), 'file') == 0 ...
            || strcmp(QSMfile{ii},'QSM_xQSM2_meanEcho') || strcmp(QSMfile{ii},'QSM_iQSM2_meanEcho')
        
        fprintf('Saving %s file %s with index %i \n', subj{ii}, QSMfile{ii}, Order(ii));

        % load nifti file
        QSMfile_root = [examPath '/swan_qsm/HDBET_allQSM/'];
        QSMpath = [QSMfile_root '/' QSMfile{ii} '.nii.gz'];
        nii = load_untouch_nii(QSMpath);
        img = double(nii.img);
        
        if strcmp(QSMfile{ii},'QSM_QSMGAN_meanEcho')
            img = img/0.5684;
        end
        
        img = img + 0.15;
        img(img < 0) = 0; img(img > 0.3) = 0.3;
        img = flip(flip(rot90(int16(img/0.3*4096))),3);
        
        mkdir([DCMout_root '/' examID '_ax']);
        
        uid = dicomuid;
        meta.Modality = 'MR';
        meta.SeriesNumber = ii;
        meta.SeriesInstanceUID = uid;
        meta.SliceThickness = nii.hdr.dime.pixdim(4);
        meta.SpacingBetweenSlices = meta.SliceThickness;
        meta.PixelSpacing = nii.hdr.dime.pixdim(2:3);
        meta.ImageOrientationPatient = [1 0 0 0 1 0];
        meta.PatientName.FamilyName = sprintf('Revision120_Case%03i', Order(ii));
        meta.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        meta.StudyDescription = [examID '_ax'];
        meta.SeriesDescription = sprintf('%i', Order(ii));
        
        z_base = nii.hdr.hist.qoffset_z;
        for ss = 1:size(img,3)
            meta.InstanceNumber = ss;
            meta.SliceLocation = z_base + (ss-1)*meta.SliceThickness;
            dcm_filename = sprintf('%s/%s/IM%04i.dcm',DCMout_root,[examID '_ax'],ss);
            dicomwrite(img(:,:,ss), dcm_filename, meta,...
                'CreateMode', 'copy');
        end
        
    end
    
end
