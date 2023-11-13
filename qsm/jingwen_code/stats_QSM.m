function [QSMdata, QSMstats] = ...
    stats_QSM(QSMfile_root, QSMfile_list, ROIfile, CSFfile, Maskfile, ...
    ROIlist, flag_erode)

if nargin < 7
    flag_erode = 0;
end
se = strel('disk',1);

% load ROI files
nii = load_nii(ROIfile);
QSM_ROI = double(nii.img);

nii = load_nii(CSFfile);
QSM_LatVen = double(nii.img);

nii = load_nii(Maskfile);
brain_mask = double(nii.img);

QSM_ROI(QSM_LatVen > 0) = 0;
QSM_ROI(brain_mask == 0) = 0;

% loop through QSM files
QSMstats(length(QSMfile_list)) = struct;
QSMdata(length(QSMfile_list)) = struct;
for nn = 1:length(QSMfile_list)
    
    fprintf('# Extracting ROIs for map: %s \n', QSMfile_list{nn});
    nii = load_nii([QSMfile_root '/' QSMfile_list{nn} '.nii.gz']);
    QSMmap = double(nii.img);
    
    % imagesc(QSMmap(:,:,60),[-0.15 0.15]); colormap('gray'); pause;
    
    QSMstats(nn).QSMname = QSMfile_list{nn};
    QSMdata(nn).QSMname = QSMfile_list{nn};
    QSMstats(nn).QSMtable = table;
    QSMdata(nn).QSMtable = table;
    
    for rr = 1:size(ROIlist,1)
        ROIindex = ROIlist{rr,2};
        ROIname = ROIlist(rr,1);
        
        ROImask = QSM_ROI == ROIindex;
        if flag_erode; ROImask = imerode(ROImask,se); end
        ROIdata = QSMmap(ROImask);
        [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    end
    
    ROIname = {'Thalamus_L'};
    ROImask = (QSM_ROI == 205) | (QSM_ROI == 207) | (QSM_ROI == 209) | (QSM_ROI == 211) ...
        | (QSM_ROI == 213) | (QSM_ROI == 215);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = QSMmap(ROImask);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    
    ROIname = {'Thalamus_R'};
    ROImask = (QSM_ROI == 206) | (QSM_ROI == 208) | (QSM_ROI == 210) | (QSM_ROI == 212) ...
        | (QSM_ROI == 214) | (QSM_ROI == 216);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = QSMmap(ROImask);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    
    ROIname = {'SubstantiaNigra_L'};
    ROImask = (QSM_ROI == 197) | (QSM_ROI == 199);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = QSMmap(ROImask);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    
    ROIname = {'SubstantiaNigra_R'};
    ROImask = (QSM_ROI == 198) | (QSM_ROI == 200);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = QSMmap(ROImask);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    
    ROIname = {'LatVen'};
    ROIdata = QSMmap(QSM_LatVen > 0);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    
    ROIname = {'WholeBrain'};
    ROIdata = QSMmap(brain_mask > 0);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn);
    
end