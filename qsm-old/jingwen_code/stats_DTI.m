function [QSMdata, QSMstats] = ...
    stats_DTI(FAfile, ROIfile, CSFfile, Maskfile, ROIlist, flag_erode, nn)

if nargin < 6
    flag_erode = 0;
end
se = strel('disk',1);

% load ROI files
nii = load_nii(ROIfile);
ROI = double(nii.img);

nii = load_nii(CSFfile);
QSM_LatVen = double(nii.img);

nii = load_nii(Maskfile);
brain_mask = double(nii.img);

ROI(QSM_LatVen > 0) = 0;
ROI(brain_mask == 0) = 0;

% loop through QSM files
QSMstats = struct;
QSMdata = struct;

% load FA
nii = load_nii(FAfile);
FAmap = double(nii.img);

QSMstats(1).QSMtable = table;
QSMdata(1).QSMtable = table;

for rr = 1:size(ROIlist,1)
    ROIindex = ROIlist{rr,2};
    ROIname = ROIlist(rr,1);
    
    ROImask = ROI == ROIindex;
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = FAmap(ROImask);
    [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);
end

ROIname = {'Thalamus_L'};
ROImask = (ROI == 205) | (ROI == 207) | (ROI == 209) | (ROI == 211) ...
    | (ROI == 213) | (ROI == 215);
if flag_erode; ROImask = imerode(ROImask,se); end
ROIdata = FAmap(ROImask);
[QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);

ROIname = {'Thalamus_R'};
ROImask = (ROI == 206) | (ROI == 208) | (ROI == 210) | (ROI == 212) ...
    | (ROI == 214) | (ROI == 216);
if flag_erode; ROImask = imerode(ROImask,se); end
ROIdata = FAmap(ROImask);
[QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);

ROIname = {'SubstantiaNigra_L'};
ROImask = (ROI == 197) | (ROI == 199);
if flag_erode; ROImask = imerode(ROImask,se); end
ROIdata = FAmap(ROImask);
[QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);

ROIname = {'SubstantiaNigra_R'};
ROImask = (ROI == 198) | (ROI == 200);
if flag_erode; ROImask = imerode(ROImask,se); end
ROIdata = FAmap(ROImask);
[QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);

ROIname = {'LatVen'};
ROIdata = FAmap(QSM_LatVen > 0);
[QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);

ROIname = {'WholeBrain'};
ROIdata = FAmap(brain_mask > 0);
[QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, 1);
