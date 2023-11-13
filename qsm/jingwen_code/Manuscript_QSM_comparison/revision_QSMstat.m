function [QSMstats] = ...
    revision_QSMstat(QSMmap, QSM_ROI, QSM_LatVen, brain_mask, ROIlist, flag_erode)

if nargin < 7
    flag_erode = 0;
end

offset = 0;

QSMstats = table;

for rr = 1:size(ROIlist,1)
    ROIindex = ROIlist{rr,2};
    ROIname = ROIlist(rr,1);
    
    ROImask = QSM_ROI == ROIindex;
    [ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
    
    Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
    QSMstats = [QSMstats; Troi];
end

ROIname = {'Thalamus_L'};
ROImask = (QSM_ROI == 205) | (QSM_ROI == 207) | (QSM_ROI == 209) | (QSM_ROI == 211) ...
    | (QSM_ROI == 213) | (QSM_ROI == 215);
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
QSMstats = [QSMstats; Troi];

ROIname = {'Thalamus_R'};
ROImask = (QSM_ROI == 206) | (QSM_ROI == 208) | (QSM_ROI == 210) | (QSM_ROI == 212) ...
    | (QSM_ROI == 214) | (QSM_ROI == 216);
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
QSMstats = [QSMstats; Troi];

ROIname = {'SubstantiaNigra_L'};
ROImask = (QSM_ROI == 197) | (QSM_ROI == 199);
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
QSMstats = [QSMstats; Troi];

ROIname = {'SubstantiaNigra_R'};
ROImask = (QSM_ROI == 198) | (QSM_ROI == 200);
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
QSMstats = [QSMstats; Troi];

ROIname = {'LatVen'};
ROImask = QSM_LatVen > 0;
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
QSMstats = [QSMstats; Troi];

ROIname = {'WholeBrain'};
ROImask = brain_mask > 0;
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
QSMstats = [QSMstats; Troi];

end

function [ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flagErode)

se = strel('sphere',1);

if flagErode
    ROImask = imerode(ROImask,se);
end
ROIdata = QSMmap(ROImask);

ROIdata(ROIdata == 0) = [];
ROIdata(isnan(ROIdata)) = [];
ROIdata(isinf(ROIdata)) = [];

ROImean = mean(ROIdata-offset);
ROIstd = std(ROIdata-offset);
ROImedian = median(ROIdata-offset);
ROImad = mad(ROIdata-offset,1);

end