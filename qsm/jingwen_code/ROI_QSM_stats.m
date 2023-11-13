function [QSMdata, QSMstats] = ROI_QSM_stats(ROIdata, ROIname, QSMdata, QSMstats, nn)

    ROIdata(ROIdata == 0) = [];
    ROIdata(isnan(ROIdata)) = [];
    ROIdata(isinf(ROIdata)) = [];
    
    ROImean = mean(ROIdata);
    ROIstd = std(ROIdata);
    ROImedian = median(ROIdata);
    ROImad = mad(ROIdata,1);
    
    Troi = table(ROIname,{ROIdata});
    QSMdata(nn).QSMtable = [QSMdata(nn).QSMtable; Troi];
    
    Troi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
    QSMstats(nn).QSMtable = [QSMstats(nn).QSMtable; Troi];
    
end