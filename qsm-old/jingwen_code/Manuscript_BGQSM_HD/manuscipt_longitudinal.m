clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

% dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_Longitudinal_1006.mat';
dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_Longitudinal_1228_iQSM_erode.mat';
img_root = '/working/lupolab/jingwen/001_QSM/temp';

%% Load data

load(dataPath);

for ii = 12 % 1:length(HDsubj_long)
    if ~isnan(HDsubj_long(ii).scanday(3))
        HDsubj_long(ii).scanday(2) = HDsubj_long(ii).scanday(3);
        HDsubj_long(ii).age(2) = HDsubj_long(ii).age(3);
        HDsubj_long(ii).imData.QSMmadTP2 = HDsubj_long(ii).imData.QSMmadTP3;
        HDsubj_long(ii).imData.QSMmedianTP2 = HDsubj_long(ii).imData.QSMmedianTP3;
        HDsubj_long(ii).imData.VolumeTP2 = HDsubj_long(ii).imData.VolumeTP3;
    end
end

age = cell2mat({HDsubj_long.age}');
sexlist = categorical([HDsubj_long.sex]',{'F','M'});
CAG = [HDsubj_long.cag]';
CAPS = age.*(CAG-33.66)/432.3326;

%% demographics

ageOnset = zeros(length(CAG),1);
for ii = 1:length(CAG)
    if isnan(CAG(ii))
        ageOnset(ii) = nan;
    else
        ageOnset(ii) = medianYearOnset(CAG(ii));
    end
end

AOO = ageOnset;
YTO = AOO - age;

statusList = [HDsubj_long.status]';
statusList(strcmp(statusList,'PM') & YTO(:,1) < 15) = {'PM near'};
statusList(strcmp(statusList,'PM') & YTO(:,1) > 15) = {'PM far'};
statusList(strcmp(statusList,'EM')) = {'Manifest'};

statusList = categorical(statusList,{'HC','PM near','PM far','Manifest'});

%% Imaging metrics

Nroi = size(HDsubj_long(1).imData,1)/2;

fieldname = 'VolumeTP1';
[volMat1, ~, ~] = extractMetric(HDsubj_long, fieldname);

fieldname = 'QSMmedianTP1';
[qsmMat1, ~, ~] = extractMetric(HDsubj_long, fieldname);

fieldname = 'VolumeTP2';
[volMat2, ~, ~] = extractMetric(HDsubj_long, fieldname);

fieldname = 'QSMmedianTP2';
[qsmMat2, ~, ~] = extractMetric(HDsubj_long, fieldname);

fieldname = 'VolumeTP3';
[volMat3, ~, ~] = extractMetric(HDsubj_long, fieldname);

fieldname = 'QSMmedianTP3';
[qsmMat3, ~, ~] = extractMetric(HDsubj_long, fieldname);

ROI_name = {'CN','PU','GPe','GPi','TH','RN','SN','STN','DN','ST'};

% volMatAll = cat(3, volMat1, volMat2, volMat3);
% qsmMatAll = cat(3, qsmMat1, qsmMat2, qsmMat3);

volMatAll = cat(3, volMat1, volMat2);
qsmMatAll = cat(3, qsmMat1, qsmMat2);

%% plot one subject

% subjInd = 9;
% volSubj = 0.5*(table2array(HDsubj_long(9).imData(1:2:end,2:3:end)) + ...
%     table2array(HDsubj_long(9).imData(2:2:end,2:3:end)));
% qsmSubj = 0.5*(table2array(HDsubj_long(9).imData(1:2:end,3:3:end)) + ...
%     table2array(HDsubj_long(9).imData(2:2:end,3:3:end)));
% 
% figure('position', [100 0 1200 600]);
% subplot(1,2,1);
% bar(volSubj);
% xticks(1:length(ROI_name)); xticklabels(ROI_name);
% subplot(1,2,2);
% bar(qsmSubj);
% xticks(1:length(ROI_name)); xticklabels(ROI_name);

%% plot

set(0,'DefaultAxesFontSize', 14);

for indROI = 1:9
    
    figure('position', [100 0 1200 600]);
    subplot(231);
    p = plotPaired1yr(qsmMatAll, statusList, indROI)
    ylabel('Susceptibility (ppm)'); box off
    subplot(232);
    plotChange1yr(qsmMatAll, statusList, indROI);
    ylabel('\Delta Susceptibility (ppm)'); box off; % ylim([-0.015 0.015]);
    
    subplot(234);
    p = plotPaired1yr(volMatAll, statusList, indROI)
    ylabel('Corrected Volume (mL)'); box off
    subplot(235);
    plotChange1yr(volMatAll, statusList, indROI);
    ylabel('\Delta Corrected Volume (mL)'); box off; % ylim([-0.35 0.1]);
    
    subplot(133);
    plot2Dpaired1yrGray(volMatAll, qsmMatAll, statusList, indROI);
    xlabel('Corrected Volume (mL)');
    ylabel('Susceptibility (ppm)');
    title(ROI_name{indROI});
    
    pause(1); export_fig([img_root '/longitudinal_' ROI_name{indROI}], '-png','-transparent'); close;
    
end

%% scan interval

scanday = cell2mat({HDsubj_long.scanday}');
fprintf('Scan interval mean %.2f SD %.2f \n', median(scanday(:,2)), mad(scanday(:,2),1));

%% helper function

function [] = plot2Dpaired1yrGray(volMatAll, qsmMatAll, statusList, indROI)

volAll = squeeze(volMatAll(:,indROI,:));
qsmAll = squeeze(qsmMatAll(:,indROI,:));

xrange = prctile(volAll(:),[0 100]);
yrange = prctile(qsmAll(:),[0 100]);

xrange = xrange + [-0.05 0.05]*(xrange(2) - xrange(1));
yrange = yrange + [-0.05 0.05]*(yrange(2) - yrange(1));

status = categorical({'HC','PM far','PM near','Manifest'});

for ii = 1:length(statusList)
    plot(volAll(ii,:),qsmAll(ii,:),'k-'); hold on;
    
    switch statusList(ii)
        case 'HC'
            s1 = scatter(volAll(ii,1),qsmAll(ii,1),50,'o',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0 0.4470 0.7410]); hold on;
            s2 = scatter(volAll(ii,2),qsmAll(ii,2),50,'o','filled',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0 0.4470 0.7410], ...
                'MarkerFaceColor',[0 0.4470 0.7410]);
        case 'PM far'
            s3 = scatter(volAll(ii,1),qsmAll(ii,1),50,'s',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0.8500 0.3250 0.0980]); hold on;
            scatter(volAll(ii,2),qsmAll(ii,2),50,'s','filled',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0.8500 0.3250 0.0980], ...
                'MarkerFaceColor',[0.8500 0.3250 0.0980]);
        case 'PM near'
            s4 = scatter(volAll(ii,1),qsmAll(ii,1),50,'d',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0.9290 0.6940 0.1250]); hold on;
            scatter(volAll(ii,2),qsmAll(ii,2),50,'d','filled',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0.9290 0.6940 0.1250],...
                'MarkerFaceColor',[0.9290 0.6940 0.1250]);
        case 'Manifest'
            s5 = scatter(volAll(ii,1),qsmAll(ii,1),50,'^',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0.4940 0.1840 0.5560]); hold on;
            scatter(volAll(ii,2),qsmAll(ii,2),50,'^','filled',...
                'MarkerFaceAlpha',0.75,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
                'MarkerFaceColor',[0.4940 0.1840 0.5560]);
    end
    
end

legend([s1 s3 s4 s5 s2],{'HC','PM far','PM near','Manifest','follow-up'},...
    'box','off','location','best');

ylim(yrange); xlim(xrange);

end

function [] = plot2Dpaired1yr(volMatAll, qsmMatAll, statusList, indROI)

volAll = squeeze(volMatAll(:,indROI,:));
qsmAll = squeeze(qsmMatAll(:,indROI,:));

xrange = prctile(volAll(:),[0 100]);
yrange = prctile(qsmAll(:),[0 100]);

xrange = xrange + [-0.05 0.05]*(xrange(2) - xrange(1));
yrange = yrange + [-0.05 0.05]*(yrange(2) - yrange(1));

status = categorical({'HC','PM far','PM near','Manifest'});

for ii = 1:length(statusList)
    l = plot(volAll(ii,:),qsmAll(ii,:),'-'); hold on;
    
    switch statusList(ii)
        case 'HC'
            l.Color = [0 0.4470 0.7410];
        case 'PM far'
            l.Color = [0.8500 0.3250 0.0980];
        case 'PM near'
            l.Color = [0.9290 0.6940 0.1250];
        case 'Manifest'
            l.Color = [0.4940 0.1840 0.5560];
    end
    
end

for ii = 1:length(unique(status))
    set(gca,'ColorOrderIndex',ii);
    s1 = scatter(volAll(statusList == status(ii),1),qsmAll(statusList == status(ii),1),50,'o','filled',...
        'MarkerFaceAlpha',0.75); hold on;
    set(gca,'ColorOrderIndex',ii);
    s2 = scatter(volAll(statusList == status(ii),2),qsmAll(statusList == status(ii),2),60,'s','filled',...
        'MarkerFaceAlpha',0.75); hold on;
end

legend([s1 s2],{'baseline','follow-up'},'box','off','location','best');

ylim(yrange); xlim(xrange);

end

% function [] = plot2Dpaired(volMatAll, qsmMatAll, statusList, indROI)
% 
% volAll = squeeze(volMatAll(:,indROI,:));
% qsmAll = squeeze(qsmMatAll(:,indROI,:));
% 
% status = categorical({'HC','PM far','PM near','Manifest'});
% 
% for ii = 1:length(statusList)
%     if isnan(volAll(ii,2))
%         l = plot(volAll(ii,[1 3]),qsmAll(ii,[1 3]),'-'); hold on;
%     else
%         l = plot(volAll(ii,:),qsmAll(ii,:),'-'); hold on;
%     end
%     
%     switch statusList(ii)
%         case 'HC'
%             l.Color = [0 0.4470 0.7410];
%         case 'PM far'
%             l.Color = [0.8500 0.3250 0.0980];
%         case 'PM near'
%             l.Color = [0.9290 0.6940 0.1250];
%         case 'Manifest'
%             l.Color = [0.4940 0.1840 0.5560];
%     end
%     
% end
% 
% for ii = 1:length(unique(status))
%     set(gca,'ColorOrderIndex',ii);
%     s1 = scatter(volAll(statusList == status(ii),1),qsmAll(statusList == status(ii),1),'o',...
%         'MarkerFaceColor','auto'); hold on;
%     set(gca,'ColorOrderIndex',ii);
%     s2 = scatter(volAll(statusList == status(ii),2),qsmAll(statusList == status(ii),2),'s',...
%         'MarkerFaceColor','auto'); hold on;
%     set(gca,'ColorOrderIndex',ii);
%     s3 = scatter(volAll(statusList == status(ii),3),qsmAll(statusList == status(ii),3),'d',...
%         'MarkerFaceColor','auto'); hold on;
% end
% 
% legend([s1 s2 s3],{'baseline','1 yr','2 yr'},'box','off','location','best');
% 
% end

function [] = plotChange1yr(volMatAll, statusList, indROI)

volAll = squeeze(volMatAll(:,indROI,:));
volChange = volAll(:,2) - volAll(:,1);

meanVolChange = [mean(volChange(statusList == 'HC',:)) ...
    mean(volChange(statusList == 'PM far',:)) ...
    mean(volChange(statusList == 'PM near',:)) ...
    mean(volChange(statusList == 'Manifest',:)) ...
    mean(volChange(statusList ~= 'HC',:))];
stdVolChange = [std(volChange(statusList == 'HC',:),1) ...
    std(volChange(statusList == 'PM far',:),1) ...
    std(volChange(statusList == 'PM near',:),1) ...
    std(volChange(statusList == 'Manifest',:),1) ...
    std(volChange(statusList ~= 'HC',:),1)];

e = errorbar(meanVolChange, stdVolChange, 'k'); hold on;
e.LineStyle = 'none';
b = bar(meanVolChange, 0.5); hold on;
b.FaceColor = 'flat';
b.CData(1,:) = [0 0.4470 0.7410];
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.9290 0.6940 0.1250];
b.CData(4,:) = [0.4940 0.1840 0.5560];
b.CData(5,:) = [0.4660 0.6740 0.1880];

ymax = max(meanVolChange + stdVolChange);
ymin = min(meanVolChange - stdVolChange);
yrange = ymax - ymin;
ylim([ymin-0.1*yrange ymax+0.1*yrange]);

xlim([0 6]);
xticks([1:5]);
xticklabels({'HC','PM far','PM near','Manifest','All HD'});

end

function [p] = plotPaired1yr(volMatAll, statusList, indROI)

volAll = squeeze(volMatAll(:,indROI,:));

for ii = 1:length(statusList)
    switch statusList(ii)
        case 'HC'
            l1 = plot([1:2],volAll(ii,:),'-'); hold on;
            l1.Color = [0 0.4470 0.7410];
            scatter([1:2],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0 0.4470 0.7410]);
        case 'PM far'
            l2 = plot([3:4],volAll(ii,:),'-'); hold on;
            l2.Color = [0.8500 0.3250 0.0980];
            scatter([3:4],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0.8500 0.3250 0.0980]);
            l5 = plot([9:10],volAll(ii,:),'-'); hold on;
            l5.Color = [0.4660 0.6740 0.1880];
            scatter([9:10],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0.4660 0.6740 0.1880]);
        case 'PM near'
            l3 = plot([5:6],volAll(ii,:),'-'); hold on;
            l3.Color = [0.9290 0.6940 0.1250];
            scatter([5:6],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0.9290 0.6940 0.1250]);
            l5 = plot([9:10],volAll(ii,:),'-'); hold on;
            l5.Color = [0.4660 0.6740 0.1880];
            scatter([9:10],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0.4660 0.6740 0.1880]);
        case 'Manifest'
            l4 = plot([7:8],volAll(ii,:),'-'); hold on;
            l4.Color = [0.4940 0.1840 0.5560];
            scatter([7:8],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0.4940 0.1840 0.5560]);
            l5 = plot([9:10],volAll(ii,:),'-'); hold on;
            l5.Color = [0.4660 0.6740 0.1880];
            scatter([9:10],volAll(ii,:),50,'o','filled','MarkerFaceAlpha',0.75, ...
                'MarkerFaceColor',[0.4660 0.6740 0.1880]);
    end
    
end

meanHC = nanmean(volAll(statusList == 'HC',:),1);
meanHC(sum(~isnan(volAll(statusList == 'HC',:))) == 1) = nan;
plot([1:2],meanHC,'-','Color',[0 0.4470 0.7410],'LineWidth',2);
p0 = signrank(volAll(statusList == 'HC',1), volAll(statusList == 'HC',2));
% [~,p1] = ttest(volAll(statusList == 'HC',1), volAll(statusList == 'HC',2));
sigstar({[1 2]},[p0]);

meanHC = nanmean(volAll(statusList == 'PM far',:),1);
meanHC(sum(~isnan(volAll(statusList == 'PM far',:))) == 1) = nan;
plot([3:4],meanHC,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
p1 = signrank(volAll(statusList == 'PM far',1), volAll(statusList == 'PM far',2));
% [~,p1] = ttest(volAll(statusList == 'PM far',1), volAll(statusList == 'PM far',2));
sigstar({[3 4]},[p1]);

meanHC = nanmean(volAll(statusList == 'PM near',:),1);
meanHC(sum(~isnan(volAll(statusList == 'PM near',:))) == 1) = nan;
plot([5:6],meanHC,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
p2 = signrank(volAll(statusList == 'PM near',1), volAll(statusList == 'PM near',2));
% [~,p1] = ttest(volAll(statusList == 'PM near',1), volAll(statusList == 'PM near',2));
sigstar({[5 6]},[p2]);

meanHC = nanmean(volAll(statusList == 'Manifest',:),1);
meanHC(sum(~isnan(volAll(statusList == 'Manifest',:))) == 1) = nan;
plot([7:8],meanHC,'-','Color',[0.4940 0.1840 0.5560],'LineWidth',2);
p3 = signrank(volAll(statusList == 'Manifest',1), volAll(statusList == 'Manifest',2));
% [~,p1] = ttest(volAll(statusList == 'Manifest',1), volAll(statusList == 'Manifest',2));
sigstar({[7 8]},[p3]);

meanHC = nanmean(volAll(statusList ~= 'HC',:),1);
meanHC(sum(~isnan(volAll(statusList ~= 'HC',:))) == 1) = nan;
plot([9:10],meanHC,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2);
p4 = signrank(volAll(statusList ~= 'HC',1), volAll(statusList ~= 'HC',2));
% [~,p1] = ttest(volAll(statusList ~= 'HC',1), volAll(statusList ~= 'HC',2));
sigstar({[9 10]},[p4]);

xlim([0 11]);
xticks([1.5:2:9.5]);
xticklabels({'HC','PM far','PM near','Manifest','All HD'});

p = [p0 p1 p2 p3 p4];

end

% function [] = plotPaired(volMatAll, statusList, indROI)
% 
% volAll = squeeze(volMatAll(:,indROI,:));
% 
% for ii = 1:length(statusList)
%     switch statusList(ii)
%         case 'HC'
%             l1 = plot([1:3],volAll(ii,:),'o-','MarkerFaceColor','auto'); hold on;
%             l1.Color = [0 0.4470 0.7410];
%         case 'PM far'
%             l2 = plot([4:6],volAll(ii,:),'o-','MarkerFaceColor','auto'); hold on;
%             l2.Color = [0.8500 0.3250 0.0980];
%         case 'PM near'
%             l3 = plot([7:9],volAll(ii,:),'o-','MarkerFaceColor','auto'); hold on;
%             l3.Color = [0.9290 0.6940 0.1250];
%         case 'Manifest'
%             if isnan(volAll(ii,2))
%                 l4 = plot([10 12],volAll(ii,[1 3]),'o-','MarkerFaceColor','auto'); hold on;
%             else
%                 l4 = plot([10:12],volAll(ii,:),'o-','MarkerFaceColor','auto'); hold on;
%             end
%             l4.Color = [0.4940 0.1840 0.5560];
%     end
%     
% end
% 
% meanHC = nanmean(volAll(statusList == 'HC',:),1);
% meanHC(sum(~isnan(volAll(statusList == 'HC',:))) == 1) = nan;
% plot([1:3],meanHC,'-','Color',[0 0.4470 0.7410],'LineWidth',2);
% [~,p1] = ttest(volAll(statusList == 'HC',1), volAll(statusList == 'HC',2));
% [~,p2] = ttest(volAll(statusList == 'HC',1), volAll(statusList == 'HC',3));
% sigstar({[1 2],[1 3]},[p1 p2]);
% 
% meanHC = nanmean(volAll(statusList == 'PM far',:),1);
% meanHC(sum(~isnan(volAll(statusList == 'PM far',:))) == 1) = nan;
% plot([4:6],meanHC,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
% [~,p1] = ttest(volAll(statusList == 'PM far',1), volAll(statusList == 'PM far',2));
% [~,p2] = ttest(volAll(statusList == 'PM far',1), volAll(statusList == 'PM far',3));
% sigstar({[4 5],[4 6]},[p1 p2]);
% 
% meanHC = nanmean(volAll(statusList == 'PM near',:),1);
% meanHC(sum(~isnan(volAll(statusList == 'PM near',:))) == 1) = nan;
% plot([7:9],meanHC,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',2);
% [~,p1] = ttest(volAll(statusList == 'PM near',1), volAll(statusList == 'PM near',2));
% [~,p2] = ttest(volAll(statusList == 'PM near',1), volAll(statusList == 'PM near',3));
% sigstar({[7 8],[7 9]},[p1 p2]);
% 
% meanHC = nanmean(volAll(statusList == 'Manifest',:),1);
% meanHC(sum(~isnan(volAll(statusList == 'Manifest',:))) == 1) = nan;
% plot([10:12],meanHC,'-','Color',[0.4940 0.1840 0.5560],'LineWidth',2);
% [~,p1] = ttest(volAll(statusList == 'Manifest',1), volAll(statusList == 'Manifest',2));
% [~,p2] = ttest(volAll(statusList == 'Manifest',1), volAll(statusList == 'Manifest',3));
% sigstar({[10 11],[10 12]},[p1 p2]);
% 
% xlim([0 13]);
% xticks([2:3:12]);
% xticklabels({'HC','PM far','PM near','Manifest'});
% 
% end
% 
% function [] = plot2YTO(volMatAll, statusList, YTO, indROI)
% 
% volAll = squeeze(volMatAll(:,indROI,:));
% 
% for ii = 1:length(statusList)
%     switch statusList(ii)
%         case 'HC'
%             plot(YTO(ii,:),volAll(ii,:),'ko-'); hold on;
%         case 'PM far'
%             plot(YTO(ii,:),volAll(ii,:),'bo-'); hold on;
%         case 'PM near'
%             plot(YTO(ii,:),volAll(ii,:),'co-'); hold on;
%         case 'EM'
%             if isnan(volAll(ii,2))
%                 plot(YTO(ii,[1 3]),volAll(ii,[1 3]),'mo-'); hold on;
%             else
%                 plot(YTO(ii,:),volAll(ii,:),'mo-'); hold on;
%             end
%         case 'Manifest'
%             plot(YTO(ii,:),volAll(ii,:),'ro-'); hold on;
%     end
%     
% end
% 
% end

function [ageOnset] = medianYearOnset(CAG)

syms f(x)
f(x) = (1+exp(pi/sqrt(3)*(-21.54-exp(9.56-0.146*CAG)+x)./(sqrt(35.55+exp(17.72-0.327*CAG))))).^-1 - 0.5;

tmp = vpasolve(f);
ageOnset = double(tmp);

end

function [matAll, matL, matR] = extractMetric(HD_BGanalysis, fieldname)

matL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
matR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
for ii = 1:length(HD_BGanalysis)
    matL(ii,:) = HD_BGanalysis(ii).imData.(fieldname)(1:2:end);
    matR(ii,:) = HD_BGanalysis(ii).imData.(fieldname)(2:2:end);
end
matAll = (matL + matR)/2;

end
