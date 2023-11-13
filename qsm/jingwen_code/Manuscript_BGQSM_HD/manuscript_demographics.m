clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

% dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_MNI_QSM.mat';
dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iQSM_0213_erode.mat';
img_root = '/working/lupolab/jingwen/001_QSM/temp';

%% read table

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20230213.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
T(strcmp(T.status_reclass,'MISSING'),:) = [];

%% Load data

load(dataPath);

% HD_BGanalysis([HD_BGanalysis.CAG] < 39) = [];

indHC = strcmp({HD_BGanalysis.group},'HC');
indPM = strcmp({HD_BGanalysis.group},'PM');
indEM = strcmp({HD_BGanalysis.group},'EM');
indMan = strcmp({HD_BGanalysis.group},'Manifest');

statusList = nan(1,length(HD_BGanalysis));

age = [HD_BGanalysis.age];
sex = {HD_BGanalysis.sex};
CAG = [HD_BGanalysis.CAG];
volList = zeros(3,length(HD_BGanalysis));

statusList(indHC) = 1;
statusList(indPM) = 2;
statusList(indEM) = 3;
statusList(indMan) = 4;
statusList(CAG < 36) = nan;

%% demographics

CAPS = age.*(CAG-33.66)/432.3326;
CAPS(isnan(statusList)) = nan;

CAG = [HD_BGanalysis.CAG];
TMS = [HD_BGanalysis.TMS];
DCL = [HD_BGanalysis.DCL];
TFC = [HD_BGanalysis.TFC];

DART = [HD_BGanalysis.DART];
DARTtiming = [HD_BGanalysis.DARTtiming];
Flanker = [HD_BGanalysis.Flanker];
Match = [HD_BGanalysis.Match];
SetShift = [HD_BGanalysis.SetShift];

ageOnset = zeros(1,length(HD_BGanalysis));
for ii = 1:length(HD_BGanalysis)
    if isnan(CAG(ii))
        ageOnset(ii) = nan;
    else
        ageOnset(ii) = medianYearOnset(CAG(ii));
    end
end

AOO = ageOnset;
YTO = AOO - age;

% YTO = exp(4.4196 - 2.8102*CAPS);

%% create table

indPM = statusList == 2;
indEM = statusList == 3;
indMan = statusList == 4;
indPMfar = indPM & YTO > 15;
indPMnear = indPM & YTO <= 15;

statusList4 = statusList;
statusList4(indHC) = 0;
statusList4(indPMfar) = 1;
statusList4(indPMnear) = 2;
statusList4(indEM | indMan) = 3;

statusList3 = statusList;
statusList3(indHC) = 0;
statusList3(indPM) = 1;
statusList3(indEM | indMan) = 2;

status = [indHC' ~indHC' indPMfar' indPMnear' [indEM | indMan]'];

strNumber = printDemoNum(status);
strSex = printDemoSex(status, sex);
strAge = printDemo(status, age);
strCAG = printDemo(status, CAG);
strCAPS = printDemo(status, CAPS);
strYTO = printDemo(status, YTO);
strTMS = printDemo(status, TMS);
strTFC = printDemo(status, TFC);
strDCL = printDemo(status, DCL);
strDART = printDemo(status, DART);
strDARTtiming = printDemo(status, DARTtiming);
strFlanker = printDemo(status, Flanker);
strMatch = printDemo(status, Match);
strSetShift = printDemo(status, SetShift);
T = table(strNumber, strSex, strAge, strCAG, strCAPS, strYTO, strTMS, strTFC, strDCL, ...
    strDART, strDARTtiming, strFlanker, strMatch, strSetShift);
T.Properties.RowNames = {'Healthy Controls','All HD','PM far from onset','PM near onset','Manifest HD'};
Tt = rows2vars(T);
Tt.Properties.RowNames = {'Number','Sex (male/female)','Age (year)','CAG Repeat','CAPS','Year To Onset (year)'...
    'Total Motor Score','Total Functional Capacity','Diagnostic Confidence Level',...
    'DART Score','DART timing','Flanker  Score','Match Score','Set Shift Score'};
Tt.OriginalVariableNames = [];
disp(Tt);

writetable(Tt,[img_root '/demographics.xlsx'],'WriteRowNames',true);

%% helper function

function [strMet] = printDemoNum(statusList)

strMet = cell(size(statusList,2),1);
for ii = 1:size(statusList,2)
    num = sum(statusList(:,ii));
    strMet{ii} = sprintf('%i', num);
end

end

function [strMet] = printDemoSex(statusList, sex)

FMet = zeros(size(statusList,2),1);
MMet = zeros(size(statusList,2),1);
strMet = cell(size(statusList,2),1);
for ii = 1:size(statusList,2)
    FMet(ii) = sum(contains(sex(statusList(:,ii)),'F'));
    MMet(ii) = sum(contains(sex(statusList(:,ii)),'M'));
    if ii > 1
        [~,p,~,~] = prop_test([FMet(1) FMet(ii)] , [FMet(1)+MMet(1) FMet(ii)+MMet(ii)], 1);
        if p < 0.05
            strAdd = '*';
        else
            strAdd = '';
        end
    else
        strAdd = '';
    end
    strMet{ii} = sprintf('%i/%i %s', MMet(ii), FMet(ii), strAdd);
end

end

function [strMet, meanMet, stdMet] = printDemo(statusList, measure)

meanMet = zeros(size(statusList,2),1);
stdMet = zeros(size(statusList,2),1);
strMet = cell(size(statusList,2),1);
for ii = 1:size(statusList,2)
    meanMet(ii) = nanmean(measure(statusList(:,ii)));
    stdMet(ii) = nanstd(measure(statusList(:,ii)));
    if ii > 1 && ~isnan(meanMet(1))
        p = ranksum(measure(statusList(:,1)), measure(statusList(:,ii)));
        if p < 0.0001
            strAdd = '****';
        elseif p < 0.001
            strAdd = '***';
        elseif p < 0.01
            strAdd = '**';
        elseif p < 0.05
            strAdd = '*';
        else
            strAdd = '';
        end
    else
        strAdd = '';
    end
    if isnan(meanMet(ii))
        strMet{ii} = 'N/A';
    else
        strMet{ii} = sprintf('%.1f±%.1f %s', meanMet(ii), stdMet(ii), strAdd);
    end
end

end

function [ageOnset] = medianYearOnset(CAG)

syms f(x)
f(x) = (1+exp(pi/sqrt(3)*(-21.54-exp(9.56-0.146*CAG)+x)./(sqrt(35.55+exp(17.72-0.327*CAG))))).^-1 - 0.5;

tmp = vpasolve(f);
ageOnset = double(tmp);

end