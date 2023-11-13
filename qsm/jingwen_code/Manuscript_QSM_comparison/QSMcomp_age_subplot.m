function [] = QSMcomp_age_subplot(Nrow, Ncol, ii, meanStat, ageList, sexList, ...
    ROIname, ROIind, REFind, yrange)

Mind = strcmp(sexList,'M');
Find = strcmp(sexList,'F');

data = mean(squeeze(meanStat(:,ROIind:ROIind+1,ii)),2) - squeeze(meanStat(:,REFind,ii));
plot(ageList(Mind), data(Mind), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
plot(ageList(Find), data(Find), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2); hold on;
% err = mean(squeeze(stdROI(:,ROIind:ROIind+1,ii)),2);
% errorbar(ageList(Mind), data(Mind), err(Mind), 'o'); hold on;
% errorbar(ageList(Find), data(Find), err(Find), 'o');

ylim(yrange); % [-0.05 0.2]
xlim([20 80]);

[F,~] = fit(ageList(:),data(:),'poly1');
h = plot(F);
h.Color = [0.5 0.5 0.5];
[R,P] = corrcoef(ageList(:),data(:));
if ii == 1
    if P(1,2) < 0.0001
        legend({'Male', 'Female', sprintf('r = %.3f\np < 0.0001', R(1,2))},...
            'box','off','location','best');
    elseif P(1,2) < 0.001
        legend({'Male', 'Female', sprintf('r = %.3f\np < 0.001', R(1,2))},...
            'box','off','location','best');
    else
        legend({'Male', 'Female', sprintf('r = %.3f\np = %.3f', R(1,2), P(1,2))},...
            'box','off','location','best');
    end
else
    if P(1,2) < 0.0001
        legend(h,{sprintf('r = %.3f\np < 0.0001', R(1,2))},...
            'box','off','location','best');
    elseif P(1,2) < 0.001
        legend(h,{sprintf('r = %.3f\np < 0.001', R(1,2))},...
            'box','off','location','best');
    else
        legend(h,{sprintf('r = %.3f\np = %.3f', R(1,2), P(1,2))},...
            'box','off','location','best');
    end
end

if mod(ii,Ncol) == 1
    ylabel([ROIname ' Susc. (ppm)']); % Putamen
else
    ylabel('');
end
if ii > Ncol*(Nrow-1)
    xlabel('Age (years)');
else
    xlabel('');
end

end