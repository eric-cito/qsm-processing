function [ opt23, alphaopt ] = draw_freque( alpha, e1, e2, e3, fig )
% Draw the Mean mask amplitudes and Zeta funtions for the 
% susceptibility reconstruction frequency analysis.
% 
% This function outputs minimum zeta value between masks 2 and 3.
%
% Parameters:
% alpha: vector with regularization weight values
% e1, e2, e3: vectors with the mean amplitude values for each mask.
% fig: figure number to display the amplitudes and zeta functions 
%      (normalized squared difference).
%
% Last modified by Carlos Milovic, 08.07.2020 



zeta12 = zetafunc(e1,e2);
zeta13 = zetafunc(e1,e3);
zeta23 = zetafunc(e2,e3);

opt12  = find(zeta12 == min(zeta12));
opt13  = find(zeta13 == min(zeta13));
opt23  = find(zeta23 == min(zeta23));


% First plot the Mean mask amplitudes
figure(fig), subplot(1,2,1);
plot(alpha,e1,'-*r','LineWidth',1.0, 'MarkerSize',6);
hold on;
plot(alpha, e2,'-*g','LineWidth',1.0, 'MarkerSize',6);
plot(alpha, e3,'-*b','LineWidth',1.0, 'MarkerSize',6);
plot(alpha(opt12), e2(opt12),'d','Color',[0.5 0.5 0],'LineWidth',3.0, 'MarkerSize',15);
plot(alpha(opt13), e1(opt13),'d','Color',[0.5 0.0 0.5],'LineWidth',3.0, 'MarkerSize',15);
plot(alpha(opt23), e3(opt23),'d','Color',[0.0 0.5 0.5],'LineWidth',3.0, 'MarkerSize',15);
hold off;
set(gca,'FontSize',24)
set(gcf,'Color','white')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('Regularization weight')
ylabel('Mask Amplitude')
legend('A1','A2','A3','\zeta_{12}','\zeta_{13}','\zeta_{23}')
legend('Location','southwest')

% Now plot the zeta functions
figure(fig), subplot(1,2,2);
plot(alpha, zeta12,'-*','Color',[0.5 0.5 0],'LineWidth',1.0, 'MarkerSize',6);
hold on;
plot(alpha, zeta13,'-*','Color',[0.5 0.0 0.5],'LineWidth',1.0, 'MarkerSize',6);
plot(alpha, zeta23,'-*','Color',[0.0 0.5 0.5],'LineWidth',1.0, 'MarkerSize',6);
set(gca,'FontSize',24)
set(gcf,'Color','white')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%xticks('none')
xlabel('Regularization weight')
ylabel('\zeta function')
plot(alpha(opt12), zeta12(opt12),'d','Color',[0.5 0.5 0],'LineWidth',3.0, 'MarkerSize',15);
plot(alpha(opt13), zeta13(opt13),'d','Color',[0.5 0.0 0.5],'LineWidth',3.0, 'MarkerSize',15);
plot(alpha(opt23), zeta23(opt23),'d','Color',[0.0 0.5 0.5],'LineWidth',3.0, 'MarkerSize',15);
hold off;
legend('\zeta_{12}','\zeta_{13}','\zeta_{23}')
legend('Location','southeast')


disp(['Minimum Zeta23: ', num2str([alpha(opt23)])])
alphaopt = alpha(opt23);

end

