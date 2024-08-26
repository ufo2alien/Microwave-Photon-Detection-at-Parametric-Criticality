close all
clear all

addpath('C:\Users\wangj21\Dropbox (Aalto)\SPDetectorJPA\PswFreqTemp');
syms alpha Delta kappa gamma Gamma K U alpha_c mu
hbar = 1.054571817E-34; % Reduced Planck constant
k_B = 1.380649E-23; % Boltzmann constant
Temp = [0.105 0.24 0.33 0.4 0.48 0.55 0.62]; % Temperature
omega_0 = 2*pi.*6.042E9;
kappa = 4.44E6*2*pi;
gamma = 2.30E6*2*pi;
% Gamma = 2*pi*5.6E6; 2*pi*4.06E6; (kappa+gamma)/2; 2*pi*4.06E6; (kappa+gamma)/2;
Gamma = (kappa+gamma)/2;
alpha =  3.75E6*2*pi; 3.34E6*2*pi; 3.75E6*2*pi; 4.20E6*2*pi; 3.00E6*2*pi; 3.41E6*2*pi; 3.44E6*2*pi;
% alpha = 1.02*Gamma;
Delta =  +2.75E6*2*pi; +2.00E6*2*pi; +2.75E6*2*pi; -0.00E6*2*pi; -2.75E6*2*pi; +0.75E6*2*pi; +0.95E6*2*pi;
K = -0.21E3*2*pi; -0.23E3*2*pi; -0.245E3*2*pi;
mu = 6.16E-3;
I = linspace(-70, 70, 4001); % Change this in order to show a full potential
U = meta_potential(I, Delta, alpha, K, Gamma);
detuning = 2*pi.*linspace(-6E6,6E6,1.201E3);
detuning_dark = 2*pi.*linspace(0.0E6,3.5E6,3.501E2);
% alpha_c = threshold(detuning, K, Gamma);
alpha_c = shifted(detuning, K, Gamma, mu);
equaloccupation = equalprob(detuning(detuning>0), K, Gamma);
tau_p = 20E-6;
rate_const = 2*pi*5.7E6;

%%
figure();
% subplot(1,2,1);
plot(I,U./(2.*Gamma),'LineWidth',2);
% set(gca,'XTick',[], 'YTick', []);
xticks([-60 -40 -20 0 20 40 60]);
% yticks([-80 -60 -40 -20 0 20 40 60 80]);

% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
% yticks([-6 -4 -2 0 2]);
% xticks([0]);
% yticks([0]);
ax = get(gca,'TickLabel');
set(gca,'TickLabel',ax,'FontSize',20);

% hh = gca;
% xRange = hh.XLim(2) - hh.XLim(1);
% xFract = (0 - hh.XLim(1))/xRange;
% xO = hh.Position(1) + xFract*hh.Position(3);
% yRange = hh.YLim(2) - hh.YLim(1);
% yFract = (0 - hh.YLim(1))/yRange;
% yO = hh.Position(2) + yFract*hh.Position(4);
% 
% % Remove one of the labels at the orogin
% ax = gca;
% % get the position of the zero tick on y-axis:
% zeroYtick = ax.YAxis.TickValues==0;
% % remove it (tick and lable)
% ax.YAxis.TickValues(zeroYtick) = [];
% % get the position of the zero tick on x-axis:
% zeroXtick = ax.XAxis.TickValues==0;
% % remove it (tick and lable)
% ax.XAxis.TickValues(zeroXtick) = [];
% % place a new zero at the origin:
% dim = [xO-0.2*xFract^2 yO-0.2*yFract^2 0.1*xFract 0.1*yFract];
% annotation('textbox',dim,'String','0','VerticalAlignment','top',...
%     'FitBoxToText','on','LineStyle','none','FontSize',24);
% % End

% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
xlabel('${\mathcal{Q}}$','FontSize',24,'Interpreter','latex');
ylabel('${\mathcal{U}(\mathcal{Q})/(\kappa+\gamma)}$','FontSize',24,'Interpreter','latex');
% ylabel('${\mathcal{U}(\mathcal{Q})}$','FontSize',32,'Interpreter','latex');
% xlabel('${\mathcal{Q}}$','FontSize',32,'Interpreter','latex');
% ylabel('${V(\mathcal{Q})}$','FontSize',18,'Interpreter','latex');
% xlabel('${I}$','FontSize',18,'Interpreter','latex');
% ylabel('${U(I)}$','FontSize',18,'Interpreter','latex');
% xlim([-35 35]);
% yticks([-2 -1 0 1]);
% set(gca,'Yticklabel',[]) 
% set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
Q_max = Qmax(Delta, K, Gamma, alpha);
U_max = meta_potential(Q_max, Delta, alpha, K, Gamma);
% U_max = barrier(Delta, alpha, K, Gamma);
U2d_0 = second_derivative_U_0(Delta, Gamma, alpha);
U2d_max = second_derivative_U_max(Delta, Gamma, alpha);
% title({['|Q_{max}|: ' num2str(Q_max) ', U(Q_{max})/(\kappa+\gamma): '...
%     num2str(U_max./(2.*Gamma)) ','] ['U^{(2)}(0)/(\kappa+\gamma): '...
%     num2str(U2d_0./(2.*Gamma)) ', U^{(2)}(Q_{max})/(\kappa+\gamma): '...
%     num2str(U2d_max./(2.*Gamma))]});
axis tight
% box off

%%
figure();
% subplot(1,2,2);
hold on
area([-6E6 6E6],[0 0],0.75,"FaceColor",'#2c7fb8');
area([0 6E6],[1/2 1/2],0.75,'FaceColor','#7fcdbb','EdgeColor','#0072BD');
area(detuning./(2*pi),alpha_c./(2.*Gamma),0.75,'FaceColor','#edf8b1','EdgeColor','#edf8b1');
plot(detuning(detuning>0)./(2*pi),equaloccupation./(2.*Gamma),'--','LineWidth',1.5);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1E6,'FontSize',14);
xlabel('${\Delta}/(2\pi)$~(MHz)','FontSize',18,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',18,'Interpreter','latex');
xlim([-6E6 6E6]);
ylim([0.30 0.75]);
yticks([0.3 0.4 0.5 0.6 0.7]);


x_pos = Delta./(2*pi);
y_pos = alpha./(2*Gamma);

plot(x_pos,y_pos,'*','Color','#dd1c77','LineWidth',1.5);
% annotation('textbox',[0.27 0.30 0.2 0.25],'String','I',...
%     'FontSize',24,'FitBoxToText','on','EdgeColor','none',...
%     'HorizontalAlignment','center','VerticalAlignment','baseline');
% annotation('textbox',[0.418 0.72 0.2 0.25],'String','II',...
%     'FontSize',24,'FitBoxToText','on','EdgeColor','none',...
%     'HorizontalAlignment','center','VerticalAlignment','baseline');
% annotation('textbox',[0.67 0.57 0.2 0.25],'String','III',...
%     'FontSize',24,'FitBoxToText','on','EdgeColor','none',...
%     'HorizontalAlignment','center','VerticalAlignment','baseline');

plot(-2.75E6,0.4451,'o','Color','#dd1c77','LineWidth',1.5); % #1
plot(+0.75E6,0.5059,'+','Color','#dd1c77','LineWidth',1.5); % #4
plot(+2.75E6,0.5564,'x','Color','#dd1c77','LineWidth',1.5); % #5
plot(+0.00E6,0.6231,'diamond','Color','#dd1c77','LineWidth',1.5); % #3
plot(+2.00E6,0.4955,'square','Color','#dd1c77','LineWidth',1.5); % #2
plot(-4.00E6,0.6000,'^','Color','#dd1c77','LineWidth',1.5); % #6
hold off


n_T = thermalphoton(omega_0,Temp,hbar,k_B);
% for ind_T=1:length(n_T)
% Gamma_dark(:,ind_T) = dark_switch(detuning_dark,K,Gamma,alpha,n_T(ind_T));
% Omega_sw(:,ind_T) = prefactor(detuning_dark, K, Gamma, alpha, n_T(ind_T));
% Exp(:,ind_T) = exponent(detuning_dark, K, Gamma, alpha, n_T(ind_T));
% p_dark(:,ind_T) = 1-exp(-1.*Exp(:,ind_T).*tau_p.*rate_const);
% end
% 
% figure();
% hold on
% for ind_T=1:length(n_T)
% plot(detuning_dark./(2*pi),Exp(:,ind_T),'LineWidth',1.0);
% end
% xlim([1.0E6 3.4E6]);
% xlabel('${\Delta/(2\pi)}$~(Hz)','FontSize',18,'Interpreter','latex');
% ylabel('${\Gamma_\mathrm{dark}}$','FontSize',18,'Interpreter','latex');
% legend({'0.105 K','0.24 K','0.33 K','0.4 K','0.48 K','0.55 K','0.62 K'},...
%     'Location','best');
% hold off


% openfig('Psw_vs_freq_detuning_no20mKUPD31012023.fig','new');
% fig = get(gca,'Children');
% xdata = get(fig, 'XData');
% ydata = get(fig, 'YData');
% % figure();
% hold on
% for ind_T=1:length(n_T)
% plot(detuning_dark./(2*pi),p_dark(:,ind_T),'o','LineWidth',0.5);
% end
% hold off
% % ax = gca;
% % ax.XLim = [0.0E6 3.5E6];
% xlim([1.0E6 3.5E6]);
% set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
% % set(ax,'xtick',1E6.*[1.0 2.0 3.0]);
% % xticks(1E6.*[0 0.5 1.0 1.5 2 2.5 3 3.5]);
% % xlim auto


n_T_20mK = thermalphoton(omega_0,20E-3,hbar,k_B);
n_T_200mK = thermalphoton(omega_0,200E-3,hbar,k_B);
n_T_500mK = thermalphoton(omega_0,500E-3,hbar,k_B);
% [X,Y] = meshgrid(detuning(1):2*pi*0.01E6:detuning(end),...
%     Gamma.*(0.3:1E-3:0.75));
[X,Y] = meshgrid(-2*pi*0.0E6:2*pi*1E4:2*pi*2.0E6,...
    2*Gamma.*(0.500:2E-4:0.550));
Z = dark_switch(X, K, Gamma, Y, n_T_20mK);
Z(imag(Z)~=0) = NaN;

Z_fixexpt = dark_switch_fixexpt(X, K, Gamma, Y, n_T_20mK);
Z_fixexpt(imag(Z_fixexpt)~=0) = NaN;

% mask = isfinite(Z);
% Z_sta = standardizeMissing(Z,NaN);
% Z_sta = standardizeMissing(Z,inf);
% Z_filled = fillmissing(Z,'constant',0);
% TF = ismissing(Z,inf);
% Z_masked = Z(mask);
boundary = threshold(X, K, Gamma);
% boundary = shifted(X, K, Gamma,mu);


%%
figure();
hold on
% ylim([0.50 0.550]);
% plot(X./(2*pi),boundary,'Color','r');
% surf(X./(2*pi),Y./(2*Gamma),Z_sta,'EdgeColor','none');
sg = surf(X./(2*pi)./(1E6),Y./(2*Gamma),Z,'EdgeColor','none');
surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Z)).*ones(size(Z)),...
    'EdgeColor','red','LineStyle','--','LineWidth',2);
% plot3(X./(2*pi)./(1E6),boundary./(2*Gamma),Z,...
%     'LineWidth',1.5,'Color','r','LineStyle','--');
view(2);
hold off
ylim([0.50 0.550]);
% zlim([min(min(Z)) max(max(Z))]);
colorbar;
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
title(['${\Gamma_\mathrm{dark}}$' ' @ 20 mK'],...
    'FontSize',18,'Interpreter','latex');
axis tight
aa = ancestor(sg,'axes');
set(aa,'yLim',[0.50 0.55]);
% 
% figure();
% surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Z)).*ones(size(Z)),...
%     'EdgeColor','red','LineStyle','--','LineWidth',2);
% view(2);
% axis tight
% % plot((-2*pi*0.0E6:2*pi*1E4:2*pi*2.0E6)./(2*pi)./(1E6),...
% %     shifted(-2*pi*0.0E6:2*pi*1E4:2*pi*2.0E6, K, Gamma,mu)./(2*Gamma));
% ylim([0.50 0.550]);



%%
Q_max = Qmax(X, K, Gamma, Y);
Q_max(imag(Q_max)~=0) = NaN;
eta = efficiency(Gamma, kappa, Z, Q_max, n_T_20mK);
eta_invalid = zeros(size(eta));
eta_fixexpt = efficiency(Gamma, kappa, Z_fixexpt, Q_max, n_T_20mK);
Exponent1 = exponent(X, K, Gamma, Y, n_T_20mK);
Exponent1(imag(Exponent1)~=0) = NaN;
for ind_row = 1:size(Exponent1,1)
    for ind_col = 1:size(Exponent1,2)
        if Y(ind_row,ind_col)-boundary(ind_row,ind_col)>0
            Exponent1(ind_row,ind_col) = NaN;
        else
            Exponent1(ind_row,ind_col)=Exponent1(ind_row,ind_col);
        end
    end
end
for ind_row = 1:size(eta,1)
    for ind_col = 1:size(eta,2)
%         if eta(ind_row,ind_col)>1 || ( Y(ind_row,ind_col)-...
%                 boundary(ind_row,ind_col)>-2*pi*2.4E4...
%                 && Y(ind_row,ind_col)-boundary(ind_row,ind_col)<+2*pi*4.20E3 )
        if eta(ind_row,ind_col)>1 || ((boundary(ind_row,ind_col)-...
                Y(ind_row,ind_col)<2*pi*3.5E4)...
                && (Y(ind_row,ind_col)-boundary(ind_row,ind_col)<+2*pi*1.00E3))
            eta(ind_row,ind_col) = NaN;
            eta(ind_row,ind_col) = 1;
            eta_invalid(ind_row,ind_col) = 1.0;
        else
%             [row_b, col_b] = find(~isnan(eta(ind_row,ind_col))&&abs(eta(ind_row,ind_col)-1)<0.1);
            eta(ind_row,ind_col)=eta(ind_row,ind_col);
            eta_invalid(ind_row,ind_col) = NaN;
        end
    end
end
% redChannel = eta > 1;
% greenChannel = eta < 1;
% blueChannel = greenChannel;
% colors = double(cat(3, redChannel, greenChannel, blueChannel));
figure();
hold on
sq = surf(X./(2*pi)./(1E6),Y./(2*Gamma),eta,'EdgeColor','none');
% sq = contour3(X./(2*pi)./(1E6),Y./(2*Gamma),eta,'LevelList',0.5,...
%     'Color','r','LineWidth',3);
% CD = sq.CData;
surf(X./(2*pi)./(1E6),boundary./(2*Gamma),1.*ones(size(Z)),...
    'EdgeColor',"#D95319",'LineStyle','--','LineWidth',2);
contour3(X./(2*pi)./(1E6),Y./(2*Gamma),Exponent1+4,'LevelList',0.5,...
    'Color',[0.4 0.4 0.4],'LineWidth',1);
% s_1 = surf(X./(2*pi)./(1E6),Y./(2*Gamma),eta_invalid,...
%     'FaceColor',[0.9769 0.9839 0.0805],'EdgeColor','none');

% % For plotting sparser mesh
% %%Extract X,Y and Z data from surface plot
% X_1=s_1.XData;
% Y_1=s_1.YData;
% Z_1=s_1.ZData;
% %%Divide the lengths by the number of lines needed
% xlength = size(Z_1,2);
% ylength = size(Z_1,1);
% xnumlines = 50; % 50 lines
% ynumlines = 50; % 50 partitions
% xspacing = round(xlength/xnumlines);
% yspacing = round(ylength/ynumlines);
% %%Plot the mesh lines 
% % Plotting lines in the X-Z plane
% % hold on
% for i = 1:yspacing:ylength
%   mesh([X_1(i,:);X_1(i,:)], [Y_1(i,:);Y_1(i,:)], [Z_1(i,:);Z_1(i,:)]);
% end
% % Plotting lines in the Y-Z plane
% for i = 1:xspacing:xlength
%   mesh([X_1(:,i),X_1(:,i)], [Y_1(:,i),Y_1(:,i)], [Z_1(:,i),Z_1(:,i)]);
% end
% % hold off
% % End plotting sparser mesh

scatter3(0.7,0.51,1,200,'r+','LineWidth',2.0);
colorbar;
ax = get(gca,'TickLabel');
set(gca,'TickLabel',ax,'FontSize',18);
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',20,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',20,'Interpreter','latex');
% title({['${\eta}$ @ 20 mK' ', the gray line is near $\eta=0.5$.']},...
%     'FontSize',18,'Interpreter','latex');
hold off
view(2);
axis tight
aa = ancestor(sq,'axes');
set(aa,'yLim',[0.50 0.55]);
set(aa,'zLim',[0 1]);
yticks([0.5 0.51 0.52 0.53 0.54 0.55]);

% % With fixed exponent of Gamma_dark
% figure();
% surf(X./(2*pi)./(1E6),Y./(2*Gamma),eta_fixexpt,'EdgeColor','none');
% colorbar;
% view(2);
% axis tight
% aa = ancestor(sq,'axes');
% set(aa,'yLim',[0.50 0.55]);
% yticks([0.5 0.51 0.52 0.53 0.54 0.55]);
% xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
% ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
% title({['${\eta}$' ' @ 20 mK, for which']...
%     ['$-\frac{2U(\mathcal{Q}_\mathrm{max})}{(n_{T}+1/2)(\kappa+\gamma)}=-4$ '...
%     'in $\Gamma_\mathrm{dark}$ is used.']},...
%     'FontSize',18,'Interpreter','latex');

% figure();
% contour3(X./(2*pi)./(1E6),Y./(2*Gamma),eta_fixexpt,'LevelList',1,'Color','red');

%%
Omega_sw = prefactor(X, K, Gamma, Y, n_T_20mK);
Omega_sw(imag(Omega_sw)~=0) = NaN;
figure();
hold on
so = surf(X./(2*pi)./(1E6),Y./(2*Gamma),Omega_sw,'EdgeColor','none');
surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Omega_sw)).*...
    ones(size(Z)),'EdgeColor','red','LineStyle','--','LineWidth',2);
colorbar;
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
title(['Prefactor ${\Omega_\mathrm{sw}}$' ' @ 20 mK'],'FontSize',18,...
    'Interpreter','latex');
hold off
view(2);
axis tight
aa = ancestor(so,'axes');
set(aa,'yLim',[0.50 0.55]);


%%
Exponential = exponential(X, K, Gamma, Y, n_T_20mK);
Exponential(imag(Exponential)~=0) = NaN;
figure();
hold on
surf(X./(2*pi)./(1E6),Y./(2*Gamma),Exponential,'EdgeColor','none');
se = surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Exponential)).*...
    ones(size(Z)),'EdgeColor','red','LineStyle','--','LineWidth',2);
colorbar;
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
title(['Exponential part of ${\Gamma_\mathrm{dark}}$' ' @ 20 mK'],...
    'FontSize',18,'Interpreter','latex');
hold off
view(2);
axis tight
aa = ancestor(se,'axes');
set(aa,'yLim',[0.50 0.55]);


%%
Exponent = exponent(X, K, Gamma, Y, n_T_20mK);
Exponent(imag(Exponent)~=0) = NaN;
figure();
hold on
for ind_row = 1:size(Exponent,1)
    for ind_col = 1:size(Exponent,2)
        if Y(ind_row,ind_col)-boundary(ind_row,ind_col)>0
            Exponent(ind_row,ind_col) = NaN;
        else
            Exponent(ind_row,ind_col)=Exponent(ind_row,ind_col);
        end
    end
end
surf(X./(2*pi)./(1E6),Y./(2*Gamma),Exponent,'EdgeColor','none');
contour3(X./(2*pi)./(1E6),Y./(2*Gamma),Exponent,'LevelList',-3.5,...
    'Color',[0.4 0.4 0.4],'LineWidth',1);
sept = surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Exponent)).*...
    ones(size(Z)),'EdgeColor','red','LineStyle','--','LineWidth',2);
colorbar;
ax = get(gca,'TickLabel');
set(gca,'TickLabel',ax,'FontSize',18);
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',20,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',20,'Interpreter','latex');
% title(['Exponent of the exponential part of ${\Gamma_\mathrm{dark}}$' ' @ 20 mK'],...
%     'FontSize',18,'Interpreter','latex');
hold off
view(2);
axis tight
aa = ancestor(sept,'axes');
set(aa,'yLim',[0.50 0.55]);
yticks([0.5 0.51 0.52 0.53 0.54 0.55]);

%%
Trial = trial_constUpp0(X, K, Gamma, Y, n_T_20mK);
Trial(imag(Trial)~=0) = NaN;
figure();
hold on
surf(X./(2*pi)./(1E6),Y./(2*Gamma),Trial,'EdgeColor','none');
st = surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Trial)).*...
    ones(size(Z)),'EdgeColor','red','LineStyle','--','LineWidth',2);
colorbar;
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
title(['${\Gamma_\mathrm{dark}}$ with const ${U^{(2)}(Q_0)}$' ' @ 20 mK'],...
    'FontSize',18,'Interpreter','latex');
hold off
view(2);
axis tight
aa = ancestor(st,'axes');
set(aa,'yLim',[0.50 0.55]);


%%
Trial_sw = prefactor_trial(X, K, Gamma, Y, n_T_20mK);
Trial_sw(imag(Trial_sw)~=0) = NaN;
figure();
hold on
surf(X./(2*pi)./(1E6),Y./(2*Gamma),Trial_sw,'EdgeColor','none');
sp = surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Trial_sw)).*...
    ones(size(Z)),'EdgeColor','red','LineStyle','--','LineWidth',2);
colorbar;
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
title(['${\Omega_\mathrm{sw}}$ with const ${U^{(2)}(Q_0)}$' ' @ 20 mK'],...
    'FontSize',18,'Interpreter','latex');
hold off
view(2);
axis tight
aa = ancestor(sp,'axes');
set(aa,'yLim',[0.50 0.55]);


%%
Q_max = Qmax(X, K, Gamma, Y);
Q_max(imag(Q_max)~=0) = NaN;
Q_0 = zeros(size(Z));
figure();
hold on
si = surf(X./(2*pi)./(1E6),boundary./(2*Gamma),max(max(Q_max)).*...
    ones(size(Z)),'EdgeColor','red','LineStyle','--','LineWidth',2);
% if Y>shifted(X, K, Gamma, mu)
    surf(X./(2*pi)./(1E6),Y./(2*Gamma),abs(Q_max),'EdgeColor','none');
% else
%     surf(X./(2*pi)./(1E6),Y./(2*Gamma),Q_0.^2,'EdgeColor','none');
% end
% if Y<shifted(0, K, Gamma, mu) | (Y<shifted(X, K, Gamma, mu)&X<0)
%     surf(X./(2*pi)./(1E6),Y./(2*Gamma),Q_0.^2,'EdgeColor','none');
% elseif Y>shifted(X, K, Gamma, mu)
%     surf(X./(2*pi)./(1E6),Y./(2*Gamma),Q_min.^2,'EdgeColor','none');
% else
%     surf(X./(2*pi)./(1E6),Y./(2*Gamma),Q_0.^2,'EdgeColor','none');
% end
% surf(X./(2*pi)./(1E6),Y./(2*Gamma),Q_0.^2,'EdgeColor','none');
colorbar;
xlabel('${\Delta/(2\pi)}$ (MHz)','FontSize',16,'Interpreter','latex');
ylabel('${\alpha/(\kappa+\gamma)}$','FontSize',16,'Interpreter','latex');
title('$|{Q}_\mathrm{max}|$','FontSize',18,...
    'Interpreter','latex');
hold off
view(2);
axis tight
% yticks([0.3 0.4 0.5 0.6 0.7]);
aa = ancestor(si,'axes');
set(aa,'yLim',[0.5 0.55]);




%% Another phase choice for alpha
% function V=meta_potential(I, Delta, alpha, K, Gamma)
% V = 1./(imag(alpha)+Gamma).*(-1/2.*(abs(alpha).^2-Gamma.^2-Delta.^2).*I.^2+...
%     3.*K.*Delta.*I.^4+6.*K.^2.*I.^6);
% end
% 
% function pump_critical=threshold(Delta, K, Gamma)
% pump_critical = sqrt(Delta.^2+Gamma.^2);
% end
% 
% function pump_shifted=shifted(Delta, K, Gamma, mu) % With pump-induced shift into account
% pump_shifted = Gamma.*(1./(sqrt(2)*mu)).*...
%     sqrt(1+2*mu.*Delta./Gamma-sqrt(1+4*mu.*(Delta./Gamma-mu)));
% end
% 
% function Q_max_abs=quadrature(Delta, K, Gamma, alpha)
% Q_max_abs = 1/sqrt(6.*abs(K)).*sqrt(Delta-sqrt((abs(alpha)).^2-Gamma.^2));
% end
% 
% function U_pp_0=second_derivative_U_0(Delta, Gamma, alpha)
% U_pp_0 = 1./(imag(alpha)+Gamma).*(Gamma.^2+Delta.^2-(abs(alpha)).^2);
% end
% 
% function U_pp_max=second_derivative_U_max(Delta, Gamma, alpha)
% U_pp_max = 4./(imag(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
%     Delta.*sqrt((abs(alpha)).^2-Gamma.^2));
% end
% 
% function U_Q_max=barrier(Delta, K, Gamma, alpha)
% U_Q_max = 1./(imag(alpha)+Gamma).*1/(36.*abs(K)).*...
%     (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
%     (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2));
% end
% 
% function n_T=thermalphoton(omega_0,Temp,hbar,k_B)
% n_T = 1./(exp(hbar*omega_0./(k_B.*Temp))-1);
% end
% 
% function Gamma_dark=dark_switch(Delta, K, Gamma, alpha, n_T)
% Gamma_dark = 2.*sqrt(1./(imag(alpha)+Gamma).*...
%     (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
%     abs(4./(imag(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
%     Delta.*sqrt((abs(alpha)).^2-Gamma.^2))))./(2*pi).*...
%     exp(-2.*1./(imag(alpha)+Gamma).*1/(36.*abs(K)).*...
%     (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
%     (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2))./(2.*Gamma.*(n_T+0.5)));
% end









%% Phase choice of pure imaginary alpha
function V=meta_potential(I, Delta, alpha, K, Gamma)
V = 1./(abs(alpha)+Gamma).*(-1/2.*(abs(alpha).^2-Gamma.^2-Delta.^2).*I.^2+...
    3.*K.*Delta.*I.^4+6.*K.^2.*I.^6);
end

function pump_critical=threshold(Delta, K, Gamma)
pump_critical = sqrt(Delta.^2+Gamma.^2);
end

function pump_equal=equalprob(Delta, K, Gamma)
pump_equal = sqrt((Delta./2).^2+Gamma.^2);
end

function pump_shifted=shifted(Delta, K, Gamma, mu) % With pump-induced shift into account
pump_shifted = Gamma.*(1./(sqrt(2)*mu)).*...
    sqrt(1+2*mu.*Delta./Gamma-sqrt(1+4*mu.*(Delta./Gamma-mu)));
end

function Q_max_abs=Qmax(Delta, K, Gamma, alpha)
Q_max_abs = 1/sqrt(6.*abs(K)).*sqrt(Delta-sqrt((abs(alpha)).^2-Gamma.^2));
end

function Q_min_abs=Qmin(Delta, K, Gamma, alpha)
Q_min_abs = 1/sqrt(6.*abs(K)).*sqrt(Delta+sqrt((abs(alpha)).^2-Gamma.^2));
end

function U_pp_0=second_derivative_U_0(Delta, Gamma, alpha)
U_pp_0 = 1./(abs(alpha)+Gamma).*(Gamma.^2+Delta.^2-(abs(alpha)).^2);
end

function U_pp_max=second_derivative_U_max(Delta, Gamma, alpha)
U_pp_max = 4./(abs(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
    Delta.*sqrt((abs(alpha)).^2-Gamma.^2));
end

function U_Q_max=barrier(Delta, K, Gamma, alpha)
U_Q_max = 1./(abs(alpha)+Gamma).*1/(36.*abs(K)).*...
    (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
    (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2));
end

function n_T=thermalphoton(omega_0,Temp,hbar,k_B)
n_T = 1./(exp(hbar*omega_0./(k_B.*Temp))-1);
end

function Gamma_dark=dark_switch(Delta, K, Gamma, alpha, n_T)
Gamma_dark = 2.*sqrt(1./(abs(alpha)+Gamma).*...
    (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
    abs(4./(abs(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
    Delta.*sqrt((abs(alpha)).^2-Gamma.^2))))./(2*pi).*...
    exp(-2.*1./(abs(alpha)+Gamma).*1/(36.*abs(K)).*...
    (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
    (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2))./(2.*Gamma.*(n_T+0.5)));
end

function Gamma_dark_fixexpt=dark_switch_fixexpt(Delta, K, Gamma, alpha, n_T)
Gamma_dark_fixexpt = 2.*sqrt(1./(abs(alpha)+Gamma).*...
    (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
    abs(4./(abs(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
    Delta.*sqrt((abs(alpha)).^2-Gamma.^2))))./(2*pi).*...
    exp(-4.0);
end

function Omega_sw=prefactor(Delta, K, Gamma, alpha, n_T)
Omega_sw = 2.*sqrt(1./(abs(alpha)+Gamma).*...
    (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
    abs(4./(abs(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
    Delta.*sqrt((abs(alpha)).^2-Gamma.^2))))./(2*pi);
% Omega_sw = 2.*sqrt(1./(abs(alpha)+Gamma).*...
%     (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
%     abs(-4./(abs(alpha)+Gamma).*(sqrt((abs(alpha)).^2-Gamma.^2).*...
%     (Delta-sqrt((abs(alpha)).^2-Gamma.^2)))))./(2*pi);
end

function Gamma_expl=exponential(Delta, K, Gamma, alpha, n_T) % No prefactor
Gamma_expl = exp(-2.*1./(abs(alpha)+Gamma).*1/(36.*abs(K)).*...
    (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
    (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2))./(2.*Gamma.*(n_T+0.5)));
end

function Gamma_expt=exponent(Delta, K, Gamma, alpha, n_T)
Gamma_expt = -2.*1./(abs(alpha)+Gamma).*1/(36.*abs(K)).*...
    (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
    (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2))./(2.*Gamma.*(n_T+0.5));
end

function Exp=exponential_Dykman(Delta, K, Gamma, alpha, n_T) % Approximated exponential
Exp = exp(-1./(2.*n_T+1).*1./(12.*abs(K)./(Gamma)).*...
    sqrt(abs(alpha./Gamma).^2-1).*...
    (Delta./(Gamma)-sqrt(abs(alpha./Gamma).^2-1)).^2);
end

function Gamma_dark=dark_Dykman(Delta, K, Gamma, alpha, n_T)
Gamma_dark = 2.*sqrt(1./(abs(alpha)+Gamma).*...
    (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
    abs(-4./(abs(alpha)+Gamma).*(sqrt((abs(alpha)).^2-Gamma.^2).*...
    (Delta-sqrt((abs(alpha)).^2-Gamma.^2)))))./(2*pi).*...
    exp(-1./(2.*n_T+1).*1./(12.*abs(K)./(Gamma)).*...
    sqrt(abs(alpha./Gamma).^2-1).*...
    (Delta./(Gamma)-sqrt(abs(alpha./Gamma).^2-1)).^2);
end

% eta_th = 1
function eta=efficiency(Gamma, kappa, Gamma_dark, Q_max, n_T)
eta = 1.*4*kappa.*Gamma_dark./(2*Gamma)^2.*Q_max.^2./(n_T+1/2).^2;
end


function trial=trial_constUpp0(Delta, K, Gamma, alpha, n_T)
trial = 2.*sqrt(1E-3.*(2*Gamma).*...
    abs(4./(abs(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
    Delta.*sqrt((abs(alpha)).^2-Gamma.^2))))./(2*pi).*...
    exp(-2.*1./(abs(alpha)+Gamma).*1/(36.*abs(K)).*...
    (Delta-sqrt((abs(alpha)).^2-Gamma.^2)).^2.*...
    (Delta+2.*sqrt((abs(alpha)).^2-Gamma.^2))./(2.*Gamma.*(n_T+0.5)));
end


function trial_sw=prefactor_trial(Delta, K, Gamma, alpha, n_T)
trial_sw = 2.*sqrt(1E-3.*(2*Gamma).*...
    abs(4./(abs(alpha)+Gamma).*((abs(alpha)).^2-Gamma.^2-...
    Delta.*sqrt((abs(alpha)).^2-Gamma.^2))))./(2*pi);
% Omega_sw = 2.*sqrt(1./(abs(alpha)+Gamma).*...
%     (Gamma.^2+Delta.^2-(abs(alpha)).^2).*...
%     abs(-4./(abs(alpha)+Gamma).*(sqrt((abs(alpha)).^2-Gamma.^2).*...
%     (Delta-sqrt((abs(alpha)).^2-Gamma.^2)))))./(2*pi);
end