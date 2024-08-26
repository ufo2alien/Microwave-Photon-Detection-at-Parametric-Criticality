clear all
% syms alpha delta Gamma mu
Gamma = (4.44 + 2.30)*1E6*pi; % (kappa+gamma)/2
omega0=6.1669e9;
gamma0=8.7e-3;
F=0.362*pi;
% B =(B0/gamma0)*cos(F)^3/sin(F)^2;
B =(Gamma/omega0/gamma0)/tan(F)^2;
X = linspace(-6E6,6E6,201)*2*pi;

Y_1 = 1./(sqrt(2).*B).*sqrt(1+2.*B.*X./Gamma-sqrt(1-4.*B.*(-X./Gamma+B)));
% Y_1 = sqrt((X./Gamma).^2+1);
%1./(sqrt(2).*B).*sqrt(1-2.*B.*-X./Gamma-sqrt(1-4.*B.*(-X./Gamma+B)));
X_1 = X;
figure (1);
area(X_1/2/pi,Y_1/2,.75,'FaceColor','#EDB120');
% colororder([0 0.5 1; 0.5 0 1]);
xlabel('${\Delta}$~(MHz)','FontSize',18,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',18,'Interpreter','latex');
 xlim([-5E6 5E6]);
% ylim([0.48 0.52]);
ylim([0.3 0.75]);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1E6)
%%
figure(2);
X_2 = X;
Y_2 = sqrt((X./Gamma).^2+1);
area(X_2,Y_2./2,0.75);
xlabel('${\Delta}$~(MHz)','FontSize',18,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',18,'Interpreter','latex');
xlim([-6E6 6E6]);
ylim([0.4 0.7]);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1E6)