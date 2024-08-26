clear all
syms epsilon delta Gamma beta
Gamma = 6.74*pi*1E6/(2*pi); % Much closer to the experimental diagram using Gamma/(2pi)
beta = 0.0230;
hold on
% f = @(delta,epsilon) delta - sqrt((epsilon).^2-(Gamma).^2) + beta.*epsilon.^2./Gamma;
% fp = fimplicit(f,[-10 10 -50 50]./Gamma);
X = linspace(-6E6,6E6,12E3);
Y = sqrt((1-2.*X.*beta./Gamma-sqrt(1-4.*X.*beta./Gamma-4.*(beta).^2))./(2.*(beta./Gamma).^2));
X_1 = -X;
Y_1 = Y./(2.*Gamma);
area(X_1,Y_1);
hold off
xlabel('${\Delta}$~(Hz)','FontSize',18,'Interpreter','latex');
ylabel('${|\alpha(\Delta)|/\left(\kappa+\gamma\right)}$','FontSize',18,'Interpreter','latex');
% xlim([-5E6 5E6]);
ylim([0.48 0.52]);
% ylim([0 .75]);