clear all
syms alpha delta Gamma mu
Gamma = (4.44 + 2.30)*1E6; % (kappa+gamma)/2
mu = 0.0115;
X = linspace(-6E6,6E6,1.2E4);

Y_1 = 1./(sqrt(2).*mu).*sqrt(1+2.*mu.*X./Gamma-sqrt(1+4.*mu.*(X./Gamma-mu)));
X_1 = X;
area(X_1,Y_1./2,0.75);
xlabel('${\Delta}$~(MHz)','FontSize',18,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',18,'Interpreter','latex');
xlim([-6E6 6E6]);
% ylim([0.48 0.52]);
ylim([0.4 0.7]);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1E6)

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