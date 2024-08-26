syms alpha Delta Q kappa gamma
kappa = 2;
gamma = 2;
% alpha = sqrt((delta)^2+(kappa+gamma)^2/4);
X = linspace(-5,5,100);
Y = sqrt((X).^2+(kappa+gamma)^2/4);
subplot(2,1,1);
area(X,Y,6);
ylim([0 6]);
% set(gca,'XTick',[], 'YTick', []);
xlabel('${\Delta}$','FontSize',18,'Interpreter','latex')
ylabel('${|\alpha|}$','FontSize',18,'Interpreter','latex')
subplot(2,1,2);
f = @(x,y) y - sqrt((x+12*(0.008).*y.^2+60*(0e-4).*y.^4).^2+(kappa+gamma)^2/4);
fimplicit(f,[-8 6 0 6])
xlabel('${\Delta}$','FontSize',18,'Interpreter','latex')
ylabel('${|\alpha|}$','FontSize',18,'Interpreter','latex')