close all
clear all

hbar = 1.054571817E-34; % Reduced Planck constant
k_B = 1.380649E-23; % Boltzmann constant
Temp = [0.105 0.24 0.33 0.4 0.48 0.55 0.62]; % Temperature
omega_0 = 2*pi.*6.042E9;
kappa = 4.44E6*2*pi;
gamma = 2.30E6*2*pi;
% Gamma = 2*pi*5.6E6; 2*pi*4.06E6; (kappa+gamma)/2; 2*pi*4.06E6; (kappa+gamma)/2;
Gamma = (kappa+gamma)/2;
alpha =  3.00E6*2*pi; 0.60*(kappa+gamma); 3.41E6*2*pi; 3.34E6*2*pi; 3.75E6*2*pi; 4.20E6*2*pi; 3.00E6*2*pi; 3.41E6*2*pi; 3.44E6*2*pi;
% alpha = 1.02*Gamma;
Delta = -2.75E6*2*pi; +0.75E6*2*pi; +2.00E6*2*pi; +2.75E6*2*pi; -0.00E6*2*pi; -2.75E6*2*pi; +0.75E6*2*pi; +0.95E6*2*pi;
K = -0.21E3*2*pi; -0.23E3*2*pi; -0.245E3*2*pi;
mu = 5.2E-3; 0.085; 6.16E-3;
I = linspace(-25, 25, 4001); % Change this in order to show a full potential
detuning = 2*pi.*linspace(-6E6,6E6,1.201E3);
alpha_c = threshold(detuning, K, Gamma);
alpha_s = shifted(detuning, K, Gamma, mu);
equaloccupation = equalprob(detuning(detuning>0), K, Gamma);

%%
figure();
% subplot(1,2,2);
hold on
% area([-6E6 6E6],[0 0],0.75,"FaceColor",'#2c7fb8');
% area([0 6E6],[1/2 1/2],0.75,'FaceColor','#7fcdbb','EdgeColor','#0072BD');
% area(detuning./(2*pi),alpha_c./(2.*Gamma),0.75,'FaceColor','#edf8b1','EdgeColor','#edf8b1');
plot(detuning./(2*pi),alpha_c./(2.*Gamma),'LineWidth',2);
plot(detuning./(2*pi),alpha_s./(2.*Gamma),'--','LineWidth',2);
% plot(detuning(detuning>0)./(2*pi),equaloccupation./(2.*Gamma),'--','LineWidth',1.5);
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/1E6,'FontSize',16);
xlabel('${\Delta}/(2\pi)$~(MHz)','FontSize',18,'Interpreter','latex');
ylabel('${|\alpha|/(\kappa+\gamma)}$','FontSize',18,'Interpreter','latex');
xlim([-6E6 6E6]);
ylim([0.30 0.75]);
yticks([0.3 0.4 0.5 0.6 0.7]);
hold off
grid on
box on

x_pos = Delta./(2*pi);
y_pos = alpha./(2*Gamma);










%%
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