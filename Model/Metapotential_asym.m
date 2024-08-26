% SPD paper
% syms g I Q mu
% % alpha - pump strength, V - effective potential, Delta - detuning,
% % K - Kerr nonlinearity, L - sextic nonlinearity, eta - coherent state,
% % phi - phase angle of eta
close all
clear all
% syms alpha V eta Delta K theta L W mu...
%     eta_plus_squared eta_minus_squared V_depth alpha_0 alpha_t alpha_c
syms alpha V eta Delta K theta L mu...
    eta_plus_squared eta_minus_squared V_depth alpha_0 alpha_t alpha_c
kappa = 4.44E6*2*pi;
gamma = 2.30E6*2*pi;
% alpha = 2.0E6*2*pi; % 2 or 6 or 5 or 46
alpha = 2.0E6*2*pi;
Delta = 0.7E6*2*pi; % 6 or -2 or 8 or 50
K = -0.23E3*2*pi; % -0.01
L = 5.7E-3*2*pi; % 1E-5
x = linspace(-100, 100, 101); % -15 to 15
y = linspace(-100, 100, 101); % -15 to 15
[X, Y] = ndgrid(x, y);
eta = X + 1i*Y;
% phi = 0:pi/4:pi;
% b = 5*exp(phi.*1i).*1E3;
b = 40*exp(0*1i).*1E3;
% phi = angle(b);
theta = angle(eta);
V  = (abs(eta)).^2.*(Delta+alpha.*cos(2.*theta)+6*K.*(abs(eta)).^2+20*L.*(abs(eta)).^4); % Simplified meta-potential
alpha_0 = Delta-3*K^2/(5*L); %% Zero-amplitude pumping strength
alpha_t = Delta-9*K^2/(20*L); %% Transition pumping strength
alpha_c = Delta; %% Critical pumping strength

figure(1);
% % Surface plot of effective potential V
% sgtitle(['$\alpha$ = ',num2str(alpha),', ','$\Delta$ = ',num2str(Delta), ',','$V_{depth} = $',num2str(V_depth)], 'Interpreter','latex');
subplot(2,2,1);
surf(real(eta),imag(eta),V./(kappa+gamma));
axis tight;
% set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
zlabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% % Real axis
subplot(2,2,2);
% figure(2);
Y_idx = find(y >= 0, 1, 'first');
plot(X(:,Y_idx), V(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
% % Imaginary axis
subplot(2,2,3);
% figure(3);
X_idx = find(x >= 0, 1, 'first');
plot(Y(X_idx,:), V(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% set(gca,'XTick',[], 'YTick', []);
subplot(2,2,4);
% figure(4);
contour(real(eta),imag(eta),V./(kappa+gamma),50,'ShowText','on');
axis tight;
% set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
zlabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');

eta_plus_squared = 1/2*(-(K)/(5*L)+sqrt(K^2/(25*L^2)-1/(15*L)*(Delta-alpha))); %% Minima/Wells; location of minima, i.e., sqrt(eta_plus_squared), corresponds to coherent oscillation amplitude
eta_minus_squared = 1/2*(-(K)/(5*L)-sqrt(K^2/(25*L^2)-1/(15*L)*(Delta-alpha))); %% Maxima/Hills
V_depth = (Delta-alpha)*eta_minus_squared+6*K.*(eta_minus_squared).^2+20*L.*(eta_minus_squared).^3; %% Excitation energy
sgtitle(['$\alpha$ = ',num2str(alpha),', ', '$\alpha_0$ = ',num2str(alpha_0),', ', '$\alpha_t$ = ',num2str(alpha_t),', ','$\Delta$ = ',num2str(Delta), ', ','$V_{depth}$ = ',num2str(V_depth), ', ','Maxima at $\pm$',num2str(sqrt(eta_minus_squared)), ', ','Minima at $\pm$',num2str(sqrt(eta_plus_squared))], 'Interpreter','latex');
% figure(2);
% surf(real(eta),imag(eta),W);

mu = 0.024;
% W = (abs(eta)).^2.*(Delta-2.*mu./(kappa+gamma).*alpha.^2+alpha.*cos(2.*theta)+...
%     6*K.*(abs(eta)).^2+20*L.*(abs(eta)).^4)+2*abs(b).*abs(eta).*sqrt(kappa).*sin(phi-theta); % Full meta potential with probe field added

% figure(2);
% % % Surface plot of effective potential W
% % sgtitle(['$\alpha$ = ',num2str(alpha),', ','$\Delta$ = ',num2str(Delta), ',','$V_{depth} = $',num2str(V_depth)], 'Interpreter','latex');
% subplot(2,2,1);
% surf(real(eta),imag(eta),W./(kappa+gamma));
% axis tight;
% % set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
% xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
% ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
% zlabel('${W(\eta)}$','FontSize',18,'Interpreter','latex');
% % % Real axis
% subplot(2,2,2);
% % figure(2);
% Y_idx = find(y >= 0, 1, 'first');
% plot(X(:,Y_idx), W(:,Y_idx)./(kappa+gamma),'LineWidth',2);
% axis tight;
% xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
% ylabel('${W(\eta)}$','FontSize',18,'Interpreter','latex');
% set(gca,'XTick',[], 'YTick', []);
% % % Imaginary axis
% subplot(2,2,3);
% % figure(3);
% X_idx = find(x >= 0, 1, 'first');
% plot(Y(X_idx,:), W(X_idx,:)./(kappa+gamma),'LineWidth',2);
% axis tight;
% xlabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
% ylabel('${W(\eta)}$','FontSize',18,'Interpreter','latex');
% % set(gca,'XTick',[], 'YTick', []);
% subplot(2,2,4);
% % figure(4);
% contour(real(eta),imag(eta),W./(kappa+gamma),50,'ShowText','on');
% axis tight;
% % set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
% xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
% ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
% zlabel('${W(\eta)}$','FontSize',18,'Interpreter','latex');
% 
% figure(5);
% % % Real axis
% Y_idx = find(y >= 0, 1, 'first');
% plot(X(:,Y_idx), W(:,Y_idx)./(kappa+gamma),'LineWidth',2);
% axis tight;
% xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
% ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% set(gca,'XTick',[], 'YTick', []);

figure(2);
% % Surface plot of effective potential W
sgtitle(['Cross section of ${V(\eta)}$ along $\Im{(\eta)}=0$ plane,',' with $\tilde{b}=4\times10^4e^{i\varphi}$ Hz$^{1/2}$'],'FontSize',18,'Interpreter','latex');
subplot(3,3,1);
W1 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 0)./(kappa+gamma);
plot(X(:,Y_idx), W1(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=0}$','FontSize',18,'Interpreter','latex');
subplot(3,3,2);
W2 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, pi/4)./(kappa+gamma);
plot(X(:,Y_idx), W2(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,3);
W3 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, pi/2)./(kappa+gamma);
Y_idx = find(y >= 0, 1, 'first');
plot(X(:,Y_idx), W3(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=\pi/2}$','FontSize',18,'Interpreter','latex');
subplot(3,3,4);
W4 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 3*pi/4)./(kappa+gamma);
X_idx = find(x >= 0, 1, 'first');
plot(X(:,Y_idx), W4(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=3\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,5);
W5 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, pi)./(kappa+gamma);
plot(X(:,Y_idx), W5(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=\pi}$','FontSize',18,'Interpreter','latex');
subplot(3,3,6);
W6 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 5*pi/4)./(kappa+gamma);
plot(X(:,Y_idx), W6(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=5\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,7);
W7 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 3*pi/2)./(kappa+gamma);
plot(X(:,Y_idx), W7(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=3\pi/2}$','FontSize',18,'Interpreter','latex');
subplot(3,3,8);
W8 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 7*pi/4)./(kappa+gamma);
plot(X(:,Y_idx), W8(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=7\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,9);
W9 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 2*pi)./(kappa+gamma);
plot(X(:,Y_idx), W9(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=2\pi}$','FontSize',18,'Interpreter','latex');

figure(3);
sgtitle(['Cross section of ${V(\eta)}$ along $\Im{(\eta)}=0$ plane,',' with $\tilde{b}=4\times10^4e^{i\varphi}$ Hz$^{1/2}$'],'FontSize',18,'Interpreter','latex');
% % Surface plot of effective potential V
% sgtitle(['$\alpha$ = ',num2str(alpha),', ','$\Delta$ = ',num2str(Delta), ',','$V_{depth} = $',num2str(V_depth)], 'Interpreter','latex');
% W1 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 0)./(kappa+gamma);
subplot(2,2,1);
surf(real(eta),imag(eta),W1./(kappa+gamma));
axis tight;
% set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
zlabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% % Real axis
subplot(2,2,2);
% figure(2);
Y_idx = find(y >= 0, 1, 'first');
plot(X(:,Y_idx), W1(:,Y_idx)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
% % Imaginary axis
subplot(2,2,3);
% figure(3);
X_idx = find(x >= 0, 1, 'first');
plot(Y(X_idx,:), W1(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% set(gca,'XTick',[], 'YTick', []);
subplot(2,2,4);
% figure(4);
contour(real(eta),imag(eta),W1./(kappa+gamma),50,'ShowText','on');
axis tight;
% set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
zlabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');

figure(4);
% % Surface plot of effective potential W
sgtitle(['Cross section of ${V(\eta)}$ along $\Re{(\eta)}=0$ plane,',' with $\tilde{b}=4\times10^4e^{i\varphi}$ Hz$^{1/2}$'],'FontSize',18,'Interpreter','latex');
subplot(3,3,1);
% W1 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 0)./(kappa+gamma);
plot(Y(X_idx,:), W1(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=0}$','FontSize',18,'Interpreter','latex');
subplot(3,3,2);
% W2 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, pi/4)./(kappa+gamma);
plot(Y(X_idx,:), W2(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,3);
% W3 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, pi/2)./(kappa+gamma);
Y_idx = find(y >= 0, 1, 'first');
plot(Y(X_idx,:), W3(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=\pi/2}$','FontSize',18,'Interpreter','latex');
subplot(3,3,4);
% W4 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 3*pi/4)./(kappa+gamma);
X_idx = find(x >= 0, 1, 'first');
plot(Y(X_idx,:), W4(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=3\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,5);
% W5 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, pi)./(kappa+gamma);
plot(Y(X_idx,:), W5(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=\pi}$','FontSize',18,'Interpreter','latex');
subplot(3,3,6);
% W6 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 5*pi/4)./(kappa+gamma);
plot(Y(X_idx,:), W6(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=5\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,7);
% W7 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 3*pi/2)./(kappa+gamma);
plot(Y(X_idx,:), W7(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=3\pi/2}$','FontSize',18,'Interpreter','latex');
subplot(3,3,8);
% W8 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 7*pi/4)./(kappa+gamma);
plot(Y(X_idx,:), W8(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=7\pi/4}$','FontSize',18,'Interpreter','latex');
subplot(3,3,9);
% W9 = meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, 2*pi)./(kappa+gamma);
plot(Y(X_idx,:), W9(X_idx,:)./(kappa+gamma),'LineWidth',2);
axis tight;
xlabel('${\Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
title('${\varphi=2\pi}$','FontSize',18,'Interpreter','latex');


function W=meta_potential(eta, Delta, kappa, gamma, mu, alpha, b, theta, K, L, phi)
W = (abs(eta)).^2.*(Delta-2.*mu./(kappa+gamma).*alpha.^2+alpha.*cos(2.*theta)+...
    6*K.*(abs(eta)).^2+20*L.*(abs(eta)).^4)+2*abs(b).*abs(eta).*sqrt(kappa).*sin(phi-theta);
end