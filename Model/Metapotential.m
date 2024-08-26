


% SPD paper
% syms g I Q mu
% % alpha - pump strength, V - effective potential, Delta - detuning,
% % K - Kerr nonlinearity, L - sextic nonlinearity, eta - coherent state,
% % phi - phase angle of eta
syms alpha V eta Delta K phi L W eta_plus_squared eta_minus_squared V_depth alpha_0 alpha_t alpha_c
alpha = 0.95*(coup+losses)/2;%1.000-2.995e-1; % 2 or 6 or 5 or 46
Delta =0.7e6; % 6 or -2 or 8 or 50
K = -1.79e3; % -0.01
L = 3.818e-1; % 1E-5
x = linspace(-50, 50, 301); % -15 to 15
y = linspace(-50, 50, 301); % -15 to 15
[X, Y] = ndgrid(x, y);
eta = X + 1i*Y;
phi = angle(eta);
V  = (abs(eta)).^2.*(Delta+((coup+losses)/2-alpha.*cos(2.*phi))+6*K.*(abs(eta)).^2+20*L.*(abs(eta)).^4); % Add a Pi phase to the argument of cosine function.
alpha_0 = Delta-3*K^2/(5*L); %% Zero-amplitude pumping strength
alpha_t = Delta-9*K^2/(20*L); %% Transition pumping strength
alpha_c = Delta; %% Critical pumping strength
grid on
% W = (abs(eta)).^2.*(Delta+alpha.*cos(2.*phi)+K*alpha^2.*(abs(eta)).^2.*cos(4.*phi)+L*alpha^3.*(abs(eta)).^4.*cos(6.*phi));

% % Surface plot of effective potential V
% sgtitle(['$\alpha$ = ',num2str(alpha),', ','$\Delta$ = ',num2str(Delta), ',','$V_{depth} = $',num2str(V_depth)], 'Interpreter','latex');
figure (248)
subplot(2,2,1);
surf(real(eta),imag(eta),V);
axis tight;
set(findobj('Type','Surface'),'EdgeColor','none')
% set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
zlabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% % Real axis
subplot(2,2,2);
% figure(2);
Y_idx = find(y >= 0, 1, 'first');
plot(X(:,Y_idx), V(:,Y_idx),'LineWidth',2);
axis tight;
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
set(gca,'XTick',[], 'YTick', []);
% % Imaginary axis
subplot(2,2,3);
% figure(3);
X_idx = find(x >= 0, 1, 'first');
plot(Y(X_idx,:), V(X_idx,:),'LineWidth',2);
axis tight;
xlabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% ylim([-3,3])
% set(gca,'XTick',[], 'YTick', []);
subplot(2,2,4);
% figure(4);
contour(real(eta),imag(eta),V,50,'ShowText','on');
axis tight;
% set(gca,'XTick',[], 'YTick', [], 'ZTick',[]);
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${Im(\eta)}$','FontSize',18,'Interpreter','latex');
zlabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');

eta_plus_squared = 1/2*(-(K)/(5*L)+sqrt(K^2/(25*L^2)-1/(15*L)*(Delta-alpha))); %% Minima/Wells; location of minima, i.e., sqrt(eta_plus_squared), corresponds to coherent oscillation amplitude
eta_minus_squared = 1/2*(-(K)/(5*L)-sqrt(K^2/(25*L^2)-1/(15*L)*(Delta-alpha))); %% Maxima/Hills
V_depth = (Delta-alpha)*eta_minus_squared+6*K.*(eta_minus_squared).^2+20*L.*(eta_minus_squared).^3; %% Excitation energy
sgtitle(['$\alpha$ = ',num2str(alpha),', ', '$\alpha_0$ = ',num2str(alpha_0),', ', '$\alpha_t$ = ',num2str(alpha_t),', ','$\Delta$ = ',num2str(Delta), ', ','$V_{depth}$ = ',num2str(V_depth), ', ','Maxima at $\pm$',num2str(sqrt(eta_minus_squared)), ', ','Minima at $\pm$',num2str(sqrt(eta_plus_squared))], 'Interpreter','latex');
 figure(2);
Y_idx = find(y >= 0, 1, 'first');
plot(X(:,Y_idx), V(:,Y_idx),'LineWidth',2);
axis tight;
xlabel('${Re(\eta)}$','FontSize',18,'Interpreter','latex');
ylabel('${V(\eta)}$','FontSize',18,'Interpreter','latex');
% figure(2);
% surf(real(eta),imag(eta),W);
