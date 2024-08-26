fa=6.042e9;
kb = 1.380649e-23;
Teff=20e-3;
h_p = 6.62607004e-34;
kappa=(6.7e6);
Delta=linspace(0,10e6,101);
fp=1.05;
mu_b=sqrt(fp^2-1);
dU=linspace(0,1,101);mu_b.*(Delta/kappa-mu_b).^2;
tau=20e-6;
Gamma=fa.*exp(-dU/coth((h_p*fa)/(kb*Teff)));
Psw=1-exp(-Gamma*tau);
ss=figure(135)
hold on
plot(dU,Psw)
%%
% fit_func0 = @(C,f)(1-exp(-tau.*C(1).*exp(-C(2).*f.^2/coth((h_p*fa)/(kb*Teff)))));
% 
%     opts = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
%     opts.OptimalityTolerance=1e-10;
%     opts.MaxIterations=30000;
%     opts.FunctionTolerance=1e-10;
%     opts.StepTolerance=1e-10;
% %     ph=0;
%     x0 = [1e9;1e-9]; 
%     [CC,resnorm] = lsqcurvefit(fit_func0,x0,f,prob20mK.',[],[])
%     ss=figure(1235)
%     clf
%     plot(f,prob20mK);
%     hold on
%     plot(f,fit_func0(CC,f));
%% plotting Psw
clearvars Psw
Teff=[105e-3 240e-3 330e-3 410e-3 480e-3 550e-3 620e-3];%[125e-3 250e-3 320e-3 390e-3 460e-3 520e-3 580e-3];
C1=1;
C2=0.1e6;0.295e6;
C3=1.01;
Delta=linspace(0,10e6,1001);
k=2.55e7*1.38;
K=0.23e3*2*pi;
tmeas=20e-6;
% fig=[];
if ~isempty(fig)
close(fig)
end
% fig=figure
fig=openfig('C:\Users\KVANTTIadmin\Dropbox (Aalto)\SPDetectorJPA\PswFreqTemp\Psw_vs_freq_detuning_no20mK.fig','reuse');
for T_i=1:length(Teff)
    fp=sqrt(C3^2-1);
    PRE=1;fp.*sqrt((2*pi*Delta/k).^2-fp.^2);
%     Gamma=10*Delta*fp./(3*fa);%ZORIN
%     Teffa(T_i)=(1*h_p.*fa^2.*coth((h_p*fa)/(2*kb.*Teff(T_i))).*pi)./(4*50.*(4e-6.*cos(.375)).^2);%ZORIN
%     Psw(T_i,:)=1-exp(-tmeas*fa.*2./sqrt(2).*2*pi*2*Delta/fa.*fp.*exp(-Gamma./Teffa(T_i)));%ZORIN
Psw(T_i,:)=1-exp(-C1*tmeas*k.*PRE.*exp(-fp.*(2*pi*Delta/k-fp).^2./(2*coth((h_p*fa)/(2*kb.*Teff(T_i)))*(6*K/k))));%coth((h_p*fa)/(kb.*Teff(T_i)))(1+(h_p*fa).^-1/(kb.*Teff(T_i)).^-1))  (1+(kb.*Teff(T_i))./(h_p*fa))

hold on
plot(Delta+0e6,Psw(T_i,:),'--','Linewidth',4)

end
legend([num2str(Teff.') repmat(' K',length(Teff),1)]);
hold off
%% dU
figure (556)
clearvars dU
% for T_i=1:length(Teff)
dU(T_i,:)=fp.*(2*pi*Delta/k-fp).^2/4;%(1+(h_p*fa).^-1/(kb.*Teff(T_i)).^-1))  (1+(kb.*Teff(T_i))./(h_p*fa))
% hold on
plot(Delta+0.0e6,dU(T_i,:),'--')
% end
xlabel('Frequency detuning','interpreter','latex')
ylabel('barrier height, dU, in $\bar{n}$','interpreter','latex')
