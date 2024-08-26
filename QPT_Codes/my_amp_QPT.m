
%% Data loading and fitting
clear all;

cd('C:\Users\petrovk2\YandexDisk\AALTO_PROJECTS\QWJPA_v2\');%CHANGE PATHS AS YOU NEED
% addpath('C:\Users\petrovk2\YandexDisk\MATLAB');
addpath('C:\Users\petrovk2\YandexDisk\AALTO_PROJECTS\functions_for_calc');
pathToFile = ['C:\Users\petrovk2\Dropbox (Aalto)\temp\Gain_12.084202104232356GHz.mat'];

ang = 0.230;
offs = -0.52e6;
% offs=0.893e6;
% offs = 0.893e6;
[v_mag,f,om0,mag,pha,Gain]=BasicFitting(pathToFile,offs,ang,0);%Use for estimations of coupling rates; kappa related to signal port, gamma - linear dissipation port

%%
clearvars -except coup losses detun magn pha v_mag mag freqs_exper pows_exper pump Gain om0 ph SignalGeneratorPumpPt EIG varI_ave varQ_ave ang offs;
coup = real((v_mag(1))); % kappa coupling rate
losses = real((v_mag(2))); % gamma coupling rate
k = (coup+losses); % total coupling rate
detun=linspace(-6e6,6e6,21); % some imaginary detuning value
Amp = [0 linspace(0.3,0.75,41)*real(k)];   [ (linspace(0,0.4,21)) (linspace(0.41,0.65,21)) linspace(0.66,1,21)].*1*real(k);%Normalized pump Amplitude Sweep
fpoints = 2e3;%Number of spectral points
fstart = -5e6;
fend = 5e6;
oml = linspace(fstart,fend,fpoints);% omega/2pi
om = 2*pi*oml; 
K = -1.5e4/12;%5.5e0.*om0/12;%8e-6;1e01.65e-1*om0; %Kerr coefficient
L=-1.5e-5*K;
gamma2 = 0.*K; %Coupling to nonlinear dissipation port

if isempty(gcp('nocreate'))
    poolobj = parpool('local', 4); 
end

samplerate = (fend-fstart);
h_p=6.62e-34;
test_A = sqrt(coup);%6.34e3; %If we have non-cold case test_A=sqrt(coth(hw/2kT));
test_B= sqrt(losses);
t_s = (0:(fpoints-1)).'./samplerate; %time trace
a_in = test_A.*1/2.*(randn(1,fpoints).'+1i.*randn(1,fpoints).');
b_in = test_B.*1/2.*(randn(1,fpoints).'+1i.*randn(1,fpoints).'); 
c_in = test_A.*1/2.*(randn(1,fpoints).'+1i.*randn(1,fpoints).'); 
f_s = linspace(-samplerate/2,samplerate/2,length(t_s)).';

a = zeros(fpoints,length(Amp));
a_out = zeros(fpoints,length(Amp));

que = parallel.pool.DataQueue;
    % afterEach(que, fff=@(p) (p+1));
m=length(Amp);
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,m) '\n\n']);
fscan=fstart:(fend-fstart)/fpoints-1:fend;
tpoints=fpoints;
% opts = odeset('RelTol',1.0e-4,'AbsTol',1.5e-6);
f1 = -1e6;
f2 = 1e6;


opts = odeset('RelTol',3e-3,'AbsTol',1e-4);
for detun_i=1:length(detun)
parfor i=1:length(Amp)
    [t,ySol] = ode45(@(t,y) myode1_nl_2order(t,y,(k),(coup),K,L,Amp(i),gamma2.*0,a_in,b_in,t_s,detun(detun_i)),t_s,[0],opts);%calculates QLE for intracavity mode a(t) *exp(-1i*pi/1.3)
    a(:,i) = ySol(:,1);%+ySol(:,2);
    a_out(:,i)=(a_in(:,1)-sqrt(coup)*(a(:,i)))/sqrt(coup);%Input-output formalism
    n_res(detun_i,i)=mean(abs(a(:,i)).^2)*(k).^2;
    
    A(detun_i,i,:)=a(:,i);
    A_out(detun_i,i,:)=a_out(:,i);
    fprintf('\b|\n');
    A_out_f(detun_i,i,:) = sqrt(samplerate/(f2-f1)).*IQfiltering(squeeze(A_out(detun_i,i,:)),samplerate,0,1,[f1 f2]);
    n_out(detun_i,i)=mean(abs(A_out_f(detun_i,i,:)).^2);
    var_A(detun_i,i)=var(A_out_f(detun_i,i,:)/sqrt(2));
    mean_A(detun_i,i)=mean(A_out_f(detun_i,i,:)/sqrt(2));
end
disp(['Variance of a_out: ',  num2str(var((a_out(1:end,1))))])

end

%%

figure(334)
surf(detun,Amp/k,n_out.') 
set(findobj('Type','Surface'),'EdgeColor','none')
view(0,90);colorbar
xlabel('${\Delta}$','FontSize',18,'Interpreter','latex')
ylabel('${\alpha}/(\kappa+\gamma)$','FontSize',18,'Interpreter','latex')
ylim([0.3 0.75]);
figure(335)
subplot(2,1,1)
surf(detun,Amp/k,abs(mean_A.'))
set(findobj('Type','Surface'),'EdgeColor','none')
view(0,90);colorbar
xlabel('${\Delta}$','FontSize',18,'Interpreter','latex')
ylabel('${\alpha}/(\kappa+\gamma)$','FontSize',18,'Interpreter','latex')
ylim([0.3 0.75]);
subplot(2,1,2)
surf(detun,Amp/k,var_A.') 
set(findobj('Type','Surface'),'EdgeColor','none')
view(0,90);colorbar
xlabel('${\Delta}$','FontSize',18,'Interpreter','latex')
ylabel('${\alpha}/(\kappa+\gamma)$','FontSize',18,'Interpreter','latex')
ylim([0.3 0.75]);
%%
for detun_i=1:length(detun)
    for i=1:length(Amp)
        Al_clrit(detun_i,i)=sqrt(k^2/4+(detun(detun_i)+12.*K.*n_res(detun_i,i)).^2);
    end
end
figure (3542)
clf
surf(detun,Amp,Al_clrit.')
%% Histograms-Correlations-Squeezing
%commented below regarded case calculating a_dag
clearvars CovMat N_A L I Q PPTcrit0 B1 B11 B22 B2 a_t_f
sc=1;%0.66
% samplerate = fend-fstart;
f1 = -1e6;
f2 = 1e6;
detun_i=8;
for pow_p =41%1:round(length(Amp)*sc)

    fcut = 1;
    [a_t_f]=sqrt(samplerate/(f2-fcut)).*IQfiltering(squeeze(A_out(detun_i,pow_p,:)),samplerate,0,2,[f1 -fcut],[fcut f2]);%.*exp(0.6*-1i.*2.*pi.*detun/4.*t_s))/test_A)
%     [a_t_f_dag]=sqrt(2*fend/f2).*IQfiltering(a_out_dag(:,pow_p),samplerate,0,2,[f1 -fcut],[fcut f2]);
     B1=a_t_f(1,:).';
     B2=(a_t_f(2,:).');
%      B1d=a_t_f_dag(2,:).';
%      B2d=(a_t_f_dag(1,:).');
    for dd=1%Calculation of the optimal IQ phase
         L=[real((B1(:,1))), imag((B1(:,1))), real(B2(:,1)), imag((B2(:,1)))];
    %      L=[real((B1(:,1)+B1d)/2), real((B1(:,1)-B1d)/1/1i), real(B2(:,1)+B2d)/2, real((B2(:,1)-B2d)/2i)];
         CovMat=cov(L);
         C21=CovMat(1:2,3:4);
         if (C21(1,1)-C21(2,2))<0
            B11(:,1)=-B1(:,1).*exp(-1i*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2))));
            IQphi=1i*(pi-atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2))));
         else
            B11(:,1)=B1(:,1).*exp(-1i*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2))));
            IQphi=-1i.*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2)));
         end
         B22(:,1)=B2(:,1);
    end

% L1=[((B1(:,1)+B1d)/2), ((B1(:,1)-B1d)/1/1i), (B2(:,1)+B2d)/2, ((B2(:,1)-B2d)/1/1i)];

     L1=[real((B11(:,1))), imag((B11(:,1))), real(B22(:,1)), imag((B22(:,1)))];
     CovMat=4*cov(L1);
     if pow_p==1
         CovMat_init=CovMat;
%          CovMat=(CovMat-CovMat_init)+0*eye(4);
     elseif pow_p~=1
         CovMat=(CovMat-0*CovMat_init)+0*eye(4);%/min(min(diag(CovMat_init)));%-CovMat_init)+eye(4)
     end

     figure (555);
     clf;
     subplot(1,1,1);

     bar3(CovMat);
     view(30, 15);
     set(gca,'XTickLabel',{'I_1','Q_1','I_2','Q_2'});
     set(gca,'YTickLabel',{'I_1','Q_1','I_2','Q_2'});

     omeg=[0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0];                  
     PPTcrit(pow_p) = min(abs(eig(1i*omeg*diag([1 1 1 -1])*CovMat*diag([1 1 1 -1]))));
     N_A(pow_p) =min(eig(CovMat));                  
     title(['\nu_{eig}: ' num2str(N_A(pow_p))]);
      
%    1-mode Squeezing
     IQ_t = sqrt(samplerate/(f2-f1)).*IQfiltering(squeeze(A_out(detun_i,pow_p,:)),samplerate,0,1,[f1 f2]).'.*exp(IQphi/2);
     AA = [-1:0.01:1].*max(abs(IQ_t));
     [MM CC] =hist3([real(IQ_t) imag(IQ_t)], {AA AA});
     varI(pow_p)=var(real(IQ_t));
     varQ(pow_p)=var(imag(IQ_t));
%      [MM CC] =hist3([(a_out(:,pow_p)+a_out_dag(:,pow_p))/2-mean((a_out(:,pow_p)+a_out_dag(:,pow_p))/2) (a_out(:,pow_p)-a_out_dag(:,pow_p))/2i-mean((a_out(:,pow_p)-a_out_dag(:,pow_p))/2i)], {AA AA});
%      varI(pow_p)=var(a_out(:,pow_p)/2+a_out_dag(:,pow_p)/2);
%      varQ(pow_p)=var((a_out(:,pow_p)-a_out_dag(:,pow_p))/2i);
     figure (5551);
     clf;
%      subplot (1,2,1); 
     surface(AA, AA, MM./length(IQ_t), 'edgecolor', 'none')

     axis square
     colorbar
     title(['\alpha/\kappa: '  num2str(Amp(pow_p)/k)]);
     
%    M(pow_p)=getframe(gcf); %For video capturing
end

%% PLOT FIGURES
p_atten=0.45;
figure (123121);
clf;

plot(abs(Amp(1:round(length(Amp).*sc)))./k,[N_A(1:round(length(Amp).*sc))'],'LineWidth',3);
 load('C:\Users\petrovk2\Dropbox (Aalto)\Multimode_DataBlanks\2modeEnt&Squeezing\2modeSim&ExperData.mat', 'pump','EIG','varI_ave','varQ_ave','varI_ave_init','varQ_ave_init','h_p','MeasFreq','RBW','samplerate','PumpGenFreqCent','RBW')

% load('C:\Users\petrovk2\YandexDisk\AALTO_PROJECTS\JPA\1_pump_data\Recalculated_old_main_21_0.4175.mat', 'SignalGeneratorPumpPt','EIG','varI_ave','varQ_ave','varI_ave_init', 'varQ_ave_init','h_p','MeasFreq');
%  load('D:\scripts\02-Nov-2020\RECALCULATED_20_0.4847.mat', 'SignalGeneratorPumpPt','EIG','varI_ave','varQ_ave','h_p')
figure (123121)
hold on
% addpath('C:\Users\petrovk2\Dropbox (Aalto)\PetrovninK\functions_for_calc')
pump1=pumpCalibration(pump);
plot(pump.'/p_atten,(EIG.'+0.00),'o','LineWidth',3)
ylim([-0 12]);
xlim([-0 1]);
ylabel ('Symplectic eigenvalue, $\min\lbrace\nu_{i}\rbrace$','interpreter','latex');
xlabel('Normalized pump amplitude, $\alpha/(\kappa+\gamma)$','interpreter','latex');
plot(xlim,[1 1],'--','Color','g','LineWidth',3);
lel1=legend('Simulation','Experiment','Quantum entanglement threshold');
set(lel1,'interpreter','latex');
grid on

figure (232141)
clf;

plot(abs(Amp(1:round(length(Amp).*sc)))./k,pow2db(abs(1+4*[(varI(1:round(length(Amp).*sc))-varI(1)).', (varQ(1:round(length(Amp).*sc))-varQ(1)).'])),'LineWidth',3);%Amp(2:round(length(Amp).*sc))
hold on
plot(pump.'/p_atten,pow2db([4*(varI_ave'-varI_ave_init.')+1, 4*(varQ_ave.'-varQ_ave_init.'+0.0)+1.0]),'o','LineWidth',3)
ylim([-5 40]);%-varI(1)+1;-varQ(1)+1
xlim([0 1]);
plot(xlim,[0 0],'--','Color','g','LineWidth',2);
grid on
lel=legend('$\left< \Delta Q^2\right>$, simulation','$\left< \Delta I^2\right>$, simulation','$\left< \Delta Q^2\right>$, experiment','$\left< \Delta I^2\right>$, experiment','Quantum entanglement threshold','Intepreter','latex');
set(lel,'interpreter','latex');
xlabel('Normalized pump amplitude, $\alpha/(\kappa+\gamma)$','interpreter','latex');
ylabel ('Squeezing (dB), $\nu_{min}$','interpreter','latex');


%% VIDEO WRITER
% v = VideoWriter('squeezing-covMat-linear1','MPEG-4');
% open(v)
% writeVideo(v,M)
% close(v)
%% PLOT S11
% p_atten=-78
% figure (31)
% clf;
% surf_exp_mag=surf((freqs_exper)+0*offs,sqrt(db2pow(pows_exper-30+p_atten)./1.054e-34/om0)*1/abs(v_mag(1)+v_mag(2)),mag2db(abs((1+exp(-1i*ang*pi).*(-1+(db2mag(mag(:,2:end)-mag(:,1)).').*exp(+1i.*(pha.'+ph)))))));%-73.5
% surf_exp_mag.EdgeColor='none';
% colormap jet;
% % freezeColors;
% 
% figure(31);
% %  clf;
% xlim([1*fstart 1*fend]+0*offs);
% hold on
% s=surf(f_s.',(Amp(1:end)/abs(k)),mag2db(abs(aw(:,1:end).')));%+flip(aw_dag(:,2:end)',2
% s.EdgeColor='none';
% colormap hot;
% hold off
% ylim([Amp(1) Amp(end)]./abs(k));
% 
% figure (1242112517);
% clf;
% plot ((Amp/k),mag2db(abs(aw(ceil(length(aw)/2),1:end).')))
% hold on
% plot(sqrt(db2pow(pows_exper-30+p_atten)./1.054e-34/om0)*1/real(v_mag(1)+v_mag(2)),mag2db(abs((magn(501,:).*exp(+1i.*(pha(501,:)+ph))))))
% % xlim([0 1]);
%%
f1 = -1e6;
f2 = 1e6;
for detun_i=1:length(detun)
for i=1:length(Amp)
    A_out_f(detun_i,i,:) = sqrt(samplerate/(f2-f1)).*IQfiltering(squeeze(A_out(detun_i,i,:)),samplerate,0,1,[f1 f2]);
    n_out(detun_i,i)=mean(abs(A_out_f(detun_i,i,:)).^2);
    var_A(detun_i,i)=var(A_out_f(detun_i,i,:));
    mean_A(detun_i,i)=mean(A_out_f(detun_i,i,:));
end
detun_i

end
