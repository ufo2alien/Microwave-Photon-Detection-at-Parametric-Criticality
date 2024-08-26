%% Load data
clear all
% load(['E:\Kirill\QWJPA_v2_2\09-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106091948\ProbeCharacter2us.mat']);
% load('E:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\PoissonFitDataPulseCharact.mat');
% load('E:\Kirill\QWJPA_v2_2\10-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106102039\ProbeCharacterization0.5us.mat');
% load('E:\Kirill\QWJPA_v2_2\10-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106102331\ProbeCharacterization0.25us.mat');

load('E:\Kirill\QWJPA_v2_2\11-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106111331\ProbeCharacter1us.mat');
%% Read data
dataOn=interp1([mean(HPProbeOff(:,:)) LPProbeOn(1,1:probePowerI)],[0 db2pow(probePower(1,1:probePowerI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)],(LPProbeOnStat(66,:)),'cubic');
% dataOff=interp1([mean(HPProbeOff(:,:)) LPProbeOn(1,1:probePowerI)],[0 db2pow(probePower(1,1:probePowerI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)],(LPProbeOffStat(1,:)),'cubic')
%% Non-necessary Recall
[v_mag,resnorm] = polyfit(dataxx,datayy,44);% lsqcurvefit(fit_func0,x0,dataxx-2.9761,datayy,[1e-3,0],[2,4],opts)

figure (738)
plot(dataxx,[datayy.' polyval(v_mag,dataxx,resnorm).'])
%% 1.Analysis If detector resolves PHOTON NUMBERs. Building calibration curve <n>(<I+iQ>)
coeff_n=1;
dataxx=[mean(LPProbeOff(:,:)) LPProbeOn(1,:)]*1.00;
datayy=[0 db2pow(probePower(1,:)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)]*coeff_n;%/0.787;
figure (3)
hold on
% [p,S]=polyfit(dataxx,datayy,20);
calCurve=interp1(dataxx,datayy,dataxx,'linear');%polyval(p,dataxx,S);%(dataxx,datayy,dataxx,'linear');
plot([mean(LPProbeOff(:,:)) LPProbeOn(1,:)],calCurve,'LineWidth',4)
xlabel('Detection Parameter - Mean of quantum Quadratures Sum $<|I+iQ|>$','interpreter','latex')
ylabel(['Mean photon number in ' num2str(probeLength*1e6) ' us pulse $<n>$'],'interpreter','latex')
 set(gca,'FontSize',18);
 grid on
%     fit_func0 = @(v,n)(v(1)*((n)).^v(2))
%     opts = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
%     opts.OptimalityTolerance=1e-4;
%     opts.MaxIterations=1600;
%     opts.FunctionTolerance=1e-4;
%     opts.StepTolerance=1e-4;
% %     ph=0;
%     x0 = [2;1]; 

%% Histograms based on previous PN-resolving calibration curve. Changing ind will change photon number in pulse.
figure (737)
clf
ind=118;
numb(ind)=[db2pow(probePower(1,ind)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)]*coeff_n;
dataOn=interp1(dataxx,datayy,(LPProbeOnStat(ind,:)),'pchip');
bar(0+(0:25).',[hist(dataOn,0.5+(0:25).','pdf').'/20000]);

var(suspiciousArray)
hold on
numb(ind)
plot((0:.2:25).',poisspdf(0:.2:25,(numb(ind))).','Linewidth',3)
grid on
xlabel('$<n>$','interpreter','latex')
ylabel('Probability','interpreter','latex')
set(gca,'FontSize',18);
ll=legend(['Experimental statistic with <' num2str((numb(ind) )) '> photon in pulse'],['Poisson distribution for ' num2str((numb(ind) )) ' photon ']);
set(ll,'FontSize',14)

%% Histograms for |I+iQ|
indMax=150;
clear co
phNum=[db2pow(probePower(1,:)-101.3-30+4.25)*Energy./(h_p.*PumpGenFreqCent/2)];
stepHist=0:0.2:18;
for ind=1:indMax
         co(ind,:)=histc(LPProbeOnStat(ind,:),stepHist)/length(LPProbeOnStat(ind,:));
end

figure(5578)
s=surf(1*phNum(1:ind),stepHist,co(:,:).');
set(s,'EdgeColor','none')
colormap hot
colorbar
ylabel('Detector value')
xlabel('<n> in 1 pulse')
view(0,90)
xlim([0 10])
%% Histograms for |I+iQ| dark

phNum=[db2pow(probePower(1,:)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)];
LPProbeOffStatMerged=[];
for ind=1:50
    LPProbeOffStatMerged=[LPProbeOffStatMerged LPProbeOffStat(ind,:)];
end

%% 

clear coDark stepHist
ThresVal=linspace(0,10,201);

for ind=1:indMax
for ThresValI=1:length(ThresVal)
    stepHist(ThresValI,:)=[min(ThresVal) ThresVal(ThresValI) max(ThresVal)];
      coDark(ind,ThresValI,:)=histc(LPProbeOffStat(ind,:),stepHist(ThresValI,:))/length(LPProbeOffStat(ind,:));
end
end
% figure(559)
% plot(ThresVal,coDark(1,:,2).');

%%
% bar(0+(stepHist).',histc(LPProbeOffStatMerged(1,:),stepHist)/length(LPProbeOffStatMerged(1,:)));
% title(['phNum: ', num2str(phNum(ind))]);

clear co
for ind=1:100
    co(ind,:)=histc(LPProbeOnStat(ind,:),stepHist)/length(LPProbeOnStat(ind,:));
end
    figure(556)
    
s=surf(phNum(1:ind),stepHist,co(:,:).');
%%

figure(559)
clf
plot(ThresVal,[coDark(1,:,2).' 1-coEff(:,2) ]);%(coEff(30,:,2).'-coDark(:,2))./(1-coDark(:,2))
%  hold on
%  [a]=fitdist(LPProbeOnStat(find(phNum>=1,1,'first'),:).','Weibull')
% plot(ThresVal,cdf('Weibull',ThresVal,a.a,a.b))

%% Poisson surface cumulative p0 dependence on <n>, threshold value
indMax=101;
phNum=[db2pow(probePower(1,:)-101.3-30+t_dur(3))*Energy./(h_p.*PumpGenFreqCent/2)];
coeff_n=1;%~1.2
coeff_nPois=1.0;%~2.4-2.6
clear co coEff coDark

ThresVal=linspace(0,10,201);

for ind=1:indMax
for ThresValI=1:length(ThresVal)
    stepHist(ThresValI,:)=[min(ThresVal) ThresVal(ThresValI) max(ThresVal)];
      coDark(ind,ThresValI,:)=histc(LPProbeOffStat(ind,:),stepHist(ThresValI,:))/length(LPProbeOffStat(ind,:));%p_dark negative
end
end

for ind=1:indMax
    for ThresValI=1:length(ThresVal)
         co(ind,ThresValI,:)=histc(LPProbeOnStat(ind,:),stepHist(ThresValI,:))/length(LPProbeOnStat(ind,:));
    end
end

% for ind=1:indMax

for ThresValI=1:length(ThresVal)

   coEff(ThresValI,:)=histc(LPProbeOnStat(find(phNum>=0.33,1,'first'),:),stepHist(ThresValI,:))/length(LPProbeOnStat(82,:));
%     coEff(ind,ThresValI,:)=histc(LPProbeOnStat(ind,:),stepHist(ThresValI,:))/length(LPProbeOnStat(82,:));
    
end
% end


% [X,Y]=meshgrid(exp(-(coeff_nPois.*phNum(1:ind))),-log(1-(TPR_Mean-FPR_Mean)./(1-FPR_Mean)));%Poisson
% [ZZ,Z]=meshgrid((phNum(1:ind)),1-(TPR_Mean-FPR_Mean)./(1-FPR_Mean));%Poisson

[X,Y]=meshgrid(phNum(1:ind),-log(1-(TPR_Mean-FPR_Mean)./(1-FPR_Mean)));%Bose-Einstein  
[ZZ,Z]=meshgrid((phNum(1:ind)),1-(TPR_Mean-FPR_Mean)./(1-FPR_Mean));%Bose-Einstein
% Z=squeeze(1-coEff(:,:,1)).';

    figure(556)
    clf
s=surf(1*phNum(1:ind),ThresVal(:),co(:,:,1).');
set(s,'EdgeColor','none')
set(s,'FaceAlpha',0.6)
view(0,90)
colormap hsv
zlim([0 1])
caxis([0 1])
freezeColors;


  figure(556)
%       clf
 hold on
% s=surf(phNum(1:ind),ThresVal(:),1-(d_pr+(1-X)).*Z+d_pr.*(1-X).*Y.^2);  
% s=surf(phNum(1:ind),ThresVal(:),X.*(1-squeeze(coDark(:,:,2)).'));
% s=surf(phNum(1:ind),ThresVal(:),1-((1-X.^(ZZ)).*(1-((Z).^(ZZ)))+(squeeze(coDark(:,:,2)).')-(squeeze(coDark(:,:,2)).').*(1-X.^(ZZ)).*(1-((Z).^(ZZ)))));
% s=surf(phNum(1:indMax),ThresVal(:),1-((1-X).*(1-(squeeze(coEff(:,:,1)).'))+(squeeze(coDark(:,:,2)).')-(squeeze(coDark(:,:,2)).').*(1-X).*(1-(squeeze(coEff(:,:,1)).'))));
% scont=contour3(phNum(1:ind),ThresVal(:),(X).^(Y.^(1./(0.115.*ZZ+1))).*(squeeze(coDark(:,:,1)).'),linspace(0,1,21),'Linewidth',3,'ShowText','off');%Poisson
scont=contour3(phNum(1:ind),ThresVal(:),((1+(X).*(Y.^(1./(0.12.*X+1)))).^-1).*(squeeze(coDark(:,:,1)).'),linspace(0,1,21),'Linewidth',3,'ShowText','off');%Bose-Einstein

clabel(scont,linspace(0,0.9,10)+0.05,'FontSize',14);
%  s=surf(phNum(1:ind),ThresVal(:),(X)); 
% s=surf(phNum(1:ind),ThresVal(:),X.*(1-Y));

% s=surf(phNum(1:ind),ThresVal(:),((1-Y).*(X)));
colormap hsv
view(0,90)
 set(gca,'FontSize',18,'FontName','Segoe UI Semibold','TickLabelInterpreter','latex');
% set(s,'FaceAlpha',0.6)
% set(s,'EdgeColor','none')
xlabel('$\bar{n}$','interpreter','latex')
ylabel('$R_\mathrm{th}$','Interpreter','latex')
clrbr=colorbar
 set(clrbr,'TickLabelInterpreter','latex')
% zlabel('$p(0)$','interpreter','latex')
xlim([0 floor(phNum(indMax))]);
%% Analysis, showing how variance of poisson distribution declines from its mean (declination from Poisson -> Bose-Einstein)
clearvars mleI mlvI numb dataOn dataOn_nonNeg
for ind=1:size(LPProbeOnStat,1)

dataOn(1,:)=interp1(dataxx,datayy,(LPProbeOnStat(ind,:)),'linear');%polyval(v_mag,LPProbeOnStat(ind,:),resnorm);%interp1([mean(LPProbeOff(:,:)) LPProbeOn(1,1:probePowerI)],[0 db2pow(probePower(1,1:probePowerI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)],(LPProbeOnStat(ind,:)),'cubic');
dataOn_nonNeg=dataOn(~isnan(dataOn(1,:)));
% counts=zeros(1,200);
% centers=zeros(1,200);
[counts,centers]=hist(dataOn,0+(0:100).','pdf');
suspiciousArray=[];
for suspI=1:length(counts)
    suspiciousArray=[suspiciousArray (suspI-1).*ones(1,counts(suspI))];
end
mlvI(ind)=(var(suspiciousArray));
mleI(ind)=mean(dataOn_nonNeg(1,:));
numb(ind)=[db2pow(probePower(1,ind)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)]*coeff_n;
end
figure(54630)
plot(numb,[mleI.' mlvI.'],'LineWidth',4);
grid on
set(gca,'FontSize',14);
xlabel('$\bar{n}$','interpreter','latex')
ylabel('$\lambda$ parameter (Poisson distr. MLE)','interpreter','latex')



%% 1st Method(Photon-Non-Resolving Detector) Constructing blocks from different <n> sequences, decomposition of each sequence to M blocks
clearvars StatBlocks NClicks
ind=1;
M=5e2;
ThresMean=2.47;
for nI=1:indMax
for blockI=1:M
StatBlocks(blockI,1:floor(length(LPProbeOnStat(nI,:))/M))=LPProbeOnStat(nI,blockI*(1:floor(length(LPProbeOnStat(nI,:))/M)));
NClicks(nI,blockI)=length(find(StatBlocks(blockI,1:floor(length(LPProbeOnStat(nI,:))/M))>=ThresMean));
end
end
%% Poisson&Bose-Einstein
maxNClicks=max(max(NClicks));
clearvars co cs meanMs varMs varM1 meanM1

indMax=101;
nbins=maxNClicks;%n_bins in histogram
Th=[0 .2 .7 1 1.2 1.5 1.8 2 2.5 3 3.5 4 5 6 7 8 9 10];
ThVal=[1 1 1 1 0.99 0.96 0.9 0.82 0.74 0.65 0.55 0.38 0.24 0.13 0.1 0.04 0.02 0.01];
phNum=[db2pow(probePower(1,1:indMax)-101.3+t_dur(3)-30)*Energy./(h_p.*PumpGenFreqCent/2)];
for nI=1:indMax
    coeff_x=-log(1-(TPR_Mean(find((ThresVal>=ThresMean),1,'first'))-FPR_Mean(find((ThresVal>=ThresMean),1,'first')))./(1-FPR_Mean(find((ThresVal>=ThresMean),1,'first'))));
    power=(coeff_x.^(1./(0.12.*phNum(nI)+1)));
Pr_uni=(1-1./(1+phNum(nI).*(power)))+ coDark(1,find((ThresVal>=ThresMean),1,'first'),2)- coDark(1,find((ThresVal>=ThresMean),1,'first'),2).*(1-1./(1+phNum(nI).*(power)));%Thermal    
%Pr_uni=(1-poisspdf(0,phNum(nI)).^(power))+
%coDark(1,find((ThresVal>=ThresMean),1,'first'),2)-
%coDark(1,find((ThresVal>=ThresMean),1,'first'),2).*(1-poisspdf(0,phNum(nI)).^(power));
%Poisson
[co(nI,:)]=hist(NClicks(nI,:),0:nbins)./(M);
[cs(nI,:)]=binopdf((0:maxNClicks),floor(length(LPProbeOnStat)./M),Pr_uni);
varM1(nI)=var(NClicks(nI,:));
meanM1(nI)=mean(NClicks(nI,:));
varMs(nI)=Pr_uni.*(1-Pr_uni).*length(LPProbeOnStat)./M;
meanMs(nI)=(Pr_uni).*length(LPProbeOnStat)./M;
end

figure (7777)
subplot(1,2,1)
s=surf(phNum.*length(LPProbeOnStat)/M,0:nbins,co(:,:).');
xlabel('$\bar{n}$ in block','interpreter','latex');
ylabel(['switching events in each block, ' num2str(length(LPProbeOnStat)/M) ' max'],'interpreter','latex')
title ('experiment','interpreter','latex')
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
clrbr=colorbar;
set(clrbr,'TickLabelInterpreter','latex')
set(s,'EdgeColor','none')
colormap hot
xlim([0 floor(phNum(end)).*length(LPProbeOnStat)/M]);
view(0,90)
figure (7777)
subplot(1,2,2)
s=surf(phNum.*length(LPProbeOnStat)/M,0:nbins,cs(:,:).');
xlim([0 floor(phNum(end)).*length(LPProbeOnStat)/M]);
xlabel('$\bar{n}$ in block','interpreter','latex');
ylabel(['switching events in each block, ' num2str(length(LPProbeOnStat)/M) ' max'],'interpreter','latex')
title ('analytics','interpreter','latex')
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
clrbr2=colorbar;
set(clrbr2,'TickLabelInterpreter','latex')
set(s,'EdgeColor','none')
colormap hot
view(0,90);

figure(6667)
% hold on
clf
mv=plot(phNum*length(LPProbeOnStat)/M,[meanM1.' varM1.' meanMs.' varMs.'],'linewidth',3);
xlabel('$\bar{n}$ in block','interpreter','latex');
ylabel('$E[X],E[X^2]$','interpreter','latex');
legend(mv,'$E[X]$, experiment','$E[X^2]$, experiment','$E[X]$, analytics','$E[X^2]$, analytics','interpreter','latex');
set(gca,'Fontsize',18,'TickLabelInterpreter','latex')
xlim([0 floor(phNum(end)).*length(LPProbeOnStat)/M]);

%% 2nd Method(Photon-Non-Resolving Detector. We take one sequence with <n> <= 1 and make decomposition over <n> by varying M(number of blocks number)
clearvars StatBlocks NClicks M
ind=1;%=find(phNum<=1.01,1,'last');
% M=1e2;
ThresMean=2.47;
M=floor(logspace(2.5,4,51));
% for M=1e2:100:1e3
for mI=1:length(M)
    for blockI=1:M(mI)
        StatBlocks(blockI,1:int16(floor(length(LPProbeOnStat(ind,:))/M(mI))))=LPProbeOnStat(ind,(blockI-1)*int16(floor(length(LPProbeOnStat(ind,:))/M(mI)))+1:blockI*int16(floor(length(LPProbeOnStat(ind,:))/M(mI))));
        NClicks(mI,blockI)=length(find(StatBlocks(blockI,1:int16(floor(length(LPProbeOnStat(ind,:))/M(mI))))>=ThresMean));
    end
end

%%

maxNClicks=max(max(NClicks));
clearvars co cs coeff_x power Pr_uni

% ind=76
% coeff_n=0.92;
% phNum=[db2pow(probePower(1,1:indMax)-101.3+t_dur(3)-30)*Energy./(h_p.*PumpGenFreqCent/2)];
phNum=0.1;
ind=find(phNum<=1.01,1,'last');
coeff_x=-log(1-(TPR_Mean(find((ThresVal>=ThresMean),1,'first'))-FPR_Mean(find((ThresVal>=ThresMean),1,'first')))./(1-FPR_Mean(find((ThresVal>=ThresMean),1,'first'))));
power=(coeff_x.^(1./(0.115.*phNum(ind)+1)));
Pr_uni=(1-poisspdf(0,phNum(ind)).^(power))+ coDark(1,find((ThresVal>=ThresMean),1,'first'),2)- coDark(1,find((ThresVal>=ThresMean),1,'first'),2).*(1-poisspdf(0,phNum(ind)).^(power));
   

for mI=1:length(M)
   
    
[co(mI,:)]=hist(NClicks(mI,1:M(mI)),(0:maxNClicks))./M(mI);
[cs(mI,:)]=binopdf((0:maxNClicks),floor(length(LPProbeOnStat)./M(mI)),Pr_uni);
%co(mI,:)=co(mI,:)/max(co(mI,:) );


meanM(mI)=sum((((find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))-1).*(co(mI,(find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))));%mean(NClicks(mI,1:M(mI)));%sum((((find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))-1).*(co(mI,(find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))));%
varM(mI)=sum((((find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))-1-mean(NClicks(mI,1:M(mI)))).^2.*(co(mI,(find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))));%var(NClicks(mI,1:M(mI)));var(NClicks(mI,1:M(mI)));
meanMs(mI)=sum((((find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))-1).*(cs(mI,(find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))));%mean(NClicks(mI,1:M(mI)));%sum((((find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))-1).*(co(mI,(find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))));%
varMs(mI)=sum((((find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))-1-meanMs(mI)).^2.*(cs(mI,(find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))));%var(NClicks(mI,1:M(mI)));var(NClicks(mI,1:M(mI)));

end
%%
figure(747)
subplot(1,2,1)
s=surf(phNum(ind).*length(LPProbeOnStat)./M(1:mI),0:maxNClicks,(co(:,:)).');

colorbar;
set(s,'EdgeColor','none')
colormap hot
xlim([(phNum(ind).*length(LPProbeOnStat)/M(end)) floor(phNum(ind).*length(LPProbeOnStat)/M(1))]);
view(0,90)

xlabel('$\left<n\right>$ in block','interpreter','latex');
ylabel(['switching events in each block'],'interpreter','latex')
colorbar;
set(s,'EdgeColor','none')
 set(gca,'FontSize',16);
xlim([(phNum(ind).*length(LPProbeOnStat)/M(end)) floor(phNum(ind).*length(LPProbeOnStat)/M(1))]);
figure(747)
subplot(1,2,2)
s=surf(phNum(ind).*floor(length(LPProbeOnStat)./M(1:mI)),0:maxNClicks,(cs(:,:)).');
xlabel('$\left<n\right>$ in block','interpreter','latex');
ylabel(['switching events in each block'],'interpreter','latex')
 set(gca,'FontSize',16);
 set(s,'EdgeColor','none')
view(0,90)
xlim([(phNum(ind).*length(LPProbeOnStat)/M(end))  floor(phNum(ind).*length(LPProbeOnStat)/M(1))]);
colorbar;
figure(666)
% hold on
mv=plot(phNum(ind).*length(LPProbeOnStat(ind,:))./M(1:mI),[meanM(1:mI).' varM(1:mI).' meanMs(1:mI).' varMs(1:mI).'],'Linewidth',3)
xlim([(phNum(ind).*length(LPProbeOnStat)/M(end))  floor((phNum(ind)).*length(LPProbeOnStat)/M(1))]);
xlabel('$\left<n\right>$ in block','interpreter','latex');
legend(mv,'Mean, experiment','Variance, experiment','Mean, analytics','Variance, analytics');
 set(gca,'FontSize',16);
grid on
%% 2nd Method(Statistics for state with n~0, DARK COUNTS)
clearvars StatBlocks NClicks M LPProbeOffStatMerged
% ind=40;
% M=1e2;
M=logspace(2.5,4.5,51);
LPProbeOffStatMerged=[];
for ind=1:indMax
    LPProbeOffStatMerged=[LPProbeOffStatMerged LPProbeOffStat(ind,:)];
end
% for M=1e2:100:1e3
StatBlocks=zeros(length(M));
NClicks=zeros(length(M));
ind=1
for mI=1:length(M)
    for blockI=1:M(mI)
        StatBlocks(blockI,1:floor(length(LPProbeOffStatMerged(ind,:))/M(mI)))=LPProbeOffStatMerged(ind,((blockI-1)*floor(length(LPProbeOffStatMerged(ind,:))./M(mI))+1):blockI*floor(length(LPProbeOffStatMerged(ind,:))./M(mI)));
        NClicks(mI,blockI)=length(find(StatBlocks(blockI,1:floor(length(LPProbeOffStatMerged(ind,:))./M(mI)))>=ThresMean));
    end
end

%%

maxNClicks=max(max(NClicks));
clearvars co meanM varM

% ind=76
phNum=1*[db2pow(probePower(1,ind)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)]
for mI=1:length(M)
[co(mI,:)]=hist(NClicks(mI,1:floor(M(mI))),(0:maxNClicks))./(M(mI));%./sum(NClicks(mI,1:M(mI)));
[cs(mI,:)]=binopdf((0:maxNClicks),floor(length(LPProbeOffStatMerged)./M(mI)),0.10);
%co(mI,:)=co(mI,:)/max(co(mI,:) );
varM(mI)=var(NClicks(mI,1:floor(M(mI))));
meanM(mI)=mean(NClicks(mI,1:floor(M(mI))));
meanMs(mI)=sum((((find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))-1).*(cs(mI,(find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))));%mean(NClicks(mI,1:M(mI)));%sum((((find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))-1).*(co(mI,(find(co(mI,:),1,'first')):(find(co(mI,:),1,'last')))));%
varMs(mI)=sum((((find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))-1-meanMs(mI)).^2.*(cs(mI,(find(cs(mI,:),1,'first')):(find(cs(mI,:),1,'last')))));%var(NClicks(mI,1:M(mI)));var(NClicks(mI,1:M(mI)));
end
figure(747)
subplot(1,2,1)
s=surf(floor(length(LPProbeOffStatMerged)./M(1:mI)),0:maxNClicks,(co(:,:)).');


% set(s,'EdgeColor','none')

figure(747)
s=surf(floor(length(LPProbeOffStatMerged)./M(1:mI)),0:maxNClicks,(cs(:,:)).');
subplot(1,2,2)
xlabel('N/M, number of pulses in each block','interpreter','latex');
ylabel(['switching events in each block'],'interpreter','latex')
% title (['number of blocks M = ' num2str(M(mI))])
colorbar;
colormap hot
xlim([0 floor(phNum(ind)).*length(LPProbeOnStat)/M(1)]);
figure(666)
% hold on
plot(length(LPProbeOffStatMerged)./M(:),[meanM.' varM.' meanMs.' varMs.'])
xlabel('N/M, number of pulses in each block','interpreter','latex');
grid on
pp=polyfit(floor(length(LPProbeOffStatMerged)./M(:)),meanM.',1)

%%
for ThresValI=1:length(ThresVal)

   coEff(ThresValI,:)=histc(LPProbeOnStat(find(phNum>=1,1,'first'),:),stepHist(ThresValI,:))/length(LPProbeOnStat(82,:));
%     coEff(ind,ThresValI,:)=histc(LPProbeOnStat(ind,:),stepHist(ThresValI,:))/length(LPProbeOnStat(82,:));
    
end
 plot(ThresVal,[(1-exp(-exp(log(3)-0.33*(ThresVal-0.9)).'/3)) (coEff(:,2))/1.4])