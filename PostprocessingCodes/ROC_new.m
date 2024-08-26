%% Load data
clear all
% load(['D:\Kirill\QWJPA_v2_2\09-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106091948\ProbeCharacter2us.mat']);
load('D:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\PoissonFitDataPulseCharact.mat');
% load('D:\Kirill\QWJPA_v2_2\10-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106102039\ProbeCharacterization0.5us.mat');
% load('D:\Kirill\QWJPA_v2_2\10-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106102331\ProbeCharacterization0.25us.mat');
%%
clear TPR_Mean FPR_Mean ClickProbMeanTrue ClickProbMeanFalse
t_dur=[2.4 3.67 4.25 4.97];
phNum=[db2pow(probePower(1,:)-101.3+t_dur(3)-30)*Energy./(h_p.*PumpGenFreqCent/2)];
probePowerI=find(phNum<=1.01,1,'last');
ThresVal=linspace(0,10,2001);
phNum(probePowerI)
for ThresValI=1:length(ThresVal)
    ClickMeanTrue(ThresValI)=0;
    ClickMeanFalse(ThresValI)=0;
for cycle_i=1:N_cycles

     
      if LPProbeOnStat(probePowerI,cycle_i)>=ThresVal(ThresValI)
          ClickMeanTrue(ThresValI)=ClickMeanTrue(ThresValI)+1;
      end
      
     if LPProbeOffStat(probePowerI,cycle_i)>=ThresVal(ThresValI)
          ClickMeanFalse(ThresValI)=ClickMeanFalse(ThresValI)+1;
     end

  
    end
    ClickProbMeanTrue(ThresValI)=ClickMeanTrue(ThresValI)/N_cycles;
    ClickProbMeanFalse(ThresValI)=ClickMeanFalse(ThresValI)/N_cycles;    
end
TPR_Mean=ClickProbMeanTrue;
FPR_Mean=ClickProbMeanFalse;
AUC_mean=trapz(flip(FPR_Mean),flip(TPR_Mean));
%% Optimal thresVal plot
ROCfig=figure(11432);

hold on
plot(FPR_Mean(1,1:end),TPR_Mean(1,1:end),'LineWidth',4);
title(['ROC - Detection by Mean. AUC: ' num2str(AUC_mean)]);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
xlim([0 1]);
ylim([0 1]);
% hold on
% plot([0 1], [0 1],'LineWidth',3);
% grid on

figure (431242);
% hold on
plot(ThresVal(:),TPR_Mean(:).*(1-FPR_Mean(:)),'Linewidth',6)
ylabel('TPR (1-FPR)','interpreter','latex')
set(gca,'FontSize',28);
xlabel('$\left<|I+iQ|\right>$ threshold value','interpreter','latex')