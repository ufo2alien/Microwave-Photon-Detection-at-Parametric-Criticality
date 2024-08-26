% Probe_Detection_shPulse202106031140
    load('D:\Kirill\QWJPA_v2_2\03-Jun-2021\Probe_Detection_shPulse202106031254\Data\Probe_Detection_shPulse202106031254_1_1.26e-06_5.8.mat');
% clear all
% load('E:\Kirill\QWJPA_v2_2\03-Jun-2021\Probe_Detection_shPulse202106031152\Data\Probe_Detection_shPulse202106031152_2_2.1e-06_5.8.mat');
    
%     close all
    clear B1
        t_s=(0:nofsamples-1)/samplerate;
        %%
for probeDelayI=1:length(probeDelay)
%         t1=1e-6+0.5e-6+0.203e-6;
%     t2=t1+2e-6;
%     t1=1e-6+0.203e-6;
%     t2=t1+1.5e-6;
%     ClickVarTrue=0;
%     ClickVarFalse=0;
    ClickMeanTrue=0;
    ClickMeanTrueNeg=0;
    ClickMeanFalse=0;
    ClickMeanFalseNeg=0;

for probeState=1:2%TURN OFF-ON probe signal

IQ_raw=zeros(2,nofsamples);

    
%% Make measurement

    load(['D:\Kirill\QWJPA_v2_2\03-Jun-2021\Probe_Detection_shPulse202106031254\Data\Probe_Detection_shPulse202106031254_' num2str(probeState) '_' num2str(probeDelay(probeDelayI)) '_5.8.mat'],'IQ_raw');

    ThresMean=2.1;
    
   
    IQ_cap=(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*tot_phase).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
        ./sqrt(h_p.*(PumpGenFreqCent)/2.*2e6.*50.*(GAIN))).*exp(1i.*1*pi/6);

 
    for cycle_i=1:N_cycles
    B1(:,cycle_i)=IQ_cap(1,round((cycle_i-1)*pulseLength.*samplerate)+(round(t1.*samplerate):round(t2.*samplerate))).';%IQfiltering(IQ_cap(1,round(t1.*samplerate):round(t2.*samplerate)),samplerate,0,1,[-25e6 25e6]).';
    end
    
   if probeState==2 
   LPProbeOnStat(probeDelayI,:)=(mean(B1(:,:),1));
   
    AAA1 = [-1:0.01:1]*max(abs(LPProbeOnStat(probeDelayI,:)));
            [MM1 CC1] = hist3([real(LPProbeOnStat(probeDelayI,:)).' imag(LPProbeOnStat(probeDelayI,:)).'], {AAA1 AAA1});
            hhh=figure(2);
            clf;
            surface(AAA1, AAA1, MM1','edgecolor', 'none');

            axis square
            axis tight
   
   
   
   
    
      for cycle_i=1:N_cycles  
    %     s1=subplot(1,2,2);
    %     cla(s1);
    %     surface(AAA1, AAA1, MM1','edgecolor', 'none');
    %     axis square
    %     axis tight
    %     view(0,90);
    %     title('PROBE ON')
%       HPProbeOnStat(probeDelayI,cycle_i)=var(B1(:,cycle_i));
 
  
      
      
      if real(LPProbeOnStat(probeDelayI,cycle_i))>=ThresMean
          ClickMeanTrue=ClickMeanTrue+1;
      elseif real(LPProbeOnStat(probeDelayI,cycle_i))<=-ThresMean
          ClickMeanTrueNeg=ClickMeanTrueNeg+1;
      end
    end
    elseif probeState==1

                    %     s2=subplot(1,2,1);
                    %     cla(s2);
                    %     surface(AAA1, AAA1, MM1','edgecolor', 'none');
                    %     axis square
                    %     axis tight
                    %     view(0,90);
                    %     title('PROBE OFF')
%      HPProbeOffStat(probeDelayI,cycle_i)=var(B1(:,cycle_i));
     LPProbeOffStat(probeDelayI,:)=(mean(B1(:,:),1));
%      if HPProbeOffStat(probeDelayI,cycle_i)>=ThresVar
%           ClickVarFalse=ClickVarFalse+1;
%      end
    for cycle_i=1:N_cycles 
     if real(LPProbeOffStat(probeDelayI,cycle_i))>=ThresMean
          ClickMeanFalse=ClickMeanFalse+1;
     elseif real(LPProbeOffStat(probeDelayI,cycle_i))<=-ThresMean
          ClickMeanFalseNeg=ClickMeanFalseNeg+1;   
         end
    
    end
    % HPProbeOff(probePowerI)=pks(1);
    % LPProbeOff(probePowerI)=pks(2);
    end
    
%     ClickProbVarTrue(probeDelayI)=ClickVarTrue/N_cycles;
%     ClickProbVarFalse(probeDelayI)=ClickVarFalse/N_cycles;
    ClickProbMeanTrue(probeDelayI)=ClickMeanTrue/N_cycles;
    ClickProbMeanTrueNeg(probeDelayI)=ClickMeanTrueNeg/N_cycles;
    ClickProbMeanFalse(probeDelayI)=ClickMeanFalse/N_cycles;    
    ClickProbMeanFalseNeg(probeDelayI)=ClickMeanFalseNeg/N_cycles; 
    
   
         clearvars A0 L IQ_cap B1
%          save([TempDirPath 'Data\' MeasName, '_' num2str(probeState) '_' num2str(probePower(probePowerI))  '_' num2str(VoltageGeneratorV) '.mat']);
%      save([TempDirPath 'Data\' MeasName,'_' num2str(probeState) '_' num2str(probeDelay(probeDelayI))  '_' num2str(VoltageGeneratorV) '.mat'])
% if (db2pow(probeDelay(1,probeDelayI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)<1.2)&&(db2pow(probeDelay(1,probeDelayI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)>0.8)
%      savefig(hhh,[TempDirPath 'Figures\' MeasName,'_' num2str(probeDelay(probeDelayI)) '_' num2str(VoltageGeneratorV) '.fig'])
% end
end
% mainfig=figure(1144);
% clf;
% subplot(2,2,1);
% semilogx(db2pow(probeDelay(1,1:probeDelayI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1,1:probeDelayI).' HPProbeOff(1,1:probeDelayI).']);
% 
% % semilogx(db2pow(probePower(1,1:probePowerI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1,1:probePowerI).' LPProbeOn(1,1:probePowerI).' HPProbeOff(1,1:probePowerI).' LPProbeOff(1,1:probePowerI).']);
% xlabel('mean N number');
% grid on;
% title('ON-OFF variances');
% % title ('ON-OFF distibutions maxima (+ 10dB attenuation in lines)')
% % figure(1145)
% % clf;
% subplot(2,2,2);
% semilogx(db2pow(probePower(1,1:probeDelayI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [LPProbeOn(1,1:probeDelayI).' LPProbeOff(1,1:probeDelayI).']);
% % semilogx(db2pow(probePower(1,1:probePowerI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
%  title ('ON-OFF means');
%  legend('PROBE ON','PROBE OFF');
% % title ('ON-OFF ratios')
% xlabel('mean N number');
% grid on;
% subplot(2,2,3);
% semilogx(db2pow(probePower(1,1:probeDelayI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1,1:probeDelayI).'-HPProbeOff(1,1:probeDelayI).' LPProbeOn(1,1:probeDelayI).'-LPProbeOff(1,1:probeDelayI).']);
% % semilogx(db2pow(probePower(1,1:probePowerI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
%  title ('ON-OFF diffs of vars(1) and means(2)');
% % title ('ON-OFF ratios')
% xlabel('mean N number');
% grid on;
% subplot(2,2,4);
% semilogx(db2pow(probePower(1,1:probeDelayI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [HPProbeOffSTD(1,1:probeDelayI).' LPProbeOff(1,1:probeDelayI).']);
% % semilogx(db2pow(probePower(1,1:probeDelayI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
%  title ('OFF 2*\sigma of MEANS AND VARS');
% % title ('ON-OFF ratios')
% xlabel('mean N number');
% grid on;
probabfig=figure(1146);
% clf;
% subplot(2,2,3)
plot(probeDelay(1,1:probeDelayI), [ ClickProbMeanTrue(1,1:probeDelayI).' ClickProbMeanTrueNeg(1,1:probeDelayI).'  ClickProbMeanFalse(1,1:probeDelayI).'  ClickProbMeanFalseNeg(1,1:probeDelayI).' ClickProbMeanTrue(1,1:probeDelayI).'+ClickProbMeanTrueNeg(1,1:probeDelayI).'  ClickProbMeanFalse(1,1:probeDelayI).'+ClickProbMeanFalseNeg(1,1:probeDelayI).'],'Linewidth',3);
title ('True and False Positive Rates');
lg=legend('$p_{\mathrm{1+ \lor dark}}^{\rightarrow}$', '$p_{\mathrm{1+ \lor dark}}^{\leftarrow}$', '$p_\mathrm{dark}^{\rightarrow}$', '$p_\mathrm{dark}^{\leftarrow}$','$p_{\mathrm{1+ \lor dark}}$','$p_\mathrm{dark}$');
set(lg,'interpreter','latex');
ylabel('probability');
xlabel('delay');
set(gca,'FontSize',14);
grid on

   
end