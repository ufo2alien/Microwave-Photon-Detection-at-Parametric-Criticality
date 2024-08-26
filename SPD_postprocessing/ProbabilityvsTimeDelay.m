

%% CHECK


% Creating ZERO arrays for timetraced data A1  AAA B1 L1 B11 A0 L t_s dem IQ_raw g G h H
n_modes=2;
% IQ_raw = zeros(2,nofsamples);
% IQ_cap = zeros(1,nofsamples);
% N_cycles=1000;
% pulseLength=10e-6;
% probeLength=2e-6;%tau
probeDelay=0.203e-6+linspace(-0.2e-6,2e-6,21);
t1=0.5e-6+probeDelay;
t2=0.5e-6+0.5e-6+probeLength+probeDelay;


for probeDelayI=1:length(probeDelay)
    ClickVarTrue=0;
    ClickVarFalse=0;
    ClickMeanTrue=0;
    ClickMeanFalse=0;
for probeState=2%TURN OFF-ON probe signal
    
%% Make measurement
t_s=(0:nofsamples-1)/samplerate;
%     IQ_cap=(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*tot_phase).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
%         ./sqrt(h_p.*(PumpGenFreqCent)/2.*2e6.*50.*(GAIN)));
  

    t1=0.5e-6+probeDelay(probeDelayI);
    t2=0.5e-6+0.5e-6+probeLength+probeDelay(probeDelayI);

    
    for cycle_i=1:N_cycles
    B1(:,cycle_i)=IQ_cap(1,round((cycle_i-1)*pulseLength.*samplerate)+(round(t1.*samplerate):round(t2.*samplerate))).';%IQfiltering(IQ_cap(1,round(t1.*samplerate):round(t2.*samplerate)),samplerate,0,1,[-25e6 25e6]).';
    
    if probeState==2
    %     s1=subplot(1,2,2);
    %     cla(s1);
    %     surface(AAA1, AAA1, MM1','edgecolor', 'none');
    %     axis square
    %     axis tight
    %     view(0,90);
    %     title('PROBE ON')
      HPProbeOnStat(probeDelayI,cycle_i)=var(B1(:,cycle_i));
      LPProbeOnStat(probeDelayI,cycle_i)=abs(mean(B1(:,cycle_i)));
      if HPProbeOnStat(probeDelayI,cycle_i)>=ThresVar(25)
          ClickVarTrue=ClickVarTrue+1;
      end
      
      
      if LPProbeOnStat(probeDelayI,cycle_i)>=ThresMean(55)
          ClickMeanTrue=ClickMeanTrue+1;
      end
      
    elseif probeState==1

                    %     s2=subplot(1,2,1);
                    %     cla(s2);
                    %     surface(AAA1, AAA1, MM1','edgecolor', 'none');
                    %     axis square
                    %     axis tight
                    %     view(0,90);
                    %     title('PROBE OFF')
     HPProbeOffStat(probeDelayI,cycle_i)=var(B1(:,cycle_i));
     LPProbeOffStat(probeDelayI,cycle_i)=abs(mean(B1(:,cycle_i)));
     if HPProbeOffStat(probeDelayI,cycle_i)>=ThresVar
          ClickVarFalse=ClickVarFalse+1;
     end
     
     if LPProbeOffStat(probeDelayI,cycle_i)>=ThresMean
          ClickMeanFalse=ClickMeanFalse+1;
     end
    % HPProbeOff(probeDelayI)=pks(1);
    % LPProbeOff(probeDelayI)=pks(2);
    end
    end
    ClickProbVarTrue(probeDelayI)=ClickVarTrue/N_cycles;
    ClickProbVarFalse(probeDelayI)=ClickVarFalse/N_cycles;
    ClickProbMeanTrue(probeDelayI)=ClickMeanTrue/N_cycles;
    ClickProbMeanFalse(probeDelayI)=ClickMeanFalse/N_cycles;    
    
    if probeState==2
        HPProbeOn(1,probeDelayI)=mean(HPProbeOnStat(probeDelayI,:),2);
        LPProbeOn(1,probeDelayI)=mean(LPProbeOnStat(probeDelayI,:),2);
%         HPProbeOnSTD(1,probeDelayI)=2*std(HPProbeOnStat(probeDelayI,:),2);
%         LPProbeOnSTD(1,probeDelayI)=2*std(LPProbeOnStat(probeDelayI,:),2);
    elseif probeState==1
        LPProbeOff(1,probeDelayI)=mean(LPProbeOffStat(probeDelayI,:),2);
        HPProbeOff(1,probeDelayI)=mean(HPProbeOffStat(probeDelayI,:),2);
        HPProbeOffSTD(1,probeDelayI)=2*std(HPProbeOffStat(probeDelayI,:));
        LPProbeOffSTD(1,probeDelayI)=2*std(LPProbeOffStat(probeDelayI,:));
    end

%      save([TempDirPath 'Data\' MeasName, '_' num2str(probePower(probeDelayI))  '_' num2str(VoltageGeneratorV) '.mat'], '-regexp', '^(?!(IQ|A1|IQ_init|IQ_cap|AAA|B1|L1|B11|A0|L|t_s|dem|IQ_raw|g|G|h|H)$).')
if (db2pow(probeDelay(1,probeDelayI)-101.4-30)*Energy./(h_p.*PumpGenFreqCent/2)<1.2)&&(db2pow(probeDelay(1,probeDelayI)-101.4-30)*Energy./(h_p.*PumpGenFreqCent/2)>0.8)
%      savefig(hhh,[TempDirPath 'Figures\' MeasName,'_' num2str(probePower(probeDelayI)) '_' num2str(VoltageGeneratorV) '.fig'])
end
end

probabfig=figure(1146);
% clf;
% subplot(2,2,3)
plot(probeDelay(1,1:probeDelayI), [ ClickProbVarTrue(1,1:probeDelayI).'  ClickProbVarFalse(1,1:probeDelayI).' ClickProbMeanTrue(1,1:probeDelayI).'  ClickProbMeanFalse(1,1:probeDelayI).'],'Linewidth',3);
title ('True and False Positive Rates');
legend('Variance TPR','Variance FPR', 'Mean TPR', 'Mean FPR');
ylabel('probability');
xlabel('delay');
set(gca,'FontSize',14);
grid on
    
end




    
    
    