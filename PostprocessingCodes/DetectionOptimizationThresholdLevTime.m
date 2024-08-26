
%% CHECK


% Creating ZERO arrays for timetraced data A1  AAA B1 L1 B11 A0 L t_s dem IQ_raw g G h H
n_modes=2;
% IQ_raw = zeros(2,nofsamples);
% IQ_cap = zeros(1,nofsamples);
% N_cycles=1000;
% pulseLength=10e-6;
% probeLength=2e-6;%tau
    t1=0;
    t2=pulseLength-1e-8;

probeDelay=linspace(0e-6,4.99e-6,501);
ThresVal=linspace(0.1,20,101);

ClickTrue=zeros(length(probeDelay),length(ThresVal));
ClickFalse=zeros(length(probeDelay),length(ThresVal));
Eff=zeros(length(probeDelay),length(ThresVal));
for probeState=1:2%TURN OFF-ON probe signal
  if probeState==1
      clear IQ_cap B1
      load('E:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\Data\Probe_Detection_shPulse_PhotonNumberSweep202106080307_1_-42.5654_5.78.mat','IQ_raw')
  elseif probeState==2   
      clear IQ_cap B1
       load('E:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\Data\Probe_Detection_shPulse_PhotonNumberSweep202106080307_2_-42.5654_5.78.mat','IQ_raw')
  end    
       t_s=(0:nofsamples-1)/samplerate;
      IQ_cap=(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*tot_phase).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
            ./sqrt(h_p.*(PumpGenFreqCent)/2.*2e6.*50.*(GAIN)));

  
      for cycle_i=1:N_cycles
        B1(:,cycle_i)=abs(IQ_cap(1,1+round((cycle_i-1)*pulseLength.*samplerate)+(round(t1.*samplerate):round(t2.*samplerate)))).';%IQfiltering(IQ_cap(1,round(t1.*samplerate):round(t2.*samplerate)),samplerate,0,1,[-25e6 25e6]).';
      end

    for probeDelayI=1:length(probeDelay)


    for ThresVal_i=1:length(ThresVal)
        ClickValTrue=0;
            ClickValFalse=0;
        for cycle_i=1:N_cycles
            
         
         if probeState==1

         if B1(1+round(samplerate*probeDelay(probeDelayI)),cycle_i)>=ThresVal(ThresVal_i)
              ClickValFalse=ClickValFalse+1;
         end
         
            
         elseif probeState==2

          if ((B1(1+round(samplerate*probeDelay(probeDelayI)),cycle_i)))>=ThresVal(ThresVal_i)
              ClickValTrue=ClickValTrue+1;
          end
        end

        end
         
         if probeState==1
             ClickFalse(probeDelayI,ThresVal_i)=ClickValFalse/N_cycles;
             
         elseif probeState==2
             ClickTrue(probeDelayI,ThresVal_i)=ClickValTrue/N_cycles;
         end
         Eff(probeDelayI,ThresVal_i)=ClickTrue(probeDelayI,ThresVal_i).*(1-ClickFalse(probeDelayI,ThresVal_i));
        end

   
    
    end
 
%      save([TempDirPath 'Data\' MeasName, '_' num2str(probePower(probeDelayI))  '_' num2str(VoltageGeneratorV) '.mat'], '-regexp', '^(?!(IQ|A1|IQ_init|IQ_cap|AAA|B1|L1|B11|A0|L|t_s|dem|IQ_raw|g|G|h|H)$).')

end
 figure(667);
ss=surf(ThresVal,probeDelay,ClickFalse);set(ss,'EdgeColor','none')
figure(665);
ss2=surf(ThresVal,probeDelay,Eff);set(ss2,'EdgeColor','none')

figure(668);
ss3=surf(ThresVal,probeDelay,ClickTrue);set(ss3,'EdgeColor','none')
