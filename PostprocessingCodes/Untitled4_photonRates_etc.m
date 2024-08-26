

%% CHECK

% load('D:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\PoissonFitDataPulseCharact.mat');

% Creating ZERO arrays for timetraced data A1  AAA B1 L1 B11 A0 L t_s dem IQ_raw g G h H
% n_modes=2;
% IQ_raw = zeros(2,nofsamples);
% IQ_cap = zeros(1,nofsamples);
% N_cycles=1000;
% pulseLength=10e-6;
% probeLength=2e-6;%tau
% probeDelay=0.203e-6+linspace(-0.2e-6,2e-6,21);
% t1=0.5e-6+probeDelay;
% t2=0.5e-6+0.5e-6+probeLength+probeDelay;


for probePowerI=1:length(probePower)
%     load('D:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\Data\Probe_Detection_shPulse_PhotonNumberSweep202106080307_2_-48.71_5.78.mat','IQ_raw','probeState')
    ClickMeanTrue=0;
    ClickMeanFalse=0;   
for probeState=1:2%TURN OFF-ON probe signal
    load(['D:\Kirill\QWJPA_v2_2\09-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106091948\Data\Probe_Detection_shPulse_PhotonNumberSweep202106091948_' num2str(probeState) '_' num2str(probePower(probePowerI)) '_5.78.mat'],'IQ_raw');
%% Make measurement
t_s=(0:nofsamples-1)/samplerate;
    IQ_cap=(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*tot_phase).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
        ./sqrt(h_p.*(PumpGenFreqCent)/2.*2e6.*50.*(GAIN)));
  
ThresMean=-2.75;
  
    for cycle_i=1:N_cycles
    B1(:,cycle_i)=IQ_cap(1,round((cycle_i-1)*pulseLength.*samplerate)+(round(t1.*samplerate):round(t2.*samplerate))).';%IQfiltering(IQ_cap(1,round(t1.*samplerate):round(t2.*samplerate)),samplerate,0,1,[-25e6 25e6]).';
    
    if probeState==2

      LPProbeOnStat(probePowerI,cycle_i)=(mean(B1(:,cycle_i)));

      
      
      if imag(LPProbeOnStat(probePowerI,cycle_i))<=ThresMean
          ClickMeanTrue=ClickMeanTrue+1;
      end
      
    elseif probeState==1


     LPProbeOffStat(probePowerI,cycle_i)=(mean(B1(:,cycle_i)));

     
     if imag(LPProbeOffStat(probePowerI,cycle_i))<=ThresMean
          ClickMeanFalse=ClickMeanFalse+1;
     end

    end
    end
    
    if probeState==2
        ClickProbMeanTrue(probePowerI)=ClickMeanTrue/N_cycles;
    elseif probeState==1
        ClickProbMeanFalse(probePowerI)=ClickMeanFalse/N_cycles;    
    end
    
    

 hhh = figure (28);
    if probeState==1
        subplot(1,2,1);
        cla;
         AA = [-1:0.0025:1].*max(abs(LPProbeOffStat(probePowerI,:)));
  
    %     AA = [-1:0.01:1].*max(abs(b_in));
    [N1 C1] =hist3([real(LPProbeOffStat(probePowerI,:)).'  imag(LPProbeOffStat(probePowerI,:)).'], {AA AA});
      surf(AA, AA, N1./length(LPProbeOffStat(probePowerI,:)), 'edgecolor', 'none');
        axis tight;
        daspect([1 1 1])
    view(2);
  
    else
       subplot(1,2,2);
        cla;
         AA2 = [-1:0.0025:1].*max(abs(LPProbeOnStat(probePowerI,:)));
          [N2 C2] =hist3([real(LPProbeOnStat(probePowerI,:)).', imag(LPProbeOnStat(probePowerI,:)).'], {AA2 AA2}); 
         surf(AA2, AA2, N2./length(LPProbeOnStat(probePowerI,:)), 'edgecolor', 'none');
           axis tight;
           daspect([1 1 1])
    view(2);
    end
end

probabfig=figure(1145);
% clf;
% subplot(2,2,3)
semilogx(db2pow(probePower(1,1:probePowerI)-101.3-30 + 4.25)*Energy./(h_p.*PumpGenFreqCent/2), [ ClickProbMeanTrue(1,1:probePowerI).'  ClickProbMeanFalse(1,1:probePowerI).' (ClickProbMeanTrue(1,1:probePowerI).'-ClickProbMeanFalse(1,1:probePowerI).')./(1-ClickProbMeanFalse(1,1:probePowerI).')],'Linewidth',3);
legend('Mean TPR', 'Mean FPR', 'p1+');
ylabel('Probability','interpreter','latex')
xlabel('Mean photon number in pulse, $\left< n \right>$','interpreter','latex')
set(gca,'FontSize',36,'FontName','Segoe UI Semibold');
grid on

end




    
    
    