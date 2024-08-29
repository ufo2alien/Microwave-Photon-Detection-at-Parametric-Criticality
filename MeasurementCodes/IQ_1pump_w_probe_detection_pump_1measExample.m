
clearvars
close all
clear mex    
%%
DeviceName = 'QWJPA_v2_2';
DirPath='C:\MatlabPrograms\PetrovninK\';
cd('C:\MatlabPrograms\PetrovninK\');
addpath('C:\Matlabprograms\PetrovninK\functions_for_calc');

filename=todayfilename;
MeasName = ['measOfOscillation',filename];
TempDirPath=['G:\Kirill\' DeviceName '\' datestr(date) '\' MeasName '\'];
addpath(TempDirPath);
mkdir(TempDirPath);
mkdir([TempDirPath 'Data\']);
mkdir([TempDirPath 'Figures\']);
addpath(TempDirPath);
CurrentFileName=mfilename('fullpath');

copyfile([CurrentFileName,'.m'],[TempDirPath,MeasName,'.m']);

%  pumpState=1;
 
%% Set parameters and initalize device
atten=0;
NAdet=-0.00e6;
FIXpump=0.65;linspace(0.40,0.60,81);0.67;0.365;0.76;%[linspace(-41,-6,36) linspace(-5,7,21) linspace(7.2,13,25) linspace(13,20.5,11)]-atten;%[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
probePowerRaw=0.025;0.0668;0.0477;0.0347;logspace(-2.5,-1,21);%linspace(-40,10,41);
probePower=probeVppCalibration(probePowerRaw);%probeCalibration(probePowerRaw);
% pump = [linspace(0.5,0.8,11)]-atten;%linspace(-20,-2,10) [linspace(-41,-6,18) linspace(-5,7,21) linspace(7.2,16,25) linspace(16.5,20.5,9)]-atten;%[linspace(-20,-2,10) linspace(0,15,16)]-atten;%%linspace(-41,-6,18)  linspace(-2,18,21)-atten;[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
PumpGenFreqCent = 12.084e9;%12.0778e9;%11.995e9;
MeasFreq = 14e6;% 6052499930.9;%PumpGenFreqCent/2-6.052533066e9+PumpGenFreqCent/2;%PumpGenFreqCent/2   
RBW = 2e6; 
SGCh2 = 3;
SGCh3 = 1; %READOUT LO oscillator channel
SGCh4 = 2;%PROBE channel
n_ave = 1;
GAIN = db2pow(105-8.2); %-16.42Quadrature gain 90.7 +- 0.2dB%Quadrature noise temperature 5.11K coth(hw/2kT)=35.9
h_p = 6.62607004e-34;
VoltageGeneratorV = 5.86;0.741;0.729;0.729;%0.7392;
samplerate = 100e6;%(2^7);
kb = 1.380649e-23;
tot_phase=pi/4;

N_cycles=100;
pulseLength=20e-5;
pumpLength=18e-5;
probeLength=2e-6;%tau
probeDelay=1e-6+0.223e-6;linspace(0e-6,2.5e-6,335);0.206e-6;
probeDetuning=0*-0.13e6;%linspace(-3e6,3e6,601);


%Threshold

% ThresMean=12.94/2;3.918;6.49;4;
% ThresVar=22.23;12.94;16.9;14.5;

% load('D:\Kirill\QWJPA_v2_1\20-Apr-2021\IQ_1pump_Digitizer_w_Triggering202104200138\Data\IQ_1pump_Digitizer_w_Triggering202104200138pumpState_1_1_0.7565.mat','Noise_Gain');
%% CONFIGURE ADQ14 digitizer
boardid = 1;
nofsamples = round(pulseLength*samplerate*N_cycles);
nof_records = 1;
trig_delay=0e-3;
channels_mask = 3; %Channels 1&2
trigger = 2;%0 - SW trigger, 1 - Internal, 2 - External and 3 - Lvl trigger
MeasLength=nofsamples/samplerate;

[samplerate]=ADQ14_init(boardid,nofsamples,nof_records,channels_mask,1e9/samplerate,trigger,trig_delay);%Init and Read precise value of samplerate

%% CONFIGURE TRIGGER GENERATOR
TriggerGenerator = visa('agilent','TCPIP0::169.254.58.12::inst0::INSTR','OutputBufferSize',10^8);
old_obj = instrfind('type','visa-tcpip','RemoteHost','169.254.58.12');
if(~isempty(old_obj))
    fclose(old_obj);
end
fopen(TriggerGenerator);


[ExpPulse,Energy,ExpPulseStr] = ArbExpPulse(1,1000*probeLength,250e6,250e6*pulseLength+1,1,250e6*probeLength+1,MeasFreq+probeDetuning(1));
[ExpPulse2,Energy2,ExpPulseStr2] = ArbExpPulse(1,20000e-6,250e6,250e6*pulseLength+1,250e6*1e-6+1,250e6*(pumpLength+1e-6)+1,MeasFreq*2);
fprintf(TriggerGenerator,'ROSC:SOUR EXT');
fprintf(TriggerGenerator,'SOUR1:VOLT:UNIT VPP');
fprintf(TriggerGenerator,'SOUR2:VOLT:UNIT VPP');
fprintf(TriggerGenerator,['SOUR1:APPL:ARB']);
fprintf(TriggerGenerator,['SOUR2:APPL:ARB']);
fprintf(TriggerGenerator,['SOUR1:DATA:VOL:CLE']);
fprintf(TriggerGenerator,['SOUR2:DATA:VOL:CLE']);
% fprintf(TriggerGenerator,['SOUR1:APPL:SIN ' num2str(MeasFreq*2) ',' num2str(FIXpump) ',0']);

% fprintf(TriggerGenerator,['SOUR2:APPL:SIN ' num2str(MeasFreq)-NAdet ',-50,0.0066']);
fprintf(TriggerGenerator,['SOUR1:DATA:ARB EXPDEC2, ' ExpPulseStr2]);
fprintf(TriggerGenerator,['SOUR2:DATA:ARB EXPDEC, ' ExpPulseStr]);
pause(0.1)
fprintf(TriggerGenerator,['SOUR2:VOLT ' num2str(0.005)]);
fprintf(TriggerGenerator,['SOUR2:VOLT:OFFS 0.0066']);
fprintf(TriggerGenerator,['SOUR1:FUNC:ARB:SRAT ' num2str(250e6)]);
fprintf(TriggerGenerator,['SOUR1:FUNC:ARB EXPDEC2']);
fprintf(TriggerGenerator,['SOUR2:FUNC:ARB:SRAT ' num2str(250e6)]);
fprintf(TriggerGenerator,['SOUR2:FUNC:ARB EXPDEC']);
fprintf(TriggerGenerator,['SOUR1:VOLT ' num2str(FIXpump(1))]);
fprintf(TriggerGenerator,'OUTP1 ON');
fprintf(TriggerGenerator,'OUTP2 ON');
fprintf(TriggerGenerator,['OUTP1:LOAD 50']);
fprintf(TriggerGenerator,['OUTP2:LOAD 50']);
fprintf(TriggerGenerator,'SOUR1:BURS:MOD TRIG');
fprintf(TriggerGenerator,['SOUR1:BURS:NCYC ' num2str(N_cycles) ]);
% fprintf(TriggerGenerator,['SOUR1:BURS:NCYC ' num2str(pulseLength*MeasFreq*2) ]);
fprintf(TriggerGenerator,'SOUR2:BURS:MOD TRIG');
fprintf(TriggerGenerator,['SOUR2:BURS:NCYC ' num2str(N_cycles) ]);
% fprintf(TriggerGenerator,['SOUR2:BURS:NCYC ' num2str(probeLength*MeasFreq) ]);
fprintf(TriggerGenerator,'UNIT:ANGLE DEG');
fprintf(TriggerGenerator,'SOUR2:BURS:PHAS 0');
fprintf(TriggerGenerator,'OUTP:TRIG 1')
fprintf(TriggerGenerator,'OUTP:SYNC:SOUR CH1')
fprintf(TriggerGenerator,'OUTP:SYNC:MODE NORM')
fprintf(TriggerGenerator,'OUTP:SYNC:SOUR CH1')
% fprintf(TriggerGenerator,['SOUR1:MARK:CYCL ' num2str(1)]);%MeasFreq*2*t1
fprintf(TriggerGenerator,'TRIG1:SOUR BUS');
fprintf(TriggerGenerator,'TRIG1:DEL 0e-6');
fprintf(TriggerGenerator,['TRIG2:DEL ' num2str(probeDelay(1))]);

fprintf(TriggerGenerator,'SOUR1:BURS:STAT ON');
fprintf(TriggerGenerator,'TRIG2:SOUR BUS');
fprintf(TriggerGenerator,'SOUR2:BURS:STAT ON');

%% CONFIGURE DC generator
VoltageGenerator = visa('agilent','TCPIP0::169.254.58.15::inst0::INSTR');
old_obj = instrfind('type','visa-tcpip','RemoteHost','169.254.58.15');
        if(~isempty(old_obj))
        fclose(old_obj);
        end
fopen(VoltageGenerator);


DCCh = 1; 
DCCh2 = 2;
% fprintf(VoltageGenerator,'ROSC:SOUR:EXT');
fprintf(VoltageGenerator,['SOUR' num2str(DCCh) ':APPL:DC']);
fprintf(VoltageGenerator,['SOUR' num2str(DCCh2) ':APPL:DC']);
% fprintf(VoltageGenerator,['SOUR' num2str(DCCh2) ':APPL:SIN ' num2str(20e6) ',0.005,0']);
fprintf(VoltageGenerator,['OUTP' num2str(DCCh) ':LOAD INF']);
fprintf(VoltageGenerator,['OUTP' num2str(DCCh2) ':LOAD INF']);
fprintf(VoltageGenerator,['SOUR' num2str(DCCh) ':VOLT:OFFS ' num2str(VoltageGeneratorV)]);
fprintf(VoltageGenerator,['SOUR' num2str(DCCh2) ':VOLT:OFFS ' num2str(-4)]);
fprintf(VoltageGenerator,['OUTP' num2str(DCCh) ' ON']);%REMEMBER TO TURN ON!
fprintf(VoltageGenerator,['OUTP' num2str(DCCh2) ' ON']);

%% CONFIGURE SIGNAL GENERATOR:PREPARATION
SignalGenerator = tcpip('169.254.58.13',18); % Anapico
% SignalGenerator = visa('ni','TCPIP0::130.233.175.193::inst0::INSTR');%visa('agilent','GPIB0::9::INSTR');
old_obj = instrfind('type','tcpip','RemoteHost','169.254.58.13');
if(~isempty(old_obj))
        fclose(old_obj);
end
fopen(SignalGenerator);
% fprintf(SignalGenerator,'*RST');
% fprintf(SignalGenerator, ['ROSC:SOUR EXT']);
% fprintf(SignalGenerator, ['SOURce:LFO:STAT 0']);
fileNameCh2='Pump';
fileNameCh3='LOFreq';
fileNameCh4='Probe';
deleteListAnapico(SignalGenerator,fileNameCh2);
deleteListAnapico(SignalGenerator,fileNameCh3);
% deleteListAnapico(SignalGenerator,fileNameCh4);
createListAnapico(SignalGenerator,fileNameCh3,SGCh3,PumpGenFreqCent/2-MeasFreq,18,MeasLength,0);
createListAnapico(SignalGenerator,fileNameCh2,SGCh2,PumpGenFreqCent-2*MeasFreq,15,MeasLength,0);
% createListAnapico(SignalGenerator,fileNameCh4,SGCh4,PumpGenFreqCent/2,-70,MeasLength,0);
loadListAnapico(SignalGenerator,fileNameCh2,SGCh2);
loadListAnapico(SignalGenerator,fileNameCh3,SGCh3);
% loadListAnapico(SignalGenerator,fileNameCh4,SGCh4);
disp(query(SignalGenerator,['MEM' num2str(SGCh2) ':FILE:LIST:DATA?']));
disp(query(SignalGenerator,['MEM' num2str(SGCh3) ':FILE:LIST:DATA?']));
fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':FREQ:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':POW:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':FREQ:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':POW:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh4) ':FREQ:MODE CW']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh4) ':POW:MODE CW']);
fprintf(SignalGenerator, ['TRIG:TYPE NORM']);
fprintf(SignalGenerator, ['INIT:CONT 0']);
fprintf(SignalGenerator, ['TRIG:DEL 0e-3']);
fprintf(SignalGenerator, ['TRIG:SOUR EXT']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':LIST:COUNT 1']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':LIST:COUNT 1']);
% fprintf(SignalGenerator, ['SOUR' num2str(SGCh4) ':LIST:COUNT 1']);
fprintf(SignalGenerator, ['OUTP' num2str(SGCh2) ':STAT ', num2str(1)]); 
fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(1)]); 
fprintf(SignalGenerator, ['OUTP' num2str(SGCh4) ':STAT ', num2str(0)]); 
pause(0.5)

%% CHECK

tot_meas_all=length(FIXpump).*(n_ave);
tot_meas_done=0;
t_ar=[];
tic

% Creating ZERO arrays for timetraced data A1  AAA B1 L1 B11 A0 L t_s dem IQ_raw g G h H
n_modes=2;
IQ_raw = zeros(2,nofsamples);
IQ_cap = zeros(1,nofsamples);
% B1 = zeros(5e2+1,1);
% B11  = zeros(nofsamples,2);
% L1 = zeros(nofsamples,4);


for pumpI=1:length(FIXpump)
    
%     [ExpPulse,Energy,ExpPulseStr] = ArbExpPulse(1,1000*probeLength,250e6,250e6*pulseLength+1,1,250e6*probeLength+1,MeasFreq+probeDetuning(1));
%     fprintf(TriggerGenerator,['SOUR2:APPL:ARB']);
%     fprintf(TriggerGenerator,['SOUR2:DATA:VOL:CLE']);
%     pause(0.2);
%     fprintf(TriggerGenerator,['SOUR2:DATA:ARB EXPDEC, ' ExpPulseStr]);
%     pause(0.3)
%     fprintf(TriggerGenerator,['SOUR2:FUNC:ARB:SRAT ' num2str(250e6)]);
%     fprintf(TriggerGenerator,['SOUR2:FUNC:ARB EXPDEC']);
%     pause(0.2)
%     fprintf(TriggerGenerator,['SOUR2:VOLT:OFFS 0.0066']);
%     fprintf(TriggerGenerator,['TRIG2:DEL ' num2str(probeDelay(1))]);
%     fprintf(TriggerGenerator,'SOUR2:BURS:MOD TRIG');
%     fprintf(TriggerGenerator,['SOUR2:BURS:NCYC ' num2str(N_cycles) ]);
%     fprintf(TriggerGenerator,'TRIG2:SOUR BUS');
%     fprintf(TriggerGenerator,'SOUR2:BURS:STAT ON'); 
    fprintf(TriggerGenerator,['SOUR1:VOLT ' num2str(FIXpump(pumpI))]);
    pause(0.2)
    t1=probeDelay;
    t2=t1+probeLength+16e-6;
    ClickVarTrue=0;
    ClickVarFalse=0;
    ClickMeanTrue=0;
    ClickMeanFalse=0;
for probeState=1:1%TURN OFF-ON probe signal

    if probeState==2
fprintf(TriggerGenerator,['SOUR2:VOLT ' num2str(probePowerRaw(1))]);

fprintf(TriggerGenerator,'OUTP2 ON');
    elseif probeState==1
fprintf(TriggerGenerator,['SOUR2:VOLT ' num2str(0.001)]);
fprintf(TriggerGenerator,'OUTP2 OFF');

    end
pause(0.5)
% fprintf(NetworkAnalyzer, ['SOUR1:POW ' num2str(probePower(probePowerI))]);
% fprintf(NetworkAnalyzer, ['OUTP1 ' num2str((ii-1))]);
% pause(0.5)
IQ_raw=zeros(2,nofsamples);
for n_a=1:n_ave
    
%% Make measurement

    fprintf(SignalGenerator,'INIT');
    pause(0.5)
    IQ_raw = IQ_raw+(ADQ14_acq(boardid,nofsamples,nof_records,channels_mask,trigger,2,TriggerGenerator));
end
IQ_raw=IQ_raw./n_ave;
    t_s=(0:nofsamples-1)/samplerate;
   
    IQ_cap=(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*tot_phase).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
        ./sqrt(h_p.*(PumpGenFreqCent)/2.*2e6.*50.*(GAIN)));
  
    hhh = figure (28);
    if probeState==1
        subplot(1,2,1);
        cla;
        plot(t_s,[real(IQ_cap(1,:)).' imag(IQ_cap(1,:)).']);
        hold on
        plot([t1 t1],ylim,'--g');
         plot([t2 t2],ylim,'--g');
         hold off
        title('PROBE OFF');
    else
       subplot(1,2,2);
        cla;
        plot(t_s,[real(IQ_cap(1,:)).' imag(IQ_cap(1,:)).']);
        hold on
        plot([t1 t1],ylim,'--g');
        plot([t2 t2],ylim,'--g');
        hold off
        title('PROBE ON');
    end
    
    
%     for cycle_i=1:N_cycles
%     B1(:,cycle_i)=IQ_cap(1,round((cycle_i-1)*pulseLength.*samplerate)+(round(t1.*samplerate):round(t2.*samplerate))).';%IQfiltering(IQ_cap(1,round(t1.*samplerate):round(t2.*samplerate)),samplerate,0,1,[-25e6 25e6]).';
    
%     if probeState==2
    %     s1=subplot(1,2,2);
    %     cla(s1);
    %     surface(AAA1, AAA1, MM1','edgecolor', 'none');
    %     axis square
    %     axis tight
    %     view(0,90);
    %     title('PROBE ON')
%       HPProbeOnStat(pumpI,cycle_i)=var(B1(:,cycle_i));
%       LPProbeOnStat(pumpI,cycle_i)=abs(mean(B1(:,cycle_i)));
%       if HPProbeOnStat(pumpI,cycle_i)>=ThresVar
%           ClickVarTrue=ClickVarTrue+1;
%       end
      
      
%       if LPProbeOnStat(pumpI,cycle_i)>=ThresMean
%           ClickMeanTrue=ClickMeanTrue+1;
%       end
      
%     elseif probeState==1
% 
%                     %     s2=subplot(1,2,1);
%                     %     cla(s2);
%                     %     surface(AAA1, AAA1, MM1','edgecolor', 'none');
%                     %     axis square
%                     %     axis tight
%                     %     view(0,90);
%                     %     title('PROBE OFF')
%      HPProbeOffStat(pumpI,cycle_i)=var(B1(:,cycle_i));
%      LPProbeOffStat(pumpI,cycle_i)=abs(mean(B1(:,cycle_i)));
%      if HPProbeOffStat(pumpI,cycle_i)>=ThresVar
%           ClickVarFalse=ClickVarFalse+1;
%      end
%      
%      if LPProbeOffStat(pumpI,cycle_i)>=ThresMean
%           ClickMeanFalse=ClickMeanFalse+1;
%      end
%     % HPProbeOff(probePowerI)=pks(1);
%     % LPProbeOff(probePowerI)=pks(2);
%     end
%     end
%     ClickProbVarTrue(pumpI)=ClickVarTrue/N_cycles;
%     ClickProbVarFalse(pumpI)=ClickVarFalse/N_cycles;
%     ClickProbMeanTrue(pumpI)=ClickMeanTrue/N_cycles;
%     ClickProbMeanFalse(pumpI)=ClickMeanFalse/N_cycles;    
    
%     if probeState==2
%         HPProbeOn(1,pumpI)=mean(HPProbeOnStat(pumpI,:),2);
%         LPProbeOn(1,pumpI)=mean(LPProbeOnStat(pumpI,:),2);
% %         HPProbeOnSTD(1,probeDetuningI)=2*std(HPProbeOnStat(probePowerI,:),2);
% %         LPProbeOnSTD(1,probePowerI)=2*std(LPProbeOnStat(probePowerI,:),2);
%     elseif probeState==1
%         LPProbeOff(1,pumpI)=mean(LPProbeOffStat(pumpI,:),2);
%         HPProbeOff(1,pumpI)=mean(HPProbeOffStat(pumpI,:),2);
%         HPProbeOffSTD(1,pumpI)=2*std(HPProbeOffStat(pumpI,:));
%         LPProbeOffSTD(1,pumpI)=2*std(LPProbeOffStat(pumpI,:));
%     end
%          clearvars A0 L IQ_cap t_s B1
%          save([TempDirPath 'Data\' MeasName, '_' num2str(probeState) '_' num2str(probePower(probePowerI))  '_' num2str(VoltageGeneratorV) '.mat']);
%      save([TempDirPath 'Data\' MeasName,'_' num2str(probeState) '_' num2str(probeDetuning(probeDetuningI))  '_' num2str(VoltageGeneratorV) '.mat'])
% if (db2pow(probeDelay(1,probeDetuningI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)<1.2)&&(db2pow(probeDelay(1,probeDetuningI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2)>0.8)
%      savefig(hhh,[TempDirPath 'Figures\' MeasName,'_' num2str(probeDelay(probeDetuningI)) '_' num2str(VoltageGeneratorV) '.fig'])
    
%     end
%     if probeState==1
%         B1_ave_OFF=mean(B1,2);
%     else
%         B1_ave_ON=mean(B1,2);
%     end

% figure (676)
%  if probeState==1
%         subplot(1,2,1);
% %         cla;
%         hold on
%         plot(linspace(t1,t2,size(B1_ave_OFF,1)),abs(B1_ave_OFF));
% %         hold on
% %          hold off
%         title('PROBE OFF');
%     else
%        subplot(1,2,2);
%        hold on
%        plot(linspace(t1,t2,size(B1_ave_ON,1)),abs(B1_ave_ON));
% %         cla;
% %         plot(t_s,[real(IQ_cap(1,:)).' imag(IQ_cap(1,:)).']);
% %         hold on
% %         plot([t1 t1],ylim,'--g');
% %         plot([t2 t2],ylim,'--g');
% %         hold off
%         title('PROBE ON');
%  end
  if probeState==1
  save([TempDirPath 'Data\' MeasName,'_' num2str(probeState) '_' num2str(FIXpump(pumpI))  '_' num2str(VoltageGeneratorV) '.mat'])
  end
end
% plot(linspace(t1,t2,size(B1_ave,1)),abs(B1_ave));
% mainfig=figure(1144);
% clf;
% subplot(2,2,1);
% semilogx(db2pow(probeDelay(1,1:probeDetuningI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1,1:probeDetuningI).' HPProbeOff(1,1:probeDetuningI).']);
% 
% % semilogx(db2pow(probePower(1,1:probePowerI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1,1:probePowerI).' LPProbeOn(1,1:probePowerI).' HPProbeOff(1,1:probePowerI).' LPProbeOff(1,1:probePowerI).']);
% xlabel('mean N number');
% grid on;
% title('ON-OFF variances');
% % title ('ON-OFF distibutions maxima (+ 10dB attenuation in lines)')
% % figure(1145)
% % clf;
% subplot(2,2,2);
% semilogx(db2pow(probePower(1,1:probeDetuningI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [LPProbeOn(1,1:probeDetuningI).' LPProbeOff(1,1:probeDetuningI).']);
% % semilogx(db2pow(probePower(1,1:probePowerI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
%  title ('ON-OFF means');
%  legend('PROBE ON','PROBE OFF');
% % title ('ON-OFF ratios')
% xlabel('mean N number');
% grid on;
% subplot(2,2,3);
% semilogx(db2pow(probePower(1,1:probeDetuningI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1,1:probeDetuningI).'-HPProbeOff(1,1:probeDetuningI).' LPProbeOn(1,1:probeDetuningI).'-LPProbeOff(1,1:probeDetuningI).']);
% % semilogx(db2pow(probePower(1,1:probePowerI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
%  title ('ON-OFF diffs of vars(1) and means(2)');
% % title ('ON-OFF ratios')
% xlabel('mean N number');
% grid on;
% subplot(2,2,4);
% semilogx(db2pow(probePower(1,1:probeDetuningI)-101.3-30)*Energy./(h_p.*PumpGenFreqCent/2), [HPProbeOffSTD(1,1:probeDetuningI).' LPProbeOff(1,1:probeDetuningI).']);
% % semilogx(db2pow(probePower(1,1:probeDetuningI)-116-10)*2e-5./(h_p.*PumpGenFreqCent/2), [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
%  title ('OFF 2*\sigma of MEANS AND VARS');
% % title ('ON-OFF ratios')
% xlabel('mean N number');
% grid on;
% probabfig=figure(1146);
% clf;
% subplot(2,2,3)
% plot(probeDetuning(1,1:pumpI), [ ClickProbVarTrue(1,1:pumpI).'  ClickProbVarFalse(1,1:pumpI).' ClickProbMeanTrue(1,1:pumpI).'  ClickProbMeanFalse(1,1:pumpI).'],'Linewidth',3);
% title ('True and False Positive Rates');
% legend('Variance TPR','Variance FPR', 'Mean TPR', 'Mean FPR');
% ylabel('probability');
% xlabel('detuning');
% set(gca,'FontSize',14);
% grid on

    tim=toc;
    t_ar=[t_ar tim];
    t_av=mean(t_ar);
    tic
    tot_meas_done=tot_meas_done+1;
    tot_meas_els=tot_meas_all-tot_meas_done;
    fprintf(['One cycle of measurement took ' num2str(tim) ' s,\n Time remaining: ' num2str(t_av.*tot_meas_els/60) ' min \n ----------------------------------------- \n'])
    
end
% savefig(mainfig,[TempDirPath 'Figures\' MeasName,'_VarsMeans_' num2str(probeLength*1e6) 'us.fig']);
% savefig(probabfig,[TempDirPath 'Figures\' MeasName,'_Probabilities_' num2str(VoltageGeneratorV) 'V_' num2str(probeLength*1e6) 'us.fig']);
%   save([TempDirPath 'Data\' MeasName,'_' num2str(probeState) '_' num2str(FIXpump(pumpI))  '_' num2str(VoltageGeneratorV) '.mat'])
%% Close devices
fclose(TriggerGenerator);
fclose(SignalGenerator);

% fprintf(VoltageGenerator, 'EXON 0');
fclose(VoltageGenerator);
% instrreset;



    
    
    