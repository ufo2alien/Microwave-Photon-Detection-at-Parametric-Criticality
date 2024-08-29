
clear all
close all
    
%%
DeviceName = 'QWJPA_v2_2';
DirPath='C:\MatlabPrograms\PetrovninK\';
cd('C:\MatlabPrograms\PetrovninK\');
addpath('C:\Matlabprograms\PetrovninK\functions');

filename=todayfilename;
MeasName = ['Probe_Detection',filename];
TempDirPath=['D:\Kirill\' DeviceName '\' datestr(date) '\' MeasName '\'];
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
NAdet=0.0e6;
FIXpump=6.90;%[linspace(-41,-6,36) linspace(-5,7,21) linspace(7.2,13,25) linspace(13,20.5,11)]-atten;%[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
probePower=linspace(-40,-10,31);

% pump = [linspace(0.5,0.8,11)]-atten;%linspace(-20,-2,10) [linspace(-41,-6,18) linspace(-5,7,21) linspace(7.2,16,25) linspace(16.5,20.5,9)]-atten;%[linspace(-20,-2,10) linspace(0,15,16)]-atten;%%linspace(-41,-6,18)  linspace(-2,18,21)-atten;[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
PumpGenFreqCent = 12.084e9;%12.0778e9;%11.995e9;
MeasFreq = 14e6;% 6052499930.9;%PumpGenFreqCent/2-6.052533066e9+PumpGenFreqCent/2;%PumpGenFreqCent/2   
RBW = 1e6; 
SGCh2 = 2;
SGCh3 = 1; %READOUT LO oscillator channel
n_ave = 1;
GAIN = db2pow(105-18.2); %-16.42Quadrature gain 90.7 +- 0.2dB%Quadrature noise temperature 5.11K coth(hw/2kT)=35.9
h_p = 6.62607004e-34;
VoltageGeneratorV = 0.739;
samplerate = 50e6;%(2^7);
kb = 1.380649e-23;
tot_phase=-30.5-56.5/2-43;

% load('D:\Kirill\QWJPA_v2_1\20-Apr-2021\IQ_1pump_Digitizer_w_Triggering202104200138\Data\IQ_1pump_Digitizer_w_Triggering202104200138pumpState_1_1_0.7565.mat','Noise_Gain');
%% CONFIGURE ADQ14 digitizer
boardid = 1;
nofsamples = 5e6;
nof_records = 1;
trig_delay=0e-3;
channels_mask = 3; %Channels 1&2
trigger = 2;%0 - SW trigger, 1 - Internal, 2 - External and 3 - Lvl trigger
MeasLength=nofsamples/samplerate*2;

[samplerate]=ADQ14_init(boardid,nofsamples,nof_records,channels_mask,1e9/samplerate,trigger,trig_delay);%Init and Read precise value of samplerate

%% CONFIGURE TRIGGER GENERATOR
TriggerGenerator = visa('agilent','TCPIP0::169.254.58.12::inst0::INSTR');
old_obj = instrfind('type','visa-tcpip','RemoteHost','169.254.58.12');
if(~isempty(old_obj))
    fclose(old_obj);
end
fopen(TriggerGenerator);


fprintf(TriggerGenerator,'ROSC:SOUR EXT');
fprintf(TriggerGenerator,'SOUR1:VOLT:UNIT VPP');
fprintf(TriggerGenerator,'SOUR2:VOLT:UNIT VPP');
fprintf(TriggerGenerator,['SOUR1:APPL:SIN ' num2str(MeasFreq*2) ',0.001,0']);
fprintf(TriggerGenerator,['SOUR2:APPL:SIN ' num2str(MeasFreq*2+1e6) ',0.001,0']);
fprintf(TriggerGenerator,'OUTP1 OFF');
fprintf(TriggerGenerator,'OUTP2 OFF');
fprintf(TriggerGenerator,['OUTP1:LOAD 50']);
fprintf(TriggerGenerator,['OUTP2:LOAD 50']);
fprintf(TriggerGenerator,'SOUR1:BURS:MOD TRIG');
fprintf(TriggerGenerator,['SOUR1:BURS:NCYC ' num2str(MeasLength*MeasFreq*2) ]);
fprintf(TriggerGenerator,'SOUR2:BURS:MOD TRIG');
fprintf(TriggerGenerator,['SOUR2:BURS:NCYC ' num2str(MeasLength*MeasFreq*2) ]);
fprintf(TriggerGenerator,'UNIT:ANGLE DEG');
fprintf(TriggerGenerator,'SOUR2:BURS:PHAS 0');
fprintf(TriggerGenerator,'SOUR2:BURS:PHAS 0');
fprintf(TriggerGenerator,'OUTP:TRIG 1')
fprintf(TriggerGenerator,'TRIG1:SOUR BUS');
fprintf(TriggerGenerator,'TRIG2:SOUR BUS');
fprintf(TriggerGenerator,'SOUR1:BURS:STAT ON');
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
fprintf(VoltageGenerator,['OUTP' num2str(DCCh) ':LOAD INF']);
fprintf(VoltageGenerator,['OUTP' num2str(DCCh2) ':LOAD INF']);
fprintf(VoltageGenerator,['SOUR' num2str(DCCh) ':VOLT:OFFS ' num2str(VoltageGeneratorV)]);
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
fprintf(SignalGenerator, ['ROSC:SOUR EXT']);
fprintf(SignalGenerator, ['SOURce:LFO:STAT 0']);
fileNameCh2='Pump';
fileNameCh3='LOFreq';
deleteListAnapico(SignalGenerator,fileNameCh2);
deleteListAnapico(SignalGenerator,fileNameCh3);
createListAnapico(SignalGenerator,fileNameCh3,SGCh3,PumpGenFreqCent/2-MeasFreq,13,MeasLength,0)
createListAnapico(SignalGenerator,fileNameCh2,SGCh2,PumpGenFreqCent,FIXpump,MeasLength,0);
loadListAnapico(SignalGenerator,fileNameCh2,SGCh2);
loadListAnapico(SignalGenerator,fileNameCh3,SGCh3);
disp(query(SignalGenerator,['MEM' num2str(SGCh2) ':FILE:LIST:DATA?']));
disp(query(SignalGenerator,['MEM' num2str(SGCh3) ':FILE:LIST:DATA?']));
fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':FREQ:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':POW:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':FREQ:MODE LIST']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':POW:MODE LIST']);
fprintf(SignalGenerator, ['INIT:CONT 0']);
fprintf(SignalGenerator, ['TRIG:DEL 0e-3']);
fprintf(SignalGenerator, ['TRIG:SOUR EXT']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':LIST:COUNT 1']);
fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':LIST:COUNT 1']);
fprintf(SignalGenerator, ['OUTP' num2str(SGCh2) ':STAT ', num2str(1)]); 
fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(1)]); 
pause(0.1)




NetworkAnalyzer = visa('agilent','TCPIP0::169.254.58.11::inst0::INSTR');
old_obj = instrfind('type','visa-tcpip','RemoteHost','169.254.58.11');
        if(~isempty(old_obj))
        fclose(old_obj);
        end
set(NetworkAnalyzer, 'InputBufferSize', 2^20);
fopen(NetworkAnalyzer);
query(NetworkAnalyzer,'*IDN?')
pause(0.2)
%% CHECK

tot_meas_all=length(probePower).*(n_ave);
tot_meas_done=0;
t_ar=[];
tic

% Creating ZERO arrays for timetraced data A1  AAA B1 L1 B11 A0 L t_s dem IQ_raw g G h H
n_modes=2;
IQ_raw = zeros(2,nofsamples);
IQ_cap = zeros(1,nofsamples);
B1 = zeros(nofsamples,1);
% B11  = zeros(nofsamples,2);
% L1 = zeros(nofsamples,4);


for probePowerI=1:length(probePower)
for ii=1:2%TURN ON-OFF probe signal
for pumpState=1

% fprintf(TriggerGenerator,['SOUR1:VOLT ' num2str(FIXpump)]);
fprintf(NetworkAnalyzer, ['OUTP1 ' num2str((ii-1))]);
fprintf(NetworkAnalyzer, ['SOUR1:POW ' num2str(probePower(probePowerI))]);
fprintf(NetworkAnalyzer, 'SENS1:BAND 1000');
fprintf(NetworkAnalyzer, ['SENS1:SWEEP:TYPE CW']);
fprintf(NetworkAnalyzer, ['SOUR1:FREQ ' num2str(PumpGenFreqCent/2+NAdet)]);
fprintf(NetworkAnalyzer, 'INIT1:CONT ON');
fprintf(NetworkAnalyzer, ['SENS1:SWE:POIN ', num2str(1)]);

for n_a=1:n_ave
    
%% Make measurement

    fprintf(SignalGenerator,'INIT');
    pause(0.2)
    IQ_raw = (ADQ14_acq(boardid,nofsamples,nof_records,channels_mask,trigger,2,TriggerGenerator));
    t_s=(0:nofsamples-1)/samplerate;
    IQ_cap=(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*pi*tot_phase/180).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
        ./sqrt(h_p.*(PumpGenFreqCent)/2.*RBW.*50.*(GAIN)));
  B1=(IQfiltering(IQ_cap,samplerate,0,1,[-RBW/2 RBW/2])).';%Mode filtering

if mean(abs(B1))>5
    if mean(real(B1))>5
        B1=B1.*exp(-1i*atan(mean(imag(B1))./mean(real(B1))));
    else
        B1=-B1.*exp(-1i*atan(mean(imag(B1))./mean(real(B1))));
    end
end
% B1=B1'.*2./sqrt(h_p.*MeasFreq.*RBW.*50.*GAIN);
AAA1 = [-1:0.01:1]*max(abs(B1));
[MM1 CC1] = hist3([real(B1) imag(B1)], {AAA1 AAA1});
[pks,locs]=findpeaks(max(MM1.'),AAA1,'MinPeakDistance',0.6*max(AAA1),'MinPeakHeight',20,'SortStr','descend','Npeaks',2);
hhh=figure(37);


if ii==2
    subplot(1,2,2)
    surface(AAA1, AAA1, MM1','edgecolor', 'none');
    axis square
    axis tight
    view(0,0);
    title('PROBE ON')
HPProbeOn(probePowerI)=pks(1);
LPProbeOn(probePowerI)=pks(2);
elseif ii==1
    clf;
    subplot(1,2,1)
    surface(AAA1, AAA1, MM1','edgecolor', 'none');
    axis square
    axis tight
    view(0,0);
    title('PROBE OFF')
HPProbeOff(probePowerI)=pks(1);
LPProbeOff(probePowerI)=pks(2);
end
% varI(ii,n_a)=var(real(A1));
% varQ(ii,n_a)=var(imag(A1));
% varI_ave(ii)=mean(varI(ii,1:n_a),2);
% varQ_ave(ii)=mean(varQ(ii,1:n_a),2);
% pur(ii,1)=1/det(CovMat_calc_ave(:,:,ii));


end


end
 save([TempDirPath 'Data\' MeasName, '_' num2str(probePower(probePowerI))  '_' num2str(VoltageGeneratorV) '.mat'], '-regexp', '^(?!(IQ|A1|IQ_init|IQ_cap|AAA|B1|L1|B11|A0|L|t_s|dem|IQ_raw|g|G|h|H)$).')
       
if pumpState==1
 savefig(hhh,[TempDirPath,MeasName,'_' num2str(probePower(probePowerI)) '_' num2str(VoltageGeneratorV) '.fig'])
end

end
figure(1144)
clf;
subplot(2,2,1)
plot(probePower(1,1:probePowerI)-116, [HPProbeOn(1,1:probePowerI).' LPProbeOn(1,1:probePowerI).' HPProbeOff(1,1:probePowerI).' LPProbeOff(1,1:probePowerI).']);
title ('ON-OFF distibutions maxima')
% figure(1145)
% clf;
subplot(2,2,2)
plot(probePower(1,1:probePowerI)-116, [HPProbeOn(1:probePowerI).'./HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'./LPProbeOff(1:probePowerI).']);
title ('ON-OFF ratios')
% figure(1146)
% clf;
subplot(2,2,3)
plot(probePower(1,1:probePowerI)-116, [HPProbeOn(1:probePowerI).'-HPProbeOff(1:probePowerI).' LPProbeOn(1:probePowerI).'-LPProbeOff(1:probePowerI).']);
title ('ON-OFF diffs')

    tim=toc;
    t_ar=[t_ar tim];
    t_av=mean(t_ar);
    tic
    tot_meas_done=tot_meas_done+1;
    tot_meas_els=tot_meas_all-tot_meas_done;
    fprintf(['One cycle of measurement took ' num2str(tim) ' s,\n Time remaining: ' num2str(t_av.*tot_meas_els/60) ' min \n ----------------------------------------- \n'])
    
end
%% Close devices
fprintf(SignalGenerator,['OUTP' num2str(SGCh2) ':STAT 0']);
fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(0)]);  
fclose(TriggerGenerator);
fclose(SignalGenerator);

% fprintf(VoltageGenerator, 'EXON 0');
fclose(VoltageGenerator);
instrreset;



    
    
    