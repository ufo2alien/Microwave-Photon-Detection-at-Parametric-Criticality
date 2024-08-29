
clearvars
% close all
    
%%
DeviceName = 'QWJPA_v2_1';
DirPath='C:\MatlabPrograms\PetrovninK\';
cd('C:\MatlabPrograms\PetrovninK\');
addpath('C:\Matlabprograms\PetrovninK\functions_for_calc');
addpath('C:\Matlabprograms\PetrovninK\functions');

filename=todayfilename;
MeasName = ['IQ_1pump_Digitizer_w_Triggering',filename];
TempDirPath=['I:\Kirill\' DeviceName '\' datestr(date) '\' MeasName '\'];
mkdir(TempDirPath);
mkdir([TempDirPath 'Data\']);
mkdir([TempDirPath 'Figures\']);
addpath(TempDirPath);
CurrentFileName=mfilename('fullpath');



copyfile([CurrentFileName,'.m'],[TempDirPath,MeasName,'.m']);

%  pumpState=1;
 
%% Set parameters and initalize device
atten = 0;%Gain of pump amplifier
pumpRaw = 0.3;%[linspace(0.001,0.6,21)]-atten;%linspace(-20,-2,10) [linspace(-41,-6,18) linspace(-5,7,21) linspace(7.2,16,25) linspace(16.5,20.5,9)]-atten;%[linspace(-20,-2,10) linspace(0,15,16)]-atten;%%linspace(-41,-6,18)  linspace(-2,18,21)-atten;[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
pump=pumpCalibration(pumpRaw);
PumpGenFreqCent = 12.084e9;%12.0778e9;%11.995e9;
MeasFreq = 14e6;% 6052499930.9;%PumpGenFreqCent/2-6.052533066e9+PumpGenFreqCent/2;%PumpGenFreqCent/2   
RBW = 2e6; 
SGCh2 = 3;
SGCh3 = 1; %LO oscillator channel
n_ave = 1;
GAIN = db2pow(86.85); %-16.42Quadrature gain 90.7 +- 0.2dB%Quadrature noise temperature 5.11K coth(hw/2kT)=35.9
h_p = 6.62607004e-34;
VoltageGeneratorV = 5.78;0.7405;0.755;
samplerate = 50e6;%(2^7);
kb = 1.380649e-23;
tot_phase=-30.5-56.5/2-43;

% load('D:\Kirill\QWJPA_v2_1\20-Apr-2021\IQ_1pump_Digitizer_w_Triggering202104200138\Data\IQ_1pump_Digitizer_w_Triggering202104200138pumpState_1_1_0.7565.mat','Noise_Gain');
%% CONFIGURE ADQ14 digitizer
boardid = 1;
nofsamples = 1e7;
nof_records = 1;
trig_delay=0e-3;
channels_mask = 3; %Channels 1&2
trigger = 2;%0 - SW trigger, 1 - Internal, 2 - External and 3 - Lvl trigger
MeasLength=nofsamples/samplerate+1;

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
fprintf(TriggerGenerator,'OUTP1 ON');
fprintf(TriggerGenerator,'OUTP2 ON');
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
fprintf(VoltageGenerator,['SOUR' num2str(DCCh) ':VOLT:OFFS ' num2str(VoltageGeneratorV(1))]);
fprintf(VoltageGenerator,['SOUR' num2str(DCCh2) ':VOLT:OFFS ' num2str((-4))]);
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
createListAnapico(SignalGenerator,fileNameCh3,SGCh3,PumpGenFreqCent/2-MeasFreq,18,MeasLength,0)
createListAnapico(SignalGenerator,fileNameCh2,SGCh2,PumpGenFreqCent-2*MeasFreq,15,MeasLength,0);
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
        
%% CHECK
tot_meas_all=length(pumpRaw).*(n_ave);
tot_meas_done=0;
t_ar=[];
tic

% Creating ZERO arrays for timetraced data A1  AAA B1 L1 B11 A0 L t_s dem IQ_raw g G h H
n_modes=2;
IQ_raw = zeros(2,nofsamples);
IQ_cap = zeros(1,nofsamples);
B1 = zeros(nofsamples,2);
B11  = zeros(nofsamples,2);
L1 = zeros(nofsamples,4);

CovMat_init = zeros(4,4,n_ave,length(pump));
CovMat_init_ave = zeros(4,4,length(pump));
CovMat = zeros(4,4,n_ave,length(pump));
CovMat_ave = zeros(4,4,length(pump));
CovMat_calc_ave = zeros(4,4,length(pump));

EIG = zeros(length(pump),1);
for ii=1:length(pump)
    for pumpState=0:1

            
%         fprintf(SignalGenerator, ['OUTP' num2str(SGCh2) ':STAT ', num2str(0)]); 
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':FREQ:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':POW:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':PHASE:MODE CW']);
% 
%         fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(0)]); 
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':FREQ:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':POW:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':PHASE:MODE CW']);


        % pause(0.1)
        if pumpState==1
            fprintf(TriggerGenerator,['SOUR1:VOLT ' num2str(pumpRaw(ii))]);
%             fprintf(TriggerGenerator,'SOUR2:VOLT 0.001');
            fprintf(TriggerGenerator,'OUTP1 ON');
            fprintf(TriggerGenerator,'OUTP2 OFF');
        
        else
%             fprintf(TriggerGenerator,'SOUR1:VOLT 0.001');
%             fprintf(TriggerGenerator,'SOUR2:VOLT 0.001');
            fprintf(TriggerGenerator,'OUTP1 OFF');
            fprintf(TriggerGenerator,'OUTP2 OFF');
        end


        

%         fprintf(SignalGenerator, ['OUTP' num2str(SGCh2) ':STAT ', num2str(0)]); 
%         fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(0)]); 
%             fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':PHASE:REF']);
%             fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':PHASE:REF']);
        for n_a=1:n_ave
            
           
            fprintf(SignalGenerator,'INIT');
            pause(0.2)
            IQ_raw = (ADQ14_acq(boardid,nofsamples,nof_records,channels_mask,trigger,2,TriggerGenerator));


            %%

            t_s=(0:nofsamples-1)/samplerate;
            IQ_cap=exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*pi*tot_phase/180).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
                ./sqrt(h_p.*(PumpGenFreqCent)/2.*RBW.*50.*(GAIN));%(exp(1i.*2*pi.*t_s.*(MeasFreq)).*IQ_raw(1,:)+ exp(1i.*2*pi.*t_s.*(MeasFreq)).*1i.*(IQ_raw(2,:)))./sqrt(h_p.*(PumpGenFreqCent)/2.*1e9.*50.*GAIN);%
%            if pumpState==0
% %             A1=IQfiltering(IQ_cap,samplerate,0,1, [-RBW/2 RBW/2] );
% %             var(A1)
% %             A1=IQfiltering(exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*pi*tot_phase/180).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:)))),samplerate,0,1, [-RBW/2 RBW/2] );
% %             var(A1)/50
%            end
            %             figure(3)
%             clf;
%             plot(t_s,[IQ_raw(1,:).' IQ_raw(2,:).']);

%             N=length(IQ_cap);
%             FFT = fftshift(fft(complex(IQ_raw(1,:),IQ_raw(2,:))))./(N-1);
% 
%             f_s=linspace(-samplerate/2,samplerate/2-samplerate/N,N);
%             figure(1);
%             clf;
%             plot(f_s,abs(FFT));
%             clearvars IQ_raw

            fcut=0e1;
           
            if pumpState == 0   %Calibration Covariance Matrix Calculation
                             [B1]=sqrt((RBW)/(RBW/2-fcut)).*IQfiltering(IQ_cap,samplerate,0,n_modes, [-RBW/2 -fcut] ,[fcut RBW/2]);%Mode filtering
                             B1=B1.';
                             L1 = [real((B1(:,1))), imag((B1(:,1))), real(B1(:,2)), imag(B1(:,2))]; % 4 column time-vector
                             CovMat_init(:,:,n_a,ii)=cov(L1);
                             CovMat_init_ave(:,:,ii)=mean(CovMat_init(:,:,1:n_a,ii),3);
            %                  A1=sqrt(samplerate/(RBW)).*IQfiltering(IQ_cap,samplerate,0,1,[-RBW/2 RBW/2]+MeasFreq).';%.*exp(sum(IQphi(:))/2);[-RBW/2 -fcut]+MeasFreq ,[fcut RBW/2]+MeasFreq  
                             A1=(B1(:,1)+B1(:,2))./sqrt(2);
                             varI_init(ii,n_a)=var(real(A1));%Variance calculations  [-RBW/2 RBW/2]-MeasFreq ,[-RBW/2 RBW/2]+MeasFreq 
                             varQ_init(ii,n_a)=var(imag(A1));
                             varI_ave_init(ii)=mean(varI_init(ii,1:n_a),2);
                             varQ_ave_init(ii)=mean(varQ_init(ii,1:n_a),2);
            %                  A0=IQ_cap;


            elseif pumpState==1
               [B1]=sqrt((RBW)/(RBW/2-fcut)).*IQfiltering(IQ_cap,samplerate,0,n_modes,[-RBW/2 -fcut] ,[fcut RBW/2]);%Mode filtering
              B1=B1.';
            for dd=1:1
              L = [real((B1(:,1))), imag((B1(:,1))), real(B1(:,2)), imag(B1(:,2))];
                CovMat(:,:,n_a,ii)=cov(L);
                C21=CovMat(1:2,3:4,n_a,ii);
                 angle1(ii)=180*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2)))/pi;
                disp(['Angle is ' num2str(angle1(ii)) ]);
                if (C21(1,1)-C21(2,2))<0
                                 B11(:,1)=-B1(:,1).*exp(-1i*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2))));
                                 IQphi(dd)=1i*(pi-atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2))));
                             else
                                 B11(:,1)=B1(:,1).*exp(-1i*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2))));
                                 IQphi(dd)=-1i.*atan((C21(1,2)+C21(2,1))./(C21(1,1)-C21(2,2)));
                            end
                            B11(:,2)=B1(:,2);
                            B1=B11;
            end
%             B11=B1;
                 L1 = [real((B11(:,1))), imag((B11(:,1))), real(B11(:,2)), imag(B11(:,2))];
                 CovMat1(:,:,n_a,ii)=cov(L1);
                 CovMat_ave(:,:,ii)=mean(CovMat1(:,:,1:n_a,ii),3);
                  CovMat_calc_ave(:,:,ii)=4*(CovMat_ave(:,:,ii)-CovMat_init_ave(:,:,ii))+eye(4);%Normalization and subtraction of Preamp noise (Calibration matrix)
%                   Noise_Gain(:,ii)=diag(CovMat_calc_ave(:,:,ii));
                  omeg=[0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0];
                 Lambda1=diag([1, 1, 1, -1]);
                 test1 = Lambda1*CovMat_calc_ave(:,:,ii)*Lambda1;
                 EIG(ii)=min((eig(CovMat_calc_ave(:,:,ii))));

                 %squeezing

            % [B1]=IQfiltering(IQ_init,samplerate,0,1,[-RBW RBW]);
            % A1=sqrt(samplerate/(RBW)).*IQfiltering(IQ_cap,samplerate,0,1,[-RBW/2 RBW/2]+MeasFreq).'.*exp(sum(IQphi(:))/2);
            % B1=B1'.*2./sqrt(h_p.*MeasFreq.*RBW.*50.*GAIN);
            A1=(B11(:,1)+B11(:,2))./sqrt(2);
            AAA1 = [-1:0.01:1]*max(abs(A1));
            [MM1 CC1] = hist3([real(A1) imag(A1)], {AAA1 AAA1});
            varI(ii,n_a)=var(real(A1));
            varQ(ii,n_a)=var(imag(A1));
            varI_ave(ii)=mean(varI(ii,1:n_a),2);
            varQ_ave(ii)=mean(varQ(ii,1:n_a),2);
            pur(ii,1)=1/det(CovMat_calc_ave(:,:,ii));

            hhh=figure(2);
            clf;
            surface(AAA1, AAA1, MM1','edgecolor', 'none');

            axis square
            axis tight

            end

               hhh2=figure (51);%Showing covariance matrix
                    clf;
                    if pumpState==0
                        bar3(4*CovMat_init_ave(:,:,ii));
                        view(-30, 15)    ;
                        title('CALIBRATION');
                        set(gca,'XTickLabel',{'I_1','Q_1','I_2','Q_2'});
                        set(gca,'YTickLabel',{'I_1','Q_1','I_2','Q_2'});
                    else
                        bar3(CovMat_calc_ave(:,:,ii));
                        view(-30, 15);
                        set(gca,'XTickLabel',{'I_1','Q_1','I_2','Q_2'});
                        set(gca,'YTickLabel',{'I_1','Q_1','I_2','Q_2'});
                        title(['Pump:', num2str(pump(ii)+atten) ,' PPT 2PARTITE CRITERIA VALUE:', num2str(EIG(ii))]);
                    end



             end
            clearvars IQ A1 IQ_init IQ_cap  AAA B1 L1 B11 A0 L t_s IQ_raw t_s

         save([TempDirPath 'Data\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_' num2str(VoltageGeneratorV) '.mat'], '-regexp', '^(?!(IQ|A1|IQ_init|IQ_cap|AAA|B1|L1|B11|A0|L|t_s|dem|IQ_raw|g|G|h|H)$).')
        if pumpState==1
         savefig(hhh,[TempDirPath 'Figures\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_' num2str(VoltageGeneratorV) '.fig'])
        end
         savefig(hhh2,[TempDirPath 'Figures\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_' num2str(VoltageGeneratorV) '_COVMAT.fig'])
        %%
        
        
        if pumpState==1
        figure (114)
        clf;
        % hold on
        plot((pump(1:ii)+atten),EIG(1:ii),'Linewidth',3)

        figure(1144)
        clf;
        plot((pump(1:ii)+atten), 10*log10(4*[(varI_ave(1:ii).'-varI_ave_init(1:ii).')+1/4, (varQ_ave(1:ii).'-varQ_ave_init(1:ii).')+1/4]));
        % figure (1115)
        % plot((pump(1:ii)+atten),pur(1:ii,1),'r')

        end   
    end
    tim=toc;
    t_ar=[t_ar tim];
    t_av=mean(t_ar);
    tic
    tot_meas_done=tot_meas_done+1;
    tot_meas_els=tot_meas_all-tot_meas_done;
    fprintf(['One cycle of measurement took ' num2str(tim) ' s,\n Time remaining: ' num2str(t_av.*tot_meas_els/60) ' min \n ----------------------------------------- \n'])
    
end
%% Close devices
% fprintf(SignalGenerator,['OUTP' num2str(SGCh2) ':STAT 0']);
% fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(0)]);  
% fprintf(SignalGenerator, ['OUTP' num2str(SGCh2) ':STAT ', num2str(0)]); 
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':FREQ:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':POW:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh2) ':PHASE:MODE CW']);
% 
%         fprintf(SignalGenerator, ['OUTP' num2str(SGCh3) ':STAT ', num2str(0)]); 
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':FREQ:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':POW:MODE CW']);
%         fprintf(SignalGenerator, ['SOUR' num2str(SGCh3) ':PHASE:MODE CW']);
fclose(SignalGenerator);
fclose(TriggerGenerator)
% fprintf(VoltageGenerator, 'EXON 0');
fclose(VoltageGenerator);
instrreset;



    
    
    