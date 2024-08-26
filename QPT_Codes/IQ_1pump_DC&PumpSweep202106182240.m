
clearvars
close all
    
%%
DeviceName = 'QWJPA_v2_1';
DirPath='C:\MatlabPrograms\PetrovninK\';
cd('C:\MatlabPrograms\PetrovninK\');
addpath('C:\Matlabprograms\PetrovninK\functions_for_calc');
addpath('C:\Matlabprograms\PetrovninK\functions');

filename=todayfilename;
MeasName = ['IQ_1pump_DC&PumpSweep' filename];
TempDirPath=['G:\Kirill\' DeviceName '\' datestr(date) '\' MeasName '\'];
mkdir(TempDirPath);
mkdir([TempDirPath 'Data\']);
mkdir([TempDirPath 'Figures\']);
addpath(TempDirPath);
CurrentFileName=mfilename('fullpath');



copyfile([CurrentFileName,'.m'],[TempDirPath,MeasName,'.m']);

%  pumpState=1;
 
%% Set parameters and initalize device
atten = 0;%Gain of pump amplifier
pumpRaw = [linspace(0.20,0.7,81)]-atten;%linspace(-20,-2,10) [linspace(-41,-6,18) linspace(-5,7,21) linspace(7.2,16,25) linspace(16.5,20.5,9)]-atten;%[linspace(-20,-2,10) linspace(0,15,16)]-atten;%%linspace(-41,-6,18)  linspace(-2,18,21)-atten;[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
pump=pumpCalibration(pumpRaw);%linspace(-20,-2,10) [linspace(-41,-6,18) linspace(-5,7,21) linspace(7.2,16,25) linspace(16.5,20.5,9)]-atten;%[linspace(-20,-2,10) linspace(0,15,16)]-atten;%%linspace(-41,-6,18)  linspace(-2,18,21)-atten;[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
% phase_offs=0:10:360;
PumpGenFreqCent = 12.084e9;%12.0778e9;%11.995e9;
MeasFreq = 14e6;% 6052499930.9;%PumpGenFreqCent/2-6.052533066e9+PumpGenFreqCent/2;%PumpGenFreqCent/2   
RBW = 2e6; 
SGCh2 = 3;
SGCh3 = 1; %LO oscillator channel
n_ave = 1;
GAIN = db2pow(86.6); %Quadrature gain 91.6139dB +- 0.15116dB%Quadrature noise temperature 5.1895K coth(hw/2kT)=35.9
h_p = 6.62607004e-34;
VoltageGeneratorV = 5.8+linspace(-0.35,0.35,401);%+linspace(-0.02,0.02,51);
%(2^7);

pump_f_shift1=0;

tot_phase=0/pi*180;

rotationUseFlag=1;

% CONFIGURE ADQ14 digitizer
boardid = 1;
samplerate = 50e6;
nofsamples = 3e6;
nof_records = 1;
trig_delay=0e-3;
channels_mask = 3; %Channels 1&2
trigger = 2;%0 - SW trigger, 1 - Internal, 2 - External and 3 - Lvl trigger
MeasLength=nofsamples/samplerate+1;

[samplerate]=ADQ14_init(boardid,nofsamples,nof_records,channels_mask,1e9/samplerate,trigger,trig_delay);%Init and Read precise value of samplerate


TriggerGenerator = visa('agilent','TCPIP0::169.254.58.12::inst0::INSTR');
old_obj = instrfind('type','visa-tcpip','RemoteHost','169.254.58.12');
if(~isempty(old_obj))
    fclose(old_obj);
end
fopen(TriggerGenerator);

%CONFIGURE TRIGGER GENERATOR
fprintf(TriggerGenerator,'ROSC:SOUR:EXT');
fprintf(TriggerGenerator,'SOUR1:VOLT:UNIT VPP');
fprintf(TriggerGenerator,'SOUR2:VOLT:UNIT VPP');
fprintf(TriggerGenerator,['SOUR1:APPL:SIN ' num2str(MeasFreq*2) ',0.001,0']);
fprintf(TriggerGenerator,['SOUR2:APPL:SIN ' num2str(MeasFreq*2+1e6) ',0.001,0']);
fprintf(TriggerGenerator,'OUTP1 ON');
fprintf(TriggerGenerator,'OUTP2 OFF');
fprintf(TriggerGenerator,['OUTP1:LOAD 50']);
fprintf(TriggerGenerator,['OUTP2:LOAD 50']);
fprintf(TriggerGenerator,'SOUR1:BURS:MOD TRIG');
fprintf(TriggerGenerator,['SOUR1:BURS:NCYC ' num2str((MeasLength*MeasFreq*2)) ]);
fprintf(TriggerGenerator,'SOUR2:BURS:MOD TRIG');
fprintf(TriggerGenerator,['SOUR2:BURS:NCYC ' num2str((MeasLength*MeasFreq*2+0)) ]);
fprintf(TriggerGenerator,'UNIT:ANGLE DEG');
fprintf(TriggerGenerator,['SOUR1:BURS:PHAS ' num2str(pump_f_shift1) ]);
fprintf(TriggerGenerator,['SOUR2:BURS:PHAS ' num2str(0) ]);
fprintf(TriggerGenerator,'OUTP:TRIG 1')
fprintf(TriggerGenerator,'TRIG1:SOUR BUS');
fprintf(TriggerGenerator,'TRIG2:SOUR BUS');
fprintf(TriggerGenerator,'SOUR1:BURS:STAT ON');
fprintf(TriggerGenerator,'SOUR2:BURS:STAT ON');
fprintf(TriggerGenerator,'OUTP1 ON');
fprintf(TriggerGenerator,'OUTP2 ON');


VoltageGenerator = visa('agilent','TCPIP0::169.254.58.15::inst0::INSTR');
old_obj = instrfind('type','visa-tcpip','RemoteHost','169.254.58.15');
        if(~isempty(old_obj))
        fclose(old_obj);
        end
fopen(VoltageGenerator);

% CONFIGURE DC generator
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
createListAnapico(SignalGenerator,fileNameCh3,SGCh3,PumpGenFreqCent/2-MeasFreq,18,MeasLength+0.2,0)
createListAnapico(SignalGenerator,fileNameCh2,SGCh2,PumpGenFreqCent-2*MeasFreq,15,MeasLength+0.2,0);
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
% SignalGeneratorPumpP = pump; 
tot_meas_all=length(VoltageGeneratorV).*length(pump).*(n_ave);
tot_meas_done=0;
t_ar=[];
tic

for Voltage_i=1:length(VoltageGeneratorV)
for ii=1:length(pump)
    
    for pumpState=0:1

        fprintf(VoltageGenerator,['SOUR' num2str(DCCh) ':VOLT:OFFS ' num2str(VoltageGeneratorV(Voltage_i))]);
        if pumpState==1
            fprintf(TriggerGenerator,['SOUR1:VOLT ' num2str(pumpRaw(ii))]);
            fprintf(TriggerGenerator,'OUTP1 ON');
        else
            fprintf(TriggerGenerator,'SOUR1:VOLT 0.001');

        end

        for n_a=1:n_ave
            
%            load([TempDirPath 'Data\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '.mat'],'B1')
            fprintf(SignalGenerator,'INIT');
            pause(0.2)
            IQ_raw = (ADQ14_acq(boardid,nofsamples,nof_records,channels_mask,trigger,2,TriggerGenerator));


            %%

            t_s=(0:nofsamples-1)/samplerate;
            IQ_cap=exp(1i.*2*pi.*t_s.*(MeasFreq)+1i.*pi*tot_phase/180).*((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))-mean((IQ_raw(1,:)+1i.*exp(1i.*8/180*pi).*IQ_raw(2,:))))...
                ./sqrt(h_p.*(PumpGenFreqCent)/2.*RBW.*50.*(GAIN));            


%             N=length(IQ_cap);
%             FFT = fftshift(fft(complex(IQ_raw(1,:),IQ_raw(2,:))))./(N-1);
% 
%             f_s=linspace(-samplerate/2,samplerate/2-samplerate/N,N);
%             figure(1);
%             clf;
%             plot(f_s,abs(FFT));
%             clearvars IQ_raw

            fcut=0e1;
            n_modes=1;

            if pumpState == 0   %Calibration Covariance Matrix Calculation
            [B1]=IQfiltering(IQ_cap,samplerate,0,n_modes, [-RBW/2 RBW/2]);%Mode filtering
            B1=B1.';
%           
                           
                             L = [real((B1(:,1))), imag((B1(:,1)))];

                             CovMat_init(:,:,n_a,ii,Voltage_i)=cov(L);
                             CovMat_init_ave(:,:,ii,Voltage_i)=mean(CovMat_init(:,:,1:n_a,ii,Voltage_i),3);
                             varI_init(ii,n_a,Voltage_i)=var(real(B1));%Variance calculations  [-RBW/2 RBW/2]-MeasFreq ,[-RBW/2 RBW/2]+MeasFreq 
                             varQ_init(ii,n_a,Voltage_i)=var(imag(B1));
                             varI_ave_init(ii,Voltage_i)=mean(varI_init(ii,1:n_a,Voltage_i),2);
                             varQ_ave_init(ii,Voltage_i)=mean(varQ_init(ii,1:n_a,Voltage_i),2);
                             meanI_init(ii,n_a,Voltage_i)=abs(mean(real(B1)));
                             meanQ_init(ii,n_a,Voltage_i)=abs(mean(imag(B1)));
                             meanI_ave_init(ii,Voltage_i)=mean(meanI_init(ii,1:n_a,Voltage_i),2);
                             meanQ_ave_init(ii,Voltage_i)=mean(meanQ_init(ii,1:n_a,Voltage_i),2);
                             P_init(ii,Voltage_i)=mean(abs(B1).^2);
            %                  A1=sqrt(samplerate/(RBW)).*IQfiltering(IQ_cap,samplerate,0,1,[-RBW/2 RBW/2]+MeasFreq).';%.*exp(sum(IQphi(:))/2);[-RBW/2 -fcut]+MeasFreq ,[fcut RBW/2]+MeasFreq  
%                              A1=(B1(:,1)+B1(:,2))./sqrt(2);
%                              varI_init(ii,n_a)=var(real(A1));%Variance calculations  [-RBW/2 RBW/2]-MeasFreq ,[-RBW/2 RBW/2]+MeasFreq 
%                              varQ_init(ii,n_a)=var(imag(A1));
%                              varI_ave_init(ii)=mean(varI_init(ii,1:n_a),2);
%                              varQ_ave_init(ii)=mean(varQ_init(ii,1:n_a),2);
            %                  A0=IQ_cap;


            elseif pumpState==1
           % Define angle here..     
           [B1]=IQfiltering(IQ_cap,samplerate,0,2,[-RBW/2 0] ,[0 RBW/2]);%Mode filtering
              B1=B1.';
            for dd=1:1
              L11 = [real((B1(:,1))), imag((B1(:,1))), real(B1(:,2)), imag(B1(:,2))];
                CovMat1(:,:,n_a,ii)=cov(L11);
                C21=CovMat1(1:2,3:4,n_a,ii);
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
            clearvars B1 B11 CovMat1 C21 L11
            %   
                
                
            [B1]=IQfiltering(IQ_cap,samplerate,0,n_modes, [-RBW/2 RBW/2]);%Mode filtering
            B1=B1.'.*exp(sum(IQphi(:))/2);
           

             

            
                L = [real((B1(:,1))), imag((B1(:,1)))];

                CovMat(:,:,n_a,ii,Voltage_i)=cov(L);
                CovMat_ave(:,:,ii,Voltage_i)=mean(CovMat(:,:,1:n_a,ii,Voltage_i),3);
                CovMat_calc_ave(:,:,ii,Voltage_i)=4*(CovMat_ave(:,:,ii,Voltage_i)-CovMat_init_ave(:,:,ii,Voltage_i))+eye(2);%Normalization and subtraction of Preamp noise (Calibration matrix)
%             B1=B1;   
            if mean(abs(B1))>15
               B1=B1.*exp(-1i*atan(mean(imag(B1))./mean(real(B1))));
            end 
            AAA1 = [-1:0.01:1]*max(abs(B1));
            [MM1 CC1] = hist3([real(B1) imag(B1)], {AAA1 AAA1});
            varI(ii,n_a,Voltage_i)=var(real(B1));%max(eig(CovMat_calc_ave(:,:,ii,Voltage_i))/4);
            varQ(ii,n_a,Voltage_i)=var(imag(B1));%min(eig(CovMat_calc_ave(:,:,ii,Voltage_i))/4);
            varI_ave(ii,Voltage_i)=mean(varI(ii,1:n_a,Voltage_i),2);
            varQ_ave(ii,Voltage_i)=mean(varQ(ii,1:n_a,Voltage_i),2);
            meanI(ii,n_a,Voltage_i)=abs(mean(real(B1)));
            meanQ(ii,n_a,Voltage_i)=abs(mean(imag(B1)));
            meanI_ave(ii,Voltage_i)=mean(meanI(ii,1:n_a,Voltage_i),2);
            meanQ_ave(ii,Voltage_i)=mean(meanQ(ii,1:n_a,Voltage_i),2);
            P(ii,Voltage_i)=mean(abs(B1).^2);
            hhh=figure(2);
            clf;
            surface(AAA1, AAA1, MM1','edgecolor', 'none');

            axis square
            axis tight
            end

               hhh2=figure (51);%Showing covariance matrix
                    clf;
                    if pumpState==0
                        bar3(4*CovMat_init_ave(:,:,ii,Voltage_i));
                        view(-30, 15)    ;
                        title('CALIBRATION');
                        set(gca,'XTickLabel',{'I_1','Q_1','I_2','Q_2','I3','Q3'});
                        set(gca,'YTickLabel',{'I_1','Q_1','I_2','Q_2','I3','Q3'});
                    else
                        bar3(CovMat_calc_ave(:,:,ii,Voltage_i));
                        view(-30, 15);
                        set(gca,'XTickLabel',{'I_1','Q_1','I_2','Q_2','I3','Q3'});
                        set(gca,'YTickLabel',{'I_1','Q_1','I_2','Q_2','I3','Q3'});
                        title(['Pump: ' num2str(pump(ii)+atten) ', Voltage: ' num2str(VoltageGeneratorV(Voltage_i))]);
                    end



             end
            clearvars IQ A1 IQ_init IQ_cap AAA  L1 B11 A0 L t_s dem IQ_raw

%          save([TempDirPath 'Data\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '.mat'])
%         if pumpState==1
%         end
%          savefig(hhh2,[TempDirPath 'Figures\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '_COVMAT.fig'])
        %%
        if pumpState==1
                     savefig(hhh,[TempDirPath 'Figures\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '.fig']);

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
    %%
        if (ii>1)&&(Voltage_i>1)
            hhh6=figure (118);
             clf;

            % hold on
            s=surf((VoltageGeneratorV(1:Voltage_i).'+4+2.35)/33.5,pump(1:ii),0.5+(P(1:ii,1:Voltage_i)-P_init(1:ii,1:Voltage_i)));
            set(s,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
%             shading interp
            legend('N-number')
            colorbar
            % title('S-criteria');\
            hhh7=figure (121);
             clf;
            s1=subplot(2,2,1)
            ss1=surf((VoltageGeneratorV(1:Voltage_i).'+4+2.35)/33.5,pump(1:ii),0.25+varI_ave(1:ii,1:Voltage_i)-varI_ave_init(1:ii,1:Voltage_i));
            set(ss1,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
%             shading interp
            legend('Variance I')
            colorbar
            s2=subplot(2,2,2)
            ss2=surf((VoltageGeneratorV(1:Voltage_i).'+4+2.35)/33.5,pump(1:ii),0.25+varQ_ave(1:ii,1:Voltage_i)-varQ_ave_init(1:ii,1:Voltage_i));
            set(ss2,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
%             shading interp
            legend('Variance Q')
            colorbar
            
            s3=subplot(2,2,3)
            ss3=surf((VoltageGeneratorV(1:Voltage_i).'+4+2.35)/33.5,pump(1:ii),meanI_ave(1:ii,1:Voltage_i)-meanI_ave_init(1:ii,1:Voltage_i));
            set(ss3,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
%             shading interp
            legend('Mean I')
            colorbar
            
            s4=subplot(2,2,4)
            ss4=surf((VoltageGeneratorV(1:Voltage_i).'+4+2.35)/33.5,pump(1:ii),meanQ_ave(1:ii,1:Voltage_i)-meanQ_ave_init(1:ii,1:Voltage_i));
            set(ss4,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
%             shading interp
            legend('Mean Q')
            colorbar
            
     
            
        end

        
end
    save([TempDirPath 'Data\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '.mat']);
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



    
    
    