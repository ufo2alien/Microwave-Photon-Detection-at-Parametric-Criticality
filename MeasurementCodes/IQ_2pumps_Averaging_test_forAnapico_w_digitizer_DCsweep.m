
clearvars
close all
    
%%
DeviceName = 'QWJPA_v2_1';
DirPath='C:\MatlabPrograms\PetrovninK\';
cd('C:\MatlabPrograms\PetrovninK\');
addpath('C:\Matlabprograms\PetrovninK\functions_for_calc');
addpath('C:\Matlabprograms\PetrovninK\functions');

filename=todayfilename;
MeasName = ['IQ_1pump_Digitizer_w_Triggering_DCsweep',filename];
TempDirPath=['D:\Kirill\' DeviceName '\' datestr(date) '\' MeasName '\'];
mkdir(TempDirPath);
mkdir([TempDirPath 'Data\']);
mkdir([TempDirPath 'Figures\']);
addpath(TempDirPath);
CurrentFileName=mfilename('fullpath');



copyfile([CurrentFileName,'.m'],[TempDirPath,MeasName,'.m']);

%  pumpState=1;
 
%% Set parameters and initalize device
atten = 0;%Gain of pump amplifier
pumpRaw = [linspace(0.35,0.48,5)]-atten;%linspace(-20,-2,10) [linspace(-41,-6,18) linspace(-5,7,21) linspace(7.2,16,25) linspace(16.5,20.5,9)]-atten;%[linspace(-20,-2,10) linspace(0,15,16)]-atten;%%linspace(-41,-6,18)  linspace(-2,18,21)-atten;[linspace(-41,-6,36) linspace(-5,5,21) linspace(5.2,10,50) linspace(10.5,20.5,21)]-atten;%linspace(-70,-25,11) [linspace(-70,-25,11) linspace(-24,20,21)];%-70;%[linspace(-70,-25,21) linspace(-24,20,45)]%[-70 -55 -45 -30 -25 -20 linspace(-15,25,41)];
pump=pumpCalibration(pumpRaw);
% phase_offs=0:10:360;
PumpGenFreqCent = 12.084e9;%12.0778e9;%11.995e9;
MeasFreq = 14e6;% 6052499930.9;%PumpGenFreqCent/2-6.052533066e9+PumpGenFreqCent/2;%PumpGenFreqCent/2   
RBW = 2e6; 
SGCh2 = 3;
SGCh3 = 1; %LO oscillator channel
n_ave = 1;
GAIN = db2pow(87.28); %Quadrature gain 91.6139dB +- 0.15116dB%Quadrature noise temperature 5.1895K coth(hw/2kT)=35.9
h_p = 6.62607004e-34;
VoltageGeneratorV = 0.7405;%+linspace(-0.02,0.02,51);
%(2^7);

pump_f_shift1=0;

tot_phase=0/pi*180;

rotationUseFlag=1;

% CONFIGURE ADQ14 digitizer
boardid = 1;
samplerate = 50e6;
nofsamples = 1e6;
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
         fprintf(VoltageGenerator,['SOUR' num2str(DCCh) ':VOLT:OFFS ' num2str(VoltageGeneratorV(Voltage_i))]);
        if pumpState==1
            fprintf(TriggerGenerator,['SOUR1:VOLT ' num2str(pumpRaw(ii))]);
%             fprintf(TriggerGenerator,['SOUR2:VOLT ' num2str(pump(ii))]);
%             fprintf(TriggerGenerator,['SOUR1:BURS:PHAS ' num2str(phase_offs(phase_i)) ]);
            fprintf(TriggerGenerator,'OUTP1 ON');
%             fprintf(TriggerGenerator,'OUTP2 OFF');
        
        else
%             fprintf(TriggerGenerator,'SOUR1:VOLT 0.001');
%             fprintf(TriggerGenerator,'SOUR2:VOLT 0.001');
           
            fprintf(TriggerGenerator,'OUTP1 OFF');
%             fprintf(TriggerGenerator,'OUTP2 OFF');
        end


        
        for n_a=1:n_ave
            
           
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
                             P_init(ii,Voltage_i)=varI_ave_init(ii,Voltage_i)+varQ_ave_init(ii,Voltage_i);
            %                  A1=sqrt(samplerate/(RBW)).*IQfiltering(IQ_cap,samplerate,0,1,[-RBW/2 RBW/2]+MeasFreq).';%.*exp(sum(IQphi(:))/2);[-RBW/2 -fcut]+MeasFreq ,[fcut RBW/2]+MeasFreq  
%                              A1=(B1(:,1)+B1(:,2))./sqrt(2);
%                              varI_init(ii,n_a)=var(real(A1));%Variance calculations  [-RBW/2 RBW/2]-MeasFreq ,[-RBW/2 RBW/2]+MeasFreq 
%                              varQ_init(ii,n_a)=var(imag(A1));
%                              varI_ave_init(ii)=mean(varI_init(ii,1:n_a),2);
%                              varQ_ave_init(ii)=mean(varQ_init(ii,1:n_a),2);
            %                  A0=IQ_cap;


            elseif pumpState==1
            [B1]=IQfiltering(IQ_cap,samplerate,0,n_modes, [-RBW/2 RBW/2]);%Mode filtering
            B1=B1.';
           

             

            
                L = [real((B1(:,1))), imag((B1(:,1)))];

                CovMat(:,:,n_a,ii,Voltage_i)=cov(L);
                %                  CovMat1(:,:,n_a,ii)=cov(L1);
                CovMat_ave(:,:,ii,Voltage_i)=mean(CovMat(:,:,1:n_a,ii,Voltage_i),3);
                CovMat_calc_ave(:,:,ii,Voltage_i)=4*(CovMat_ave(:,:,ii,Voltage_i)-CovMat_init_ave(:,:,ii,Voltage_i))+eye(2);%Normalization and subtraction of Preamp noise (Calibration matrix)
               
             
            AAA1 = [-1:0.01:1]*max(abs(B1));
            [MM1 CC1] = hist3([real(B1) imag(B1)], {AAA1 AAA1});
            varI(ii,n_a,Voltage_i)=var(real(B1));
            varQ(ii,n_a,Voltage_i)=var(imag(B1));
            varI_ave(ii,Voltage_i)=mean(varI(ii,1:n_a,Voltage_i),2);
            varQ_ave(ii,Voltage_i)=mean(varQ(ii,1:n_a,Voltage_i),2);
            P(ii,Voltage_i)=varI_ave(ii,Voltage_i)+varQ_ave(ii,Voltage_i);
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

         save([TempDirPath 'Data\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '.mat'])
        if pumpState==1
         savefig(hhh,[TempDirPath 'Figures\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '.fig'])
        end
%          savefig(hhh2,[TempDirPath 'Figures\' MeasName,'pumpState_' num2str(pumpState) '_' num2str(pump(ii)+atten) '_'  num2str(VoltageGeneratorV(Voltage_i)) '_COVMAT.fig'])
        %%
        if pumpState==1
%         figure (114)
%         clf;
%         % hold on
%          plot(pump(1:ii),[EIG12(1:ii,Voltage_i)' EIG23(1:ii,phase_i)' EIG13(1:ii,phase_i)' ],'LineWidth',3)
%   
% %     xlim([0 0.5])
%     hold on
%     plot(xlim,[1 1],'--','Color','g','LineWidth',3);
%       legend('EIG1-23','EIG2-13','EIG3-12','Quantum threshold');
%         figure(1144)
%         clf;
%         plot((pump(1:ii)+atten), 10*log10(4*[(varI_ave(1:ii).'-varI_ave_init(1:ii).')+1/4, (varQ_ave(1:ii).'-varQ_ave_init(1:ii).')+1/4]));
        % figure (1115)
        % plot((pump(1:ii)+atten),pur(1:ii,1),'r')
        
%         legend('EIG 1-23','EIG 2-13','EIG 3-12');
%         H = -1:0.001:1; %optimising entanglement criterion
%         G = -1:0.001:1;
%         [h,g] = meshgrid(H,G);
%         for iter=0:1:3
%         CM=circshift(CovMat_calc_ave(:,:,ii,phase_i),2*[-iter -iter]);
% %         
% %         S1(ii,iter+1) = min(min((CM(1,1)+(h.^2).*(CM(3,3)+CM(5,5)-2.*CM(3,5))+...
% %         2*h.*(CM(1,3)-CM(1,5))+CM(2,2)+(g.^2).*(CM(4,4)+CM(6,6)+...
% %         +2.*CM(4,6))+2*g.*(CM(2,4)+CM(2,6)))./(1+abs(h.*g-h.*g))))/2;
% %     
% %         S2(ii,iter+1) = min(min((CM(1,1)+(h.^2).*(CM(3,3)+CM(5,5)+2.*CM(3,5))+...
% %         2*h.*(CM(1,3)+CM(1,5))+CM(2,2)+(g.^2).*(CM(4,4)+CM(6,6)+...
% %         -2.*CM(4,6))+2*g.*(CM(2,4)-CM(2,6)))./(1+abs(h.*g-h.*g))))/2;
% %         
% %         S3(ii,iter+1) = min(min((CM(1,1)-(h.^2).*(CM(3,3)+CM(5,5)-2.*CM(3,5))+...
% %         2*h.*(CM(1,3)-CM(1,5))+CM(2,2)+(g.^2).*(CM(4,4)+CM(6,6)+...
% %         -2.*CM(4,6))+2*g.*(CM(2,4)-CM(2,6)))./(1+abs(+h.*g+h.*g))))/2;
% %     
%         S4(ii,iter+1) = min(min((CM(1,1)+(h.^2).*(CM(3,3)+CM(5,5)+2.*CM(3,5))+...
%         2*h.*(CM(1,3)+CM(1,5))+CM(2,2)+(g.^2).*(CM(4,4)+CM(6,6)+...
%         +2.*CM(4,6))+2*g.*(CM(2,4)+CM(2,6)))./(1+abs(h.*g+h.*g))))/2;
% %         S(ii,iter+1) = min(min((CM(1,1)+(h.^2).*(CM(3,3)+CM(5,5)+2.*CM(3,5))+...
% %              2*h.*(CM(1,3)-CM(1,5))+CM(2,2)+(g.^2).*(CM(4,4)+CM(6,6)+...
% %              +2.*CM(4,6))+2*g.*(CM(2,4)+CM(2,6)))./(1+abs(h.*g+h.*g))));
% %         S(ii,iter+1)=min({ S1(ii,iter+1), S2(ii,iter+1), S3(ii,iter+1), S4(ii,iter+1)});
%         end
%         S(ii,phase_i)=min(S4(ii,:),[],2);
         %2*C(2,6))+2*g.*(C(2,4)+C(4,6))


%         hhh6=figure (118);
%          clf;
% 
%         hold on
%         plot(pump(1:ii)+atten,S(1:ii),'Linewidth',3);
%         title('S-criteria');
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
        if (ii>1)&&(Voltage_i>1)
            hhh6=figure (118);
             clf;

            % hold on
            s=surf(pump(1:ii),VoltageGeneratorV(1:Voltage_i),0.5+(P(1:ii,1:Voltage_i).'-P_init(1:ii,1:Voltage_i).'));
            set(s,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
            shading interp
            legend('N-number')
            colorbar
            % title('S-criteria');\
            hhh7=figure (121);
             clf;
            s1=subplot(1,2,1)
            ss1=surf(pump(1:ii),VoltageGeneratorV(1:Voltage_i),0.25+varI_ave(1:ii,1:Voltage_i).'-varI_ave_init(1:ii,1:Voltage_i).');
            set(ss1,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
            shading interp
            legend('Variance I')
            colorbar
            s2=subplot(1,2,2)
            ss2=surf(pump(1:ii),VoltageGeneratorV(1:Voltage_i),0.25+varQ_ave(1:ii,1:Voltage_i).'-varQ_ave_init(1:ii,1:Voltage_i).');
            set(ss2,'EdgeColor','none')
%             zlim([0 2]);
%             xlim([0 .5]);
%             caxis([0 2]);
            shading interp
            legend('Variance Q')
            colorbar
     
            
        end

        
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



    
    
    