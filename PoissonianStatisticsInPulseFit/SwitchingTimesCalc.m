 %% Load data
clear all
% load(['E:\Kirill\QWJPA_v2_2\09-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106091948\ProbeCharacter2us.mat']);
% load('E:\Kirill\QWJPA_v2_2\08-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106080307\PoissonFitDataPulseCharact.mat');
% load('E:\Kirill\QWJPA_v2_2\10-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106102039\ProbeCharacterization0.5us.mat');
% load('E:\Kirill\QWJPA_v2_2\10-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106102331\ProbeCharacterization0.25us.mat');

load('F:\Kirill\QWJPA_v2_2\11-Jun-2021\Probe_Detection_shPulse_PhotonNumberSweep202106111331\ProbeCharacter1us.mat');
%% find out indeces with value>Threshold
% ThresVal=ThresMean;
DataArray=reshape(abs(IQ_raw(1,:)+1i.*IQ_raw(2,:)),[length(IQ_raw(1,:))/samplerate/pulseLength, 500]);
%%
clearvars t1
t1=zeros(1,size(DataArray,1));
ThresVal=4e-2;
for ii=1:size(DataArray,1)
    if isempty(find(DataArray(ii,:)>ThresVal,1))==0
    t1(ii)=(find(DataArray(ii,:)>ThresVal,1)-1)./samplerate;
    end
end


hist(t1,500)

