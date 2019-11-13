function [wave,period,scale,coi,phases] = DS_RawWavelet(EEGSR,EEGchunk,UpperLimitHz,LowerLimitHz,NumScales)
% [tempwave,period,scale,coi,phases] = DS_RawWavelet(EEGSR,EEGchunk,UpperLimitHz,LowerLimitHz,NumScales)
% this function is a wrapper for the Torrence and Compo wavelet.m function,
% tempwave is power, period is the period of each level (i.e., frequency
% band), coi is the "cone of influence", and phases are the wavelet phases

dt = 1/EEGSR; %  time per sample of the eeg file (dt)

s0 = 1 / UpperLimitHz; %(e.g. period of fastest freq)

j1 = NumScales; % # of scales: frequency resolution (usually 1/4 per hz)

LowerLimitSeconds = 1 / LowerLimitHz; % usually 2 Hz
dj = log2((LowerLimitSeconds/s0))/j1;

[rawwave,period,scale,coi] = DS_wavelet(EEGchunk,dt,1,dj,s0,j1);

wave = log(abs(rawwave).^2);

phases = atan2(imag(rawwave),real(rawwave));
