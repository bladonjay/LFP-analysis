function [a] = DS_LoadWaveletFrequencyBandLFP(filebase,channel,frequencyidx,NormName)
% NP_LoadWaveletFrequencyBandLFP(filebase,channel,frequencyidx,NormName)
% if you want to load z-score normalized wavelets, enter the normalization
% name as NormName.  For example, 'nonREM'.  If you have not already
% calculated the normalization, do so with DS_NormalizeWavelet.m

if (nargin <= 3)
  normalized = 0;
else
  normalized = 1;
end

NumLevels = 65;

NP_NavDir(filebase);


fid = fopen([filebase,'_DSlfpWaveletCH',int2str(channel)]);
fseek(fid,(frequencyidx-1)*4,-1);

load([filebase,'_DSlfpWaveletInfo.mat']);

display(['loading ',num2str(1./period(frequencyidx)),' Hz frequency band']);

a = fread(fid,EEGlength,'single',64*4);

if (normalized == 1)
  load([filebase,'_DSlfpWaveletstatsCH',int2str(channel),'_',NormName,'.mat']);
  a = (a-BandMean(frequencyidx))./BandStd(frequencyidx);
end

fclose(fid);