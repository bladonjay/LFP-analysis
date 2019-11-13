function [phases] = DS_LoadWaveletPhaseBandLFP(filebase,channel,frequencyidx)
% function [phases] = DS_LoadWaveletPhaseBandLFP(filebase,channel,frequencyidx)
% loads phases for wavelet level frequencyidx
NumLevels = 65;

load([filebase,'_DSlfpWaveletInfo.mat']);

display(['loading phases for ',num2str(1./period(frequencyidx)),' Hz frequency band']);

  fid = fopen([filebase,'_DSlfpWaveletPhaseCH',int2str(channel)]);
  fseek(fid,(frequencyidx-1)*4,-1);
  
  EEGlength = NP_GetEEGLength(filebase);
  
  phases = zeros(EEGlength,1);
  
  phases = fread(fid,EEGlength,'single',64*4);


 fclose(fid);

