function [wave,period,t] = DS_LoadWaveletPhaseChunkLFP(filebase,channel,StartSample,EndSample,EEGSR,toplot)
% [wave,period,t] = NP_LoadWaveletPhaseChunkLFP(filebase,channel,StartTime,EndTime,toplot)
% this function loads the wavelet data for the recording segment specified 
% this function assumes that WholeChannelWaveletCSD has been called for the specified channel
% the wavelet, along with the associated trace, is plotted if toplot = 1
if (nargin < 6)
  toplot = 0;
end

NumLevels = 65;

NumSamples = length(StartSample:EndSample);


fid = fopen([filebase,'_DSlfpWaveletPhaseCH',int2str(channel)]);

fseek(fid,(StartSample-1)*4*NumLevels,-1);

load([filebase,'_DSlfpWaveletInfo.mat']);

for i = 1:NumSamples
  wave(1:65,i) = fread(fid,NumLevels,'single');
end
wave = wrapTo2Pi(wave);


t = StartSample:EndSample;

if (toplot == 1)
  figure;contourf((StartSample:EndSample)/EEGSR,(1./period),wave,10);shading flat;caxis([0 2*pi]);colormap hsv%colorbar;
  
end
fclose(fid);