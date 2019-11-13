function [wave,period,t] = DS_LoadWaveletChunkLFP(filebase,channel,StartSample,EndSample,NormName,EEGSR,toplot)
% [wave,period,t] = NP_LoadWaveletChunk(filebase,channel,StartTime,EndTime,toplot)
% this function loads the wavelet data for the recording segment specified (CSD)
% this function assumes that WholeChannelWaveletCSD has been called for the specified channel
% the wavelet, along with the associated trace, is plotted if toplot = 1
if (nargin < 7)
  toplot = 0;
end

NumLevels = 65;


NumSamples = length(StartSample:EndSample);

fid = fopen([filebase,'_DSlfpWaveletCH',int2str(channel)]);

fseek(fid,(StartSample-1)*4*NumLevels,-1);

load([filebase,'_DSlfpWaveletstatsCH',int2str(channel),'_',NormName,'.mat']);
load([filebase,'_DSlfpWaveletInfo.mat']);

for i = 1:NumSamples
  wave(1:65,i) = fread(fid,NumLevels,'single');
  
end

for i = 1:NumLevels
  wave(i,:) = (wave(i,:)-BandMean(i))./BandStd(i);
end

t = StartSample:EndSample;

if (toplot == 1)
  figure;contourf((StartSample:EndSample)/EEGSR,(1./period),wave,40);shading flat;caxis([0 2]);xlabel('time (sec)');ylabel('frequency (Hz)');colorbar;
  caxis([0 2]);
end
fclose(fid);