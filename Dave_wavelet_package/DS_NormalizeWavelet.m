function [] = DS_NormalizeWavelet(filebase,channel,epochs,NormName)
%DS_NormalizeWavelet 
%   filebase is the base name of the session
% channel is the channel to normalize
% epochs is an nX2 matrix, the first column being starts of epochs to use
% for normalization,the second column being the ends of these epochs.  epochs must be in 
% EEG samples
% NormName is a string to specifty which normalization epochs were used.
% for example, 'nonREM'


NumBands = 65; % CHANGE THIS ONLY IF YOU ABSOLUTELY KNOW WHAT YOU:RE DOING

load([filebase,'_DSlfpWaveletInfo.mat']); % gets HighFreq LowFreq period scale EEGlength'
GoodSamples = [];
for i = 1:size(epochs,1)
  GoodSamples = [GoodSamples,epochs(i,1):epochs(i,2)];
end
    

display('Calculating mean and std of each frequency band');
for i = 1:NumBands
  display(['calculating mean and std of the ',num2str(1./period(i)),' Hz level of the wavelet'])
  tempfreq = DS_LoadWaveletFrequencyBandLFP(filebase,channel,i);
  BandMean(i) = mean(tempfreq(GoodSamples));
  BandStd(i) = std(tempfreq(GoodSamples));
  display(['mean = ',num2str(BandMean(i))]);
  display(['std = ',num2str(BandStd(i))]);
end
savestr = ['save ',filebase,'_DSlfpWaveletstatsCH',int2str(channel),'_',NormName,'.mat BandMean BandStd'];
eval(savestr);
display('Done!!!');

figure(479);plot(1./period,BandMean);hold on;errorbar(1./period,BandMean,BandStd);axis tight;xlabel('frequency Hz');ylabel('mean power, errorbars are std');
set(gca,'XScale','log');
end


