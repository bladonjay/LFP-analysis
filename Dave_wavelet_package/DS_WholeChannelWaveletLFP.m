function [] = DS_WholeChannelWaveletLFP(filebase,channel,NumEEGChannels,EEGSR)
% NP_WholeChannelWaveletCSD(filebase,channel)
% calculates the wavelet transform for the whole channel, saves it to disk
% also caclulates the mean and std of each frequency band so they can be saved independently
% the wavelet parameters (non-adjustable) are lowfreq: 2, highfreq:300, numbands: 64 (actually 65)
% input arguments: filebase is the base name of the session, channel is the
% channel # to make the wavelet on, NumEEGChannels is the # of channels in
% the eeg file, and EEGSR is the .eeg sampling rate
% there are no outputs since the wavelet gets saved to disk.

% step one: load the channel:
% i think this is just the voltages
v = readmulti([filebase,'.eeg'],NumEEGChannels,channel);
EEGlength = length(v);



% Some key parameters
HighFreq = 300;
LowFreq = 2;
NumBands = 64;
WaveletChunkLength = 100000; 
NumBadSamples = 858; % determined for the (300,2,64) for the above parameters

StartSample = 1;
EndSample = WaveletChunkLength;
% # of samples we need to overlap by is the timepoint at which the coi hits the lowest frequency, plus one


% open up the file
fid = fopen([filebase,'_DSlfpWaveletCH',int2str(channel)],'w');
fidph = fopen([filebase,'_DSlfpWaveletPhaseCH',int2str(channel)],'w');
samples_written = [];

% calculate the wavelet piecemeal, save power and phase to files
while (EndSample < EEGlength)
  display([num2str(EndSample/EEGlength*100),'% of the channel has been processed']); 
  % v is the signal, and wavelet idx is the epoch youre using
  wavelet_idx = StartSample:EndSample;
  [tempwave,period,scale,coi,phases] = DS_RawWavelet(EEGSR,v(wavelet_idx),HighFreq,LowFreq,NumBands);
  
  if (StartSample == 1)
    FirstWriteIdx = 1;
  else
    FirstWriteIdx = NumBadSamples + 1;
  end
  
  LastWriteIdx = WaveletChunkLength-(NumBadSamples +1);
  SamplesToBeWritten = wavelet_idx(FirstWriteIdx:LastWriteIdx);
  
  for i = FirstWriteIdx:LastWriteIdx
    fwrite(fid,(single(tempwave(:,i))),'single');
    fwrite(fidph,(single(phases(:,i))),'single');
  end
  samples_written = [samples_written,SamplesToBeWritten];
  
  StartSample = StartSample + WaveletChunkLength-(2*NumBadSamples+1);
  EndSample = StartSample+WaveletChunkLength - 1;
  
end

% do the last bit
EndSample = EEGlength;
wavelet_idx = StartSample:EndSample;
[tempwave,period,scale,coi,phases] = DS_RawWavelet(EEGSR,v(wavelet_idx),HighFreq,LowFreq,NumBands);
FirstWriteIdx = NumBadSamples + 1;
LastWriteIdx = length(wavelet_idx);
SamplesToBeWritten = wavelet_idx(FirstWriteIdx:LastWriteIdx);
for i = FirstWriteIdx:LastWriteIdx
  fwrite(fid,single(tempwave(:,i)),'single');
  fwrite(fidph,single(phases(:,i)),'single');
end
samples_written = [samples_written,SamplesToBeWritten]; 

% save some important info
savestr = ['save ',filebase,'_DSlfpWaveletInfo.mat HighFreq LowFreq period scale EEGlength'];
eval(savestr);



f = 1./period;
save WaveletFrequencies.mat f -ASCII

fclose(fid);
fclose(fidph);

