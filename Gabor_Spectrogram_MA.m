function [SpecLFP, T, F, P]=Gabor_Spectrogram_MA(LFP, fs)
%fs= sampling rate, usually 1000 Hz

% Gabor spectrograms
L       = 10;
NFFT    = 2^L;
tw      = -NFFT/2+1:NFFT/2;
sigma   = .2;%[sec]
sigSamp = sigma*fs;
w       = sqrt(sqrt(2)/sigSamp)*exp(-pi*tw.*tw/sigSamp/sigSamp);
overlap = NFFT-1;
[SpecLFP, F, T, P]=spectrogram(LFP,w,overlap,NFFT,fs); %deriv gaussian windowed spectrogram
figure;imagesc(T,F,abs(SpecLFP));set(gca,'YDir','normal');
ylim([0 150])
set(gca,'YTick',[0:10:150])
grid on
xlabel('time (sec)');ylabel('Hz')
title('Gabor spectrogram CA1 LFP \theta')


end