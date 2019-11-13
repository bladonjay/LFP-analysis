%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TEST THE FUNCTIONS OUT ON KNOWN DATA %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lets test these functions out
% first lets make a wave that switches freq half way thru

xt=.001:.001:2;
% a 10 Hz oscillation
noise=rand(1,length(xt));
x=sin(2*pi*xt*40)+noise;

% then a 5 Hz oscillation
x(1001:2000)=sin(2*pi*xt(1:1000)*25)+noise(1001:2000);
plot(xt,x);

%
dt=.001;

%%%%% spectrogram %%%%%%%
% x, window, noverlap, f, fs
[s]=spectrogram(x,100,0,[2:.5:100],1000);
imagesc(abs(s));


%%%%% periodogram %%%%%%%
% now lets look how to pick up both, single taper
[pxx,f]=periodogram(x,hamming(length(x)),[2:100],1000);
figure; plot(f,pxx); title('periodogram');


%%%% multi taper periodogram%%%%
% inputs: data, T*sampfreq/2, Freqs to look at, samp freq
[pxx,f] = pmtm(x,5,[1:100],1/dt);
figure; plot(f,pxx.*log(f)); title('multitaper');

%%%% wavelet transform %%%%%%%
figure;
% data, dt, pad, dj(df), So(period of fastest oscillation 200 hz)
[wave{1},period,scale]=DS_wavelet(x,.001,2,.1,.005,55);
imagesc(abs(wave{1}).^2);
period2=1./period;
set(gca,'YTick',1:2:length(period2));
set(gca,'YTickLabel',period2([1:2:end]))


%%%%%%% PlotPowerSpectrum %%%%%%
[faxis,Sxx]=PlotPowerSpectrum(x,xt);
figure; plot(faxis,Sxx.*log(f)); title('power spect');


%%% Chronux MT spectrum %%%%%%
params.Fs=1000;
[S,F]= mtspectrumc(x,params);
figure; plot(F,S);




% scale to frequency
%f=scal2frq(scale,'cmor',.001);
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Lets just observe some data and get it organized %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('C:\Users\jhbladon\Desktop\Rapid Analysis\LEC 9-9-14.mat')
load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC 2-24-15.mat')
lfpdir='C:\Users\jhbladon\Desktop\Rapid Analysis\LFP data';

lfpfiles=dir(lfpdir);
dash=WhichDash;



sessionnum=8; % this is the session


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%look at data%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i think many of these thetas are clocking in at 8 hz for MEC

session=LEC(4);
freqband=[10 20];
data=session.lfp(:,1); ts= session.lfp(:,2);
[sessionlfp]=PeekSessionLFP(session.bf, data, ts,.001,freqband);

% I think LEC hits a big beta peak

[pxx,f]=periodogram(data,hamming(length(data)),[2:100],1000);
figure; plot(f,pxx); title('periodogram');


%%
% lets specify the time epoch, hopefully stitching wont add too much noise
for i=1:length(session.bf.sample_start)
    zerots=find(session.lfp(:,2)<session.bf.sample_start(i),1,'last');
    start=zerots-10000; finish=zerots+10000;
    start2=zerots-1000; finish2=zerots+1000;
    
    %  wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
    thisdata=session.lfp(start:finish,1);
    thisdata=thisdata-mean(thisdata);
    [wave,period,scale]=DS_wavelet(thisdata,.001,2,.1,.005,55);
    
    if i==1
        avgwave=abs(wave);
        data2=session.lfp(start2:finish2,1);
        data2=data2-mean(data2);
        ts2=session.lfp(start2:finish2,2);
    else
        tempdata2=session.lfp(start2:finish2,1);
        tempdata2=tempdata2-mean(tempdata2);
        data2=[data2;tempdata2];
        ts2=[ts2 session.lfp(start2:finish2,2)];
        avgwave=avgwave+abs(wave);
    end
    fprintf('%d done \n', i);
end
data2=data2-mean(data2);
wave=wave./30;
[pxx,f]=periodogram(data2,hamming(length(data2)),[5:100],1000);
figure; plot(f,pxx); title('periodogram');


figure; imagesc(abs(avgwave).^2);
period2=1./period;
set(gca,'YTick',1:2:length(period2));
set(gca,'YTickLabel',round(period2([1:2:end])))


params.Fs=1000;
[S,F]= mtspectrumc(data2,params);
figure; plot(F,S);
%%
% use only a little data
datapeek=data(1:1000000); tspeek=ts(1:1000000);

[pxx,f]=periodogram(data,hamming(length(data)),[1:100],1000);
% or many tapers
%[pxx,f]=pmtm(datapeek,3,[1:.5:100],1000);
figure; semilogy(f,smooth(pxx,2)); title('periodogram');
%% Maybe some cross frequency coupling?
% this will take a huge amount of time, maybe parse down data for only when
% these two freq bands are prevalent? say take the hilbert and only use
% when amp is above lets say 1.5 z
for i=1:30
    for j=1:10
        [Phaseamp(i,j),Metric(i,j)]=CrossFreqCouple(data,ts,[i+5 i+10], [j*5+40 j*5+45],0);
    end
    fprintf([num2str(i) ' \n']);
end

%%
% grab the wideband data for each trial:
%trialdur=mean(diff(session.bf.trial));
flag=session.bf.sample_start;
firstbin=-2;
lastbin=2;
% this is as big as each trial can get

for q=1:length(flag)    
    
    % get the index of the flag
    index=find(ts>flag(q),1,'first');

    %build the image
    trialdata=data(index+firstbin*1000:index+lastbin*1000);
    
    
    [wave,period,scale]=DS_wavelet(trialdata,.001,2,.1,.005,55);

    if q==1
        avgwave=abs(wave);
    else
        avgwave=avgwave+abs(wave);
    end
    
    
    fprintf('trial %d finished \n',q);
end
%wave=wave./30;
figure;
imagesc(abs(avgwave).^2);
period2=1./period;
set(gca,'YTick',1:2:length(period2));
set(gca,'YTickLabel',round(period2([1:2:end])))
set(gca,'XTick',0:500:4000);
set(gca,'XTickLabel',-2:.5:2);

%%
dt=.001; 
thetafreq=[6 9];  thetapass=thetafreq.*(2*dt);
gammafreq=[40 60];
[a,b]=butter(3,thetapass);
[c,d]=butter(3,gammafreq.*(2*dt));
% lets plot our theta amplitude and our gamma amplitude;

clear thetaamp gammaamp
for i=1:length(flag)
    
    theta=filtfilt(a,b,(trialdata(i,:)));
    thetaparts=hilbert(theta);
    thetaamp(i,:)=decimate(abs(thetaparts),50);
    gamma=filtfilt(c,d,trialdata(i,:));
    gammaparts=hilbert(gamma);
    gammaamp(i,:)=decimate(abs(gammaparts),50);


end

scales=5:.5:100;
wname='cmor2-0.5';
freq=scal2frq(scales,wname,.001);
time=trialts(1,:);
for i=1:length(flag);
    % build a spectrogram too
    %[spect,period,scale]=DS_wavelet(trialdata(i,:),.001,2,.5,.01);
    spect=cwt(trialdata(i,:),scales,wname);
    %spect=spectrogram(trialdata(i,:),300,0,[1:2:100],1000);
    if i==1
        avgspect=abs(spect);
    else
    avgspect=avgspect+abs(spect);
    fprintf('trial %d finished \n',i);
    end
end
period=1./period;

avgspect=avgspect./89;
imagesc(time,scales,avgspect);
set(gca,'YTick', scales(1:2:end));
set(gca,'YTickLabel',round(freq(1:2:end)));
timesteps=0.004:0.004:length(gammaamp)*.004;

%%
for i=1:length(flag)
    for i=1:
%%
