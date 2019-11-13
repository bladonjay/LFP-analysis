

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TEST THE FUNCTIONS OUT ON KNOWN DATA %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lets test these functions out
% first lets make a wave that switches freq half way thru

xt=.001:.001:2;
% some noise
noise=rand(1,length(xt));
% 40 Hz
x=sin(2*pi*xt*40)+noise;
% then a 25 Hz oscillation
x(1001:2000)=sin(2*pi*xt(1:1000)*25)+noise(1001:2000);

plot(xt,x);

%
dt=.001;

%%%%% spectrogram %%%%%%%
% x, window, noverlap, f, fs
[s,f,t]=spectrogram(x,100,0,[2:.5:100],1000);
imagesc(t,f,abs(s));


%%%%% periodogram %%%%%%%
% now lets look how to pick up both
[pxx,f]=periodogram(x,hamming(length(x)),[2:100],1000);
figure; plot(f,pxx); title('periodogram');

%%%% multi taper periodogram%%%%
% inputs: data, T*sampfreq/2, Freqs to look at, samp freq
[pxx,f] = pmtm(x,5,[1:100],1/dt);
figure; plot(f,pxx.*log(f)); title('multitaper');

%%%% wavelet transform %%%%%%%
figure;
% data, dt, pad, dj(df), So(period of fastest oscillation 200 hz)
[wave{1},period,scale]=DS_wavelet(x,.001,2,.1,.005);
imagesc(abs(wave{1}).^2);
period2=1./period;


%%%%%%% PlotPowerSpectrum %%%%%%
[faxis,Sxx]=PlotPowerSpectrum(x,xt);
figure; plot(faxis,Sxx.*log(f)'); title('multitaper');


% scale to frequency
%f=scal2frq(scale,'cmor',.001);

%%% Chronux MT spectrum %%%%%%
params.Fs=1000;
params.fpass=[1 90];
params.tapers=[2 3];
params.err=[2 .05];
params.pad=2;

[S,t,f,Serr]=mtspecgramc(x,[0.5 0.025],params);

%% maybe using a morlet continuous wavelet from matlab? 
dt=.001;
scales=1:1:100;
%t=ts(1:1000000); y=data(1:1000000);
fb=2; % envelope freq
fc=0.5; % .5 Hz inside freq
wname = ['cmor' num2str(fb) '-' num2str(fc)];
%waveid='db2';
f=scal2frq(scales,wname,dt);

[X]=cwt(x,scales,wname);
subplot(1,2,1);
imagesc(t,scales,real(X))
set(gca, 'YTick', scales(1:2:end))
set(gca, 'YTickLabel', round(f(1:2:end)))
xlabel('Time [s]')
ylabel('Freq [Hz]')
subplot(1,2,2);
imagesc(t,scales,abs(X))
linkaxes;
set(gca, 'YTick', scales(1:2:end))
set(gca, 'YTickLabel', round(f(1:2:end)))
xlabel('Time [s]')
ylabel('Freq [Hz]')
colorbar

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Lets just observe some data and get it organized %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('C:\Users\jhbladon\Desktop\Rapid Analysis\LEC 9-9-14.mat')
load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC 2-24-15.mat')
lfpdir='C:\Users\jhbladon\Desktop\Rapid Analysis\LFP data';

lfpfiles=dir(lfpdir);
dash=WhichDash;

k=4;  % lfp file number

sessionnum=10; % this is the session
load([lfpdir dash lfpfiles(k).name]);
LFPstruct=session.LFP.contvars{13};
%% SET UP THE LFP TRACE SO ITS CLEAN AND HAS PROPER TIMESTAMPS


% first lets clean up our LFP data;
fragments=LFPstruct.fragmentStarts;
timefrags=LFPstruct.timestamps;
% add first timestamp
LFPstruct.time=[];
dt=1/LFPstruct.ADFrequency;
% last fragment we have to treat different because it doesnt have a stop
% time
for i=1:length(fragments)-1
    % find number of datapoints recorded
    thislength=fragments(i+1)-fragments(i);
    % now step forward till right before the next frag starts
    fragts=(dt:dt: thislength *dt) + timefrags(i);
    LFPstruct.time=[LFPstruct.time fragts];
end

% adress last fragment, unique because now we account for last timestamp:
lastlength=length(LFPstruct.data)-fragments(end)+1;
lastfrag=(dt:dt:(lastlength*dt)) + fragts(end);
LFPstruct.time=[LFPstruct.time lastfrag];

ts=LFPstruct.time'; data=LFPstruct.data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% filter the LFP %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scrub timestamps and data that are railing
% I'll use the max and min, so either i will scrub two datapoints, or all
% those that rail

% So then the data is all concatenated, and where there are borders, i can
% go back with my timestamps and remove

bottomrail=min(data);
removeind=find(data==bottomrail);
data(removeind)=[]; ts(removeind)=[];

toprail=max(data);
removeind=find(data==toprail);
data(removeind)=[]; ts(removeind)=[];
% this should yield a cleaned data trace
%%
% i think  thetas  clocking mostly in at 8 hz
freqband=[6 10];

[sessionlfp]=PeekSessionLFP(session.bf, data, ts,.001,freqband);

% myname=[MEC(i).name ' LFP ' num2str(freqband)];
% title(myname)
% 
[pxx,f]=periodogram(data,hamming(length(data)),[2:100],1000);
%% or many tapers
% 
session=Session(2);
%for i=6:16
[data,ts]=CleanLFPquick(session.LFP.contvars{14});
[pxx,f]=pmtm(data(1:10000),6,[2:100],1000);
figure; semilogy(f,pxx); title(['periodogram lfp ' num2str(i)]);
%end
%%
% grab the wideband data for each trial:
trialdur=10;
flag=session.flags.treadmill;
firstbin=-trialdur/2;
lastbin=trialdur/2;
% this is as big as each trial can get
trialdata=nan(length(flag),ceil((trialdur)*1000));
trialts=nan(length(flag),ceil((trialdur)*1000));
for q=1:length(flag)
    % unfortunately i pull in a ton of
    
    timestart=flag(q)+firstbin;
    timeend=flag(q)+lastbin;
    indices=ts>timestart & ts<timeend;
    %build the image
    trialdata(q,1:sum(indices))=data(indices);
    trialts(q,1:sum(indices))=ts(indices)-flag(q);
    
    fprintf('trial %d finished \n',q);
end

imagesc(trialdata);

%%
%%% Chronux MT spectrum %%%%%%
params.Fs=1000;
params.fpass=[5 150];
params.tapers=[2 3];
params.err=[2 .05];
params.pad=2;

[S,t,f,Serr]=mtspecgramc(trialdata(1,:),[0.5 0.025],params);
%%
thetafreq=[7 10]; gammafreq=[40 80]
[a,b]=butter(3,thetafreq.*(2*dt));
[c,d]=butter(3,gammafreq.*(2*dt));
% lets plot our theta amplitude and our gamma amplitude;

    

%%
% now lets see if we can grab some phase and amplitude info from each trial

%%
% lets look at some raw traces
ts=TETFP01(:,2);
plot(ts(1:100000),TETFP21(1:100000,1));

% howabout we average a few traces to get some regional LFP

% now lets downsample to 250 HZ (im not interested in any oscillations over
% 100
timestamps2=decimate(Timestamps,4); % take a quarter of the data
LFP2=decimate(TETFP09,4);
dec_fs=timestamps2(2)-timestamps2(1);


% lets get a multitapered spectrum
 % signal, a filter, frequency looks, and samp rate(Hz)
[pxx,f]=periodogram(LFP2,hamming(length(LFP2)),[2:100],dec_fs);
figure; plot(f,pxx); title('periodogram');

% inputs: data, T*sampfreq/2, Freqs to look at, samp freq
[pxx,f] = pmtm(LFP2,3,[1:100],dec_fs);
figure; plot(f,pxx.*log(f)); title('multitaper');

% what does the power spectrum look like otherwise?
[faxis,Sxx]=PlotPowerSpectrum(LFP2,timestamps2);
figure; plot(faxid,Sxx.*log(f)); title('multitaper');


multiwf=[TETFP01(:,1) TETFP05(:,1)];
avgwf=mean(multiwf,2);
plot(ts(1:100000),avgwf(1:100000))
% broadband theta nested gamma, there are many gammas tho, so see if you
% can ID them using the pmpm or periodogram
[Phaseamp,Metric]=CrossFreqCouple(LFP2,dec_fs,[7 11], [40 80],1);

%% lets see if we can pull a trial structure out of this

% two sessions im working on, a megaman one and a jalapeno one,
% load up the session file and the LFP file to get everyhting

data=LFP(1).data;
ts=LFP(1).ts;


dt=.001;
% now lets get a trial structure and do some math
samples=Session(1).flags.board_start;


for j=1:length(samples)
    % get starting timestamp
    thistime=samples(j);
    % find which band of LFP to grab
    startind=find(ts>thistime,1);
    endind=startind+round(4/dt);
    thisLFP(j,:)=data(startind:endind);
    thisLFPtime=ts(startind:endind)-ts(startind);
    [faxis,Sxx]=PlotPowerSpectrum(thisLFP(j,:),thisLFPtime,'hanning',0,'suppress',1);
    %[s,f,t]=spectrogram(thisLFP);
    smoothSxx(j,:)=smooth(Sxx,10);
end
%%
finalSxx(:,i)=mean(smoothSxx,2);
fprintf('bin %d done \n',i);

    %%
imagesc(real(finalSxx));




%% lets try to butterworth our data and see what the resultant LFP looks like


% dt
FS=1000;
% the freq pass
fpass=[6 12];
% filter the LFP
fpass=fpass*(2/FS);
[a,b]=butter(3,fpass);
lfp=filtfilt(a,b,data);

hlfp=hilbert(lfp);
phi=angle(hlfp);
amp=abs(hlfp);
samples=Session.flags.board_up;

for j=1:samples
    thistime=samples(j)-2;
        % find which band of LFP to grab
        startind=find(ts>thistime,1);
        endind=startind+round(4/dt);
        thisLFP(j,:)=amp(startind:endind);
        thisLFPtime=ts(startind:endind)-ts(startind);
        
end

%%
% maybe average the time by power across some epoch that I can use
% every epoch looks pretty much the same






[s]=spectrogram(F10(1:200000,1),500,0,[2:100],1000);

imagesc(abs(s));
set(gca,'ydir','normal')
% ds wavelet( data, dt, pad, chg in scales, starting freq
[wave{1},period,scale]=DS_wavelet(F10(1:10000),.001,0,.01,.01,450);
figure; imagesc(abs(wave{1}));
set(gca,'YTick',(1:50:length(period)));
freqs=1./period;
set(gca,'YTickLabel', freqs(1:50:end));

set(gca,'XTick',linspace(0,10000,10));
set(gca,'XTickLabel',-2:.4:2);

% lets try to get the average spectrogram locked to an event, like digging
% to do that take
ts=TETFP01(:,2);
data=TETFP01(:,1);
% the min is .625 and the max is the same

% I guess 1000 steps should be good, that means 500 bins on each side, 500
% is .625 and -500=-.625
data=round((data*1000)*(500/625));


% get our voltages to look goo
%%
dt=.001;
scales=1:1:100;
t=ts(1:1000000); y=data(1:1000000);
waveid='cmor2-.05';
f=scal2frq(scales,waveid,dt);

X=cwt(y,scales,waveid);
imagesc(t,scales,real(X));
set(gca,'YTick', scales(1:2:end));
set(gca,'YTickLabel', round(f(1:2:end)));
ylabel('freq');

%%
samples=Session.bf.end_first;
ts=F10_ts(1):.001:F10_ts(end); % this probably isnt quite right
data=F10(1:length(ts));
dt=.001;
finalwave=[];
for j=1:length(samples)
    % get starting timestamp
    thistime=samples(j)-2;
    % find which band of LFP to grab
    startind=find(ts>thistime,1);
    endind=startind+round(4/dt); % so 4 seconds, make sure vectors same length
    thisLFP(j,:)=data(startind:endind);
    thisLFP(j,:)=thisLFP(j,:)-mean(thisLFP(j,:));
    %thisLFPtime=LFP1ts(startind:endind)-LFP1ts(startind);
    %[wave{j},f,t]=spectrogram(thisLFP(j,:),100,0,[2:.1:100],1000);
    % data, dt, pad, df(by scale not freq), highest freq
    % period(100Hz), freqs
    [wave{j},period,scale]=DS_wavelet(thisLFP(j,:),.001,0,.01,.01,450);
    %[wave{j},period,scale,coi,phases] = DS_RawWavelet(1000,thisLFP(j,:),200,5,400);
end
freqscale=1./period;
% remember the period refers to the freq, not the scale

finalwave=abs(wave{1});
for i=2:length(wave)
    finalwave=finalwave+abs(wave{i});
end
finalwave=finalwave/length(wave);
figure; imagesc(finalwave,[0 100]);
for i=1:length(finalwave(:,1))
    finalwave(i,:)=finalwave(i,:)*1/i;
end
% wave = log(abs(rawwave).^2); % to normalize

% from what i can tell this rat clocks theta at around 8.2 Hz
% PlotPowerSpectrum
% now lets see what we can find for theta phase?

%%
ts=F10(:,2);

% dt
FS=1000;
% the freq pass
fpass=[6 12];
% filter the LFP
fpass=fpass*(2/FS);
[a,b]=butter(3,fpass);
lfp=filtfilt(a,b,data);
hlfp=hilbert(lfp);
phi=angle(hlfp);
amp=abs(hlfp);

%% wonder if i can find spike phase coherence?

%spikets=Session.units(1).ts;
% ehren used a 4th order butter filter from 6-10 HZ, but id like to see
% whether 
phases=interp1(ts,data',Session.units(1).ts,'nearest');

for i=1:length(Session.units)
    phases=interp1(ts,phi,Session.units(i).ts,'nearest');
    figure;
   rose(phases,20);
end

% another way to do it is to use the chronux package: coherencycpb (cont.
% perbin)
% basically use 1000 Hz LFP data, and snap your spikes to 1000 hz and make
% them poisson 0's and 1's at that frequency, then just feed in there and
% note the params:

%spikes=hist(
params.Fs=1/dt; % 1000
params.tapers=[5 10]; %(time bandwith product and number of tapers_
params.trialave=1; % average across trials but i think you can change this
%[C,phi,S12,S1,S2,f]=coherencypb(

% like 9/10 cells are phase locked to the rising phase of the LFP, what
% math do i use to ensure coherence





