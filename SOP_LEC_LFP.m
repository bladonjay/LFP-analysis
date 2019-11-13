%% LEC LFP Intro

load('K:\LEC Projects\Scooby Sessions short.mat')

%% I only have lfp data for one session, but lets see if we can get this to work

% this is a ton of LFP data...
LFP=ReconstituteLFP(Session(2).LFP.contvars);

%LFP=CleanLFPquick(

%%
for k=1:2:10
%k=28;
data=LFP(k).data; ts=LFP(k).ts;

flag=Session(2).board_matrix(:,4);
before=4; after=16;
fs=.001;
trialLFP=[];


for i=1:40
    startind=find(ts>(flag(i)-before),1);
    endind=startind+after/fs; % ten seconds
    trialts(:,i)=ts(startind:endind);
    trialLFP(:,i)=data(startind:endind);
end
figure;

subplot(2,1,1);
imagesc(trialLFP'); colormap hsv
subplot(2,1,2);
plot([0.001:.001:length(mean(trialLFP'))*.001]-before, mean(trialLFP')); xlim([-before length(mean(trialLFP'))*.001-before]);



end
 %%
params.Fs=1000;
params.fpass=[1 150];
params.tapers=[3 5];
params.err=[2 .05];
params.pad=2;

%
[S,t,f,Serr]=mtspecgramc(trialLFP,[0.5 0.025],params);
%%
for i=1:5
[wave(:,:,i),period,scale]=DS_wavelet(trialLFP(:,i),.001,2,.1,.005);
figure; imagesc(real(wave(:,:,i)));
end

%%
for k=1:5
spect=real(S(:,:,k))';
for i=1:length(spect(:,1));
    spect(i,:)=spect(i,:)*f(i);
    %spect(i,:)=zscore(spect(i,:));
end

figure;
imagesc(t-before,f,spect);
 set(gca,'YDir','Normal');
end
%%
spect=(real(mean(S,3))');
for i=1:length(spect(:,1));
    %spect(i,:)=spect(i,:)*f(i);
    spect(i,:)=zscore(spect(i,:));
end

imagesc(t-before,f,spect); set(gca,'YDir','Normal');


%%
% lets see what each band does individually
