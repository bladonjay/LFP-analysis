%Ripple/Replay SOP

%% look for spikes in multiunit activity on the order of 100 ms

load('K:\LEC Projects\Scooby Sessions short.mat')


%% first cat all spikes

ses=Session(2);
keep=[1 2 3 4 5 6 8 9 11 12 13 16 17 19 20 21];
mua=[];
for i=1:length(keep)
    spktemp=ses.unitdata.units(keep(i)).ts;
    spktemp(:,2)=i;
    mua=[mua;spktemp];
    fprintf('unit %d done \n',i);
end
[~,ia]=sort(mua(:,1));

mua=mua(ia,:);

%%
bins=0:.02:max(mua);

[y,x]=hist(mua(:,1),bins);
%%
[lfpdata,lfptime]=CleanLFPquick(ses.LFP.contvars{5});

band=[80 120];
[gamphase,gamamp] = GetLFPBand(lfpdata,lfptime,band);
band=[6 9];
[thetaphase,thetaamp] = GetLFPBand(lfpdata,lfptime,band);


%%
plot(x,smooth(y)-mean(y));
hold on
plot(lfptime,lfpdata,'r');

%%

bins=0:lfptime(2)-lfptime(1):max(mua);
[y,x]=hist(mua(:,1),bins');

plot(x,zscore(y));
hold on
plot(lfptime,zscore(lfpdata),'r');

%% lets plot the mua's and do a mua centered spectrogram;

bins=0:.01:max(bins);
[counts,bins]=hist(mua(:,1),bins);
Scounts=SmoothMat2(counts,[6,1]);
figure;

plot(bins,counts)


hold on

plot(lfptime,zscore(SmoothMat2(gamamp,[1,6])),'r');
plot(lfptime,zscore(thetaamp),'g');
legend('mua','gamma','theta');
%%
% lets see if theres a trial structure to the mua
muamat=[];
timespan=16/0.02;
for i=1:length(ses.flags.treadmill)
    ia=find(bins>ses.flags.treadmill(i)-2,1,'first');
    muamat(i,:)=Scounts(ia:ia+timespan);
end
figure;
subplot(2,2,1);
imagesc(-2:.02:14,1:85,zscore(muamat,1,2));
subplot(2,2,2);
imagesc(-2:.02:14,1:85,zscore(muamat,1,2)>2);
subplot(2,2,3);
imagesc(-2:.02:14,1:85,zscore(muamat,1,2)>3);
subplot(2,2,4);
plot(-2:.02:14,mean(zscore(muamat,1,2)>2));
%%
%