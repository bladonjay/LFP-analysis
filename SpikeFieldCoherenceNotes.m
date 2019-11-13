%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Lets just observe some data and get it organized %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC LEC\LEC 5-27-15.mat')
%load('C:\Users\jhbladon\Desktop\Rapid Analysis\MEC LEC\MEC 2-24-15.mat')


%lfpdir='C:\Users\jhbladon\Desktop\Rapid Analysis\MEC LEC\LFP data';
lfpdir='F:\LFP\Mats\LEC';
%lfpdir='F:\LFP\Mats\MEC';
lfpfiles=dir(lfpdir);
dash=WhichDash;

%% This matches LFP's to each session
% this is just to match LFPs to the sessions
% now lets get the whole session theta phase
% recursively go through each session
Region=LEC;
for i=1:length(Region);
    % go through sessions
    % then find the LFP file that yo uwant
    for j=1:length(lfpfiles)
        % so match rat name and then match the date
        if any(strfind(lfpfiles(j).name, Region(i).name(1:4))) &&...
                any(strfind(lfpfiles(j).name, Region(i).name(end-6:end)))
            % if you match the lfp
            load([lfpdir dash lfpfiles(j).name]);
            [LFPdata,LFPts]=CleanLFPquick(LFPstruct);
            Region(i).lfp=[LFPdata LFPts];
           fprintf([lfpfiles(j).name '\n']);
        end
    end
end

%% This is just testing work

sessionnum=10; % this is the session


session=Region(sessionnum);

% First we should find where the theta peak is, maybe at specific
% timepoints in the trial, so we should stitch togethe different epochs and
% run the spectrogram.  I have a feeling that context explore theta will be
% different from that around sampling.




% now to look at that theta around sampling events
[~,~,inds]=event_spikes(session.lfp(:,2),session.bf.sample_start,1,2);
for i=1:length(inds)
    % make all thetas same length
    inds{i}=[inds{i};ones(3000-length(inds{i}),1)];
    lfpim(:,i)=Region(1).lfp(inds{i},1);
end
figure; imagesc(lfpim');

% Now lets look at the theta power at sampling:
%% This is how you would get your phases and amplitues

% butterworth filter and hilbert transform

%%
LFPdata=session.lfp(:,2);
LFPts=session.lfp(:,1);
band=[15 35];
[Phase,Amplitude1] = GetLFPBand(LFPdata,LFPts,band);

band=[8 10];
[Phase,Amplitude2] = GetLFPBand(LFPdata,LFPts,band);


% get our event spikes
[~,~,inds]=event_spikes(session.lfp(:,1),Session.flags.treadmill-2,1,2);


for i=1:length(inds)
    % make all thetas same length
    inds{i}=[inds{i};ones(3000-length(inds{i}),1)];
    
    thetaim1(:,i)=Amplitude1(inds{i},1);
    thetaim2(:,i)=Amplitude2(inds{i},1);

    spectdata(:,i)=session.lfp(inds{i},2);
    spectdata(:,i)=spectdata(:,i)-mean(spectdata(:,i));
end

figure;
imagesc(thetaim1'); set(gca,'YDir','normal')
figure; imagesc(thetaim2'); set(gca,'YDir','normal')

  
    

%%
% plot out theta power with error background;
%figure;
y=mean(thetaim1,2)';
x=[1:3000];
z=SEM(thetaim1,2)';

figure;
% get shape of patch
xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[1 .8 0],'EdgeColor','none')
pl(1) = line(x,y,'linewidth',2,'color','r');

hold on;
y=mean(thetaim2,2)';
x=[1:3000];
z=SEM(thetaim2,2)';

% get shape of patch
xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[0 .8 1],'EdgeColor','none')
pl(1) = line(x,y,'linewidth',2,'color','b');

legend('high theta',' ','low theta', ' ');

%%
% this uses the chronux multitaper spect continuous data

params.Fs=1000;
params.fpass=[1 90];
params.tapers=[2 3];
params.err=[2 .05];
params.pad=2;

[S,t,f,Serr]=mtspecgramc(spectdata,[0.5 0.025],params);

figure;
sgram1=mean(S,3)';
for i=1:length(sgram1(:,1))
    sgram1(i,:)=sgram1(i,:)*i^2;
end
imagesc(t-1,f,sgram1); set(gca,'ydir','normal');
% now self normalize each band

newS=sgram1;
for i=1:length(newS(1,:))
    newS2(:,i)=newS(:,i)-mean(newS(:,i));
end
figure;
imagesc(t-1,f,newS2);  set(gca,'ydir','normal');


%%
% now lets add the LFP coherence to each cell, probably add it in to the
% existing cells so I can just add it in to the existing dataset
region=Region;
%stats.sess=nan;
for i=2 %3:length(region)
    session=region(i);
    if ~isempty(session.lfp)
        
        dt=.001;
        % low theta is way better for this
        thetafreq=[15 35];  
        
        LFPdata=session.lfp(:,1);
        LFPts=session.lfp(:,2);
        % butterworth filter and hilbert transform
        thetapass=thetafreq.*(2*dt);
        [a,b]=butter(3,thetapass);
        theta=filtfilt(a,b,LFPdata);
        thetaparts=hilbert(theta);
        thetaamp=abs(thetaparts);
        thetaphase=angle(thetaparts);
        % now go through each cell and get the rose plot of theta phase pref:
        % we might want to filter for only times where theta power is above some
        % amplitude
        
        %%%%%%%%%%%%%%%
        % this thresholds theta
        %%%%%%%%%%%%%%%%
        % DONT DO THIS IF YOU ARE USING INTERP1
        % taking the top... 50%?
        % thetathresh=prctile(thetaamp,75);
        thetathresh=median(thetaamp); % you can play with this, may want to
        % get a histogram to figure it out
        ditch=thetaamp<thetathresh;
        thetaamp(ditch)=[]; thetaphase(ditch)=[];
        LFPdata(ditch)=[]; LFPts(ditch)=[];
        
        
        %phases=interp1(LFPts,LFP1',Session.units(1).ts,'nearest');
        %
        for j=1:length(session.units)
            % if you want to use only spikes around a specific event (cant
            % add phases to struct though
            myspikes=event_spikes(session.units(j).ts,session.bf.sample_start,2,2);
            %myspikes=session.units(j).ts(1,:);
            
            % this interpolates but maybe i just want to find the nearest
            % and threshold at a certain latenct
           % phases=interp1(LFPts,[thetaphase LFPts],myspikes,'nearest');
            phases=[];
           for k=1:length(myspikes)
                [idx]=find(LFPts>myspikes(k),1,'first');
                if LFPts(idx)-myspikes(k)<.05
                    phases(k)=thetaphase(idx);
                else
                    phases(k)=nan;
                end
           end
           phases(isnan(phases))=[];
           
           
           
            %MEC(i).units(j).phases=phases;
            figure; %subplot(1,2,1);
            %[f,x]=hist(phases,36);
            %bar(x,zscore(f)); hold on; plot(x,(cos(x))+1,'r');
            %ylim([-2.5 2.5]);
            subplot(1,2,2); rose(phases,25);
            % now the structure is the same order as MEC
            clear temp
%             temp=circ_stats(x,f');
%             stats(i).sess(1,j)=temp.mean;
%             stats(i).sess(2,j)=temp.var;
            % mean vector angle
            stats(i).sess(1,j)=circ_mean(x,f');
            % mean vector length
            stats(i).sess(2,j)=circ_r(x,f');
            %stats(i).name=session.name;
            %session.units(i).phasepref=x(max(f));
            fprintf(num2str(length(session.units)-j));
        end
    else
        for j=1:length(session.units)
            stats(i).sess(1:2,j)=[nan;nan];
        end
    end
    fprintf('\n session %d done \n',i);
    
end


% to calculate the phase prefernece, circ_r is the function, under the
% toolbox that john gave me.
%% First lets get a look at our lfp and theta

%these two cell matrices will be the timestamps and the datapoints for
%each trial
session=MEC(1);
LFPdata=session.lfp(:,1);
LFPts=session.lfp(:,2);

thisflag=session.bf.sample_start;
for i=1:length(thisflag)
    % we'll use a 4 second epoch ending at the clipboard
    expinds=LFPts<thisflag(i)+2 & LFPts>thisflag(i)-2;
    expdata{i}=LFPdata(expinds);
    expts{i}=LFPts(expinds);
end

% now we can concatenate the data
test=cat(1,expdata{:});

% Laras parameters;
% Fs=1000, fpass=1 90, tapers=2,3 err=2 0.05, and pad=1
%
%
%%%%%%%%%%%%%%%%%%%%%

for i=1:length(expdata)
    tslength(i)=length(expdata{i});
    expdata{i}=expdata{i}-mean(expdata{i});
end
tsfloor=min(tslength);
for i=1:length(expdata)
    tempdata{i}=expdata{i}(1:tsfloor);
end
trialdata=cat(2,tempdata{:});
params.Fs=1000; params.fpass=[1 90]; params.tapers=[2 3];
params.err=[2 0.05]; params.pad=1;
[S,t,f,Serr]=mtspecgramc(trialdata,[0.5 0.025], params);
Smean=mean(S,3)';

for i=1:length(f)
    sxx(i,:)=Smean(i,:)*10*log(f(i));
end

% and normalize by frequency
figure;
imagesc(t-2,f(1:30),sxx(1:30,:))
set(gca,'ydir','normal')



%%%%%%%%%%%%%%%%%
[pxx1,f1]=pmtm(test,3,[1:.5:100],1000);
figure; semilogy(f1,smooth(pxx1,2)); title(' sample periodogram');

thisflag=session.bf.sample_start;
for i=1:length(thisflag)
    % we'll use a 3 second epoch surrounding each sample
    trialinds=LFPts<thisflag(i)+2 & LFPts>thisflag(i)-1;
    trialdata{i}=LFPdata(trialinds);
    trialts{i}=LFPts(trialinds);
end

test=cat(1,trialdata{:});

[pxx2,f2]=pmtm(test,3,[1:.5:100],1000);
figure; semilogy(f2,smooth(pxx2,2)); title('sample periodogram');

% about the same
%%
% the wavelet transform of that
for i=1:length(trialdata(1,:))
[wave(:,:,i),period,scale]=DS_wavelet(trialdata(:,i),.001,0,.01,.01,450);
fprintf(['trial ' num2str(i) ' done \n']);
end
%%
% other way of calculating spike field coherence
% from mark kramers class, and it will be across all frequencies
% lets first build our spike trians and then fit them to a trial structure
spikes=Session.unitdata.units(i).ts;
K=length(session.units);
% make a huge vector of zeros
for i=1:K
    spikets=zeros(1,length(LFPts));
    spikesnap=interp1(LFPts,LFPts,spikes,'nearest');
    [~,spikepos]=intersect(LFPts,spikesnap);
    spikets(spikepos)=1;
    % this is a 1000 Hz signal if poisson like spikes
    spiketrain=spikets';
    
    % now we break our data into trials
    % first epoch will be explore
    
    thisflag=Session.flags.treadmill(1:82);
    
    % have to freq axis, lame
    freqs=[(0:5000/2-1)*(1/5), (-5000/2:-1)*(1/5)];
    S11=zeros(length(freqs),1); S12=zeros(length(freqs),1);
    S22=zeros(length(freqs),1);
    
    for j=1:length(thisflag)
        % we'll use a 5 second epoch ending at the clipboard
        T=5;
        expinds=LFPts<thisflag(j) & LFPts>thisflag(j)-5;
        expLFP{j}=LFPdata(expinds);
        expSPK{j}=spiketrain(expinds);
        expSPK{j}=expSPK{j}-mean(expSPK{j});
        
        % compute fft
        Xspk{j}=fft(expSPK{j});
        Xlfp{j}=fft(expLFP{j});
    
        % and get the spectra
        S11= S11 + dt^2/T*(Xspk{j}.*conj(Xspk{j}));
        S12= S12 + dt^2/T*(Xlfp{j}.*conj(Xspk{j}));
        S22= S22 + dt^2/T*(Xlfp{j}.*conj(Xlfp{j}));
    end
    cohr=S12.*conj(S12) ./S11 ./S22;
    figure;
    subplot(3,1,1);
    plot(fftshift(freqs),fftshift(real(S11)));
    xlabel('time'); ylabel(['neuron' num2str(i)]);
     xlim([0 100]);

    subplot(3,1,2);
     plot(fftshift(freqs),fftshift(real(S22)));
    xlabel('time'); ylabel('LFP');
    xlim([0 100]);
    
    subplot(3,1,3);
    plot(fftshift(freqs), fftshift(cohr));
    ylim([0 1]); xlim([0 50]);
    xlabel('freq hz');

end 

% now the question is how to divvy up the cells?  Make a coherence measure
% and a phase preference?

% The chronux approahc to coherence:

%%

spikets=zeros(1,length(LFPts));
spikesnap=interp1(LFPts,LFPts,session.units(i).ts,'nearest');
[~,spikepos]=intersect(LFPts,spikesnap);
spikets(spikepos)=1;
% this is a 1000 Hz signal if poisson like spikes
spiketrain=spikets';

% now we break our data into trials
% first epoch will be explore

thisflag=session.bf.clipboard;
for j=1:length(thisflag)
    % we'll use a 5 second epoch ending at the clipboard
    T=5;
    expinds=LFPts<thisflag(j) & LFPts>thisflag(j)-5;
    expLFP{j}=LFPdata(expinds);
    expSPK{j}=spiketrain(expinds);
    expSPK{j}=expSPK{j}-mean(expSPK{j});
end

% continuous, point binned
[c,~,Sxx,Syy,Sxy,f]=coherencycpb();

%%
for i=1:length(MEC)
    thisname=MEC(i).name;
    newname=thisname(1:find(thisname=='-',1,'last')+2);
    MEC(i).name=newname;
end

