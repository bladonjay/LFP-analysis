% peek Whole Session LFP Band


function [LFPsession] = PeekSessionLFP(bf, lfpdata, ts, dt, freqband)
%

% the variables needed are:

% trialdur
trialdur=mean(diff(bf.trial));

%   flag
flag=bf.clipboard_up;



%   freqband


%   bins
bins=80;  % fifty bins?

%   dt

%   binstep

binstep=(trialdur/bins);

%   timebins
timebins=[-trialdur/2:binstep:0  binstep:binstep:trialdur/2];
% theta band pass with butterworth

freqpass=freqband.*(2*dt);
% third order filter (three exponentials)
[a,b]=butter(3,freqpass);
thetalfp=filtfilt(a,b,lfpdata);


theta=hilbert(thetalfp);
phi=angle(theta);
amp=abs(theta);





% right now this is about .8 of a second

% if we center our trial over first sample
% if we center our trial over trial start
%timebins=0:binstep:trialdur;




%
% basically get the mean for each bin for each trial

thetasess=zeros(length(flag),length(timebins));
for q=1:length(flag)
    
    for r=1:length(timebins)
        timestart=flag(q)+timebins(r);
        timeend=flag(q)+timebins(r)+binstep;
        indices=ts>timestart & ts<timeend;
        %build the image
        thetasess(q,r)=nanmean(amp(indices));
    end
    
    fprintf('trial %d finished \n',q);
end


figure;
plot(ts,smooth(amp,30));
hold on
plot(flag,ones(1,length(flag))*.3,'r*');

for r=1:length(thetasess(:,1))   
    LFPsession(r,:)=zscore(smooth(thetasess(r,:),2));
end

figure
imagesc(LFPsession)
fixtix=round(linspace(0,trialdur(end),11)-trialdur(end)/2);
set(gca,'Xtick',[1:(max(bins)/10):max(bins)],'XTickLabel',[fixtix])
% set name and the lfp band

%title(myname);
binsize=num2str(trialdur(end)/bins);
xlabel(['time to divider up binsize ' binsize 'sec']); ylabel('trial number');

end


