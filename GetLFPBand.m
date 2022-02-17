function [Phase,Amplitude,Filtered] = GetLFPBand(LFPdata,LFPts,band,method)
% function [Phase,Amplitude,Filtered] = GetLFPBand(LFPdata,LFPts,band,correction)
% GetLFPBands pulls the phases and amplitudes of your data in the freq band
% you're looking for.  


% methods 0 and 1 Use a butterworth filter then takes the hilbert
% transorm.  Does not threshold the dataset. This is a rough approximation
% of the dataset.  Data should come in at 1000 Hz.



% method 2 takes the wavelet, finds the maximum amplitude at each time
% point and uses the phase from that frequency.  You will get the filtered
% signal as well



if ~exist('band','var')
    freq=[5 12];
    fprintf('using 5-12 Hz \n');
else
    freq=band;
end

% if no corr var or if its neither 0 or 1
if ~exist('method','var')
    method=1;
elseif method~=0 && method ~=1
    method=1;
end

tempdt=min(diff(LFPts));

dt=.001;
fNQ=1/2/dt;
% get params for bandpass filter

if method <2
thetapass=freq.*(2*dt);

% use a 3 pole butterworth filter
[a,b]=butter(3,thetapass);

Filtered=filtfilt(a,b,LFPdata);

% for other filtereing methods
%{
% or a digital bandpass iir of order 20
bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
         'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',band(2), ...
         'SampleRate',1000);
Filtered=filtfilt(bpFilt,LFPdata);
     
     
% or a first order fiir this takes wayyyyyy long because its order 8000
filtorder=3*fix(band(1)/dt);
fi=[0 (1-.15)*band(1)/fNQ band(1)/fNQ band(2)/fNQ (1+.15)*band(2)/fNQ 1];
mi=[0  0                  1           1            0                  0] ;
filtwts = firls(filtorder,fi,mi);             % get FIR filter coefficients
Filtered = filtfilt(filtwts,1,LFPdata);
%}
% get phase and amplitude using the hilbert trnaform
thetaparts=hilbert(Filtered);
% these are the parts
Amplitude=abs(thetaparts);
Phase=angle(thetaparts);
elseif method == 2
    % wavelet the data in 10 linear increments within your band
    wname='cmor16-.5'; % a 16 x envelope of a morlet at 0.5 hz central frequency
    freqs=linspace(band(1),band(2),10); % log spaced vector to ~250 hz
    bases=500./freqs;
    wtband=cwt(LFPdata,bases,wname);
    allAmps=abs(wtband);
    AllPhases=angle(wtband);
    [Amplitude,maxinds]=max(allAmps,[],1);
    fullinds=sub2ind(size(allPhases),maxinds,1:size(allPhases,2));
    Phase=allPhases(fullinds);
end
% and here we correct for the asymmetry of the oscillation phases, as it
% can be saw-toothed. This basically makes sure the distribution of phases
% is uniform, so each phase occurs equally often
if method==1
    %[cdfphase,newphase]=sort(Phase);
    %[mycdf,x]=ecdf(cdfphase);
    %[~,ib]=sort(newphase);
    %Phase2=mycdf(ib)*2*pi-pi;
    %Phase=Phase2;
    % first get the values
    Phase=tiedrank(Phase);
    Phase=((Phase-1)/max(Phase))*2*pi-pi;
end



end

