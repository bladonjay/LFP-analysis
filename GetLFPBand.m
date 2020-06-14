function [Phase,Amplitude,Filtered] = GetLFPBand(LFPdata,LFPts,band,correction)
% function [Phase,Amplitude,Filtered] = GetLFPBand(LFPdata,LFPts,band,correction)
% GetLFPBands pulls the phases and amplitudes of your data in the freq band
% you're looking for.  Uses a butterworth filter then takes the hilbert
% transorm.  Does not threshold the dataset. This is a rough approximation
% of the dataset.  Data should come in at 1000 Hz.


% butterworth filter and hilbert transform


if ~exist('band','var')
    freq=[5 12];
    fprintf('using 5-12 Hz \n');
else
    freq=band;
end

% if no corr var or if its neither 0 or 1
if ~exist('correction','var')
    correction=1;
elseif correction~=0 && correction ~=1
    correction=1;
end

tempdt=min(diff(LFPts));

dt=.001;
fNQ=1/2/dt;
% get params for bandpass filter
thetapass=freq.*(2*dt);

% use a 3 pole butterworth filter
[a,b]=butter(3,thetapass);

Filtered=filtfilt(a,b,LFPdata);

%{
% or a digital bandpass iir of order 20
bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
         'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',band(2), ...
         'SampleRate',1000);
Filtered=filtfilt(bpFilt,LFPdata);
     
     
% or a first order fiir
filtorder=3*round(dt*band(1));
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

% and here we correct for the asymmetry of the oscillation phases, as it
% can be saw-toothed.
if correction==1
    [cdfphase,newphase]=sort(Phase);
    [mycdf,x]=ecdf(cdfphase);
    [~,ib]=sort(newphase);
    Phase2=mycdf(ib)*2*pi-pi;
    Phase=Phase2;
end



end

