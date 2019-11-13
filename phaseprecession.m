function [spikephi,spikeloc] = phaseprecession(spikes,lfp,tracking,varargin)
%PHASEPRECESSION
% builds an xy scatter of x position by y phase
%spikes - Nx1 vector or cell array of vectors.  Spike timestamps
%lfp - Mx2 matrix.  Column 1: lfp timestamps.  Column 2: lfp signal
%tracking - Tx2 matrix.  Column 1: tracking timestamps.  Column 2:
%       linearized positoin
%
%Jon Rueckemann 2014

%Default values
fpass=[8 12];
plotscatter=true;
locmethod='pchip'; %choose 'nearest' for no spatial interpolation

%Extract input arguments
extractvarargin(varargin);

%Obtain phase angle for each lfp timestamp
lfpts=lfp(:,1);
lfp=lfp(:,2);
fpass=fpass*(2/Fs);
[a,b] = butter(3,fpass); 
lfp=filtfilt(a,b,lfp);
hlfp=hilbert(lfp);
phi=angle(hlfp);

%Interpolate instantaneous phase and location for each spike
[spikephi,spikeloc]=deal(cell(numel(spikes),1));
for s=1:numel(spikes)
    spikephi{s}=interp1(lfpts,phi,spikes,'pchip');
    spikeloc{s}=interp1(tracking(:,1),tracking(:,2),spikes,locmethod);
    if plotscatter
        figure;
        scatter(spikeloc,spikephi);
    end
end