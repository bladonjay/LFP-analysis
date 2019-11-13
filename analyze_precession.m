function [rho,p,s,b,spk_flds_pos] = analyze_precession(spk_ts,lfp_sig,lfp_phase,varargin)
%ANALYZE_PRECESSION - Analyzes phase precession
%
% Analyzes phase precession as in Hafting et al (2008)
%
% INPUT
%   spk_ts - Time of each spike
%   lfp_sig - LFP sifnal
%
% OPTIONAL PARAMETERS
%   lfp_phase ([]) - Phase of lfp_sig.  If [], calculates from lfp_sig
%
% PARAMETERS
%   spk_flds_pos ([]) - position of spikes in the fields.  If empty,
%       include other parameters to run find_fields. See doc find_fields
%       for more information.
%   lfp_fs (250) - Sampling frequency for the LFP
%   startt (0) - Start time
%   theta_range ([6 10]) - Frequency band for theta
%   spk_theta_phase ([]) - Theta phase of each spike.
%   plottit (true) - Plots the precession.
%   vid_ts ([]) - The timestamps of state.  If [], uses ((1:numel(state))-1)/fs+startt
%   
% RETURNS
%   rho - Linear circular correlation (Kempter, 2012)
%   p - Significance of rho (Kempter, 2012)
%   s - Slope of linear-circular regression (Kempter, 2012)
%   b - y-intercept of linear-circular regression
%   spk_flds_pos - Positions of the spikes in each field
%
% This code has been freely distributed by the authors as part of linear_track_precession.
%
% Release notes
%   v0.1 - 2014-03-03 Jason Climer (jason.r.climer@gmail.com)
%   v0.11 - 2014-05-05 Added vid_ts Jason Climer (jason.r.climer@gmail.com)
if ischar(lfp_phase)
    varargin = [{lfp_phase} varargin];
    lfp_phase = [];
end

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParamValue('spk_flds_pos',[]);
ip.addParamValue('lfp_fs',250);
ip.addParamValue('startt',0);
ip.addParamValue('theta_range',[6 10]);
ip.addParamValue('spk_theta_phase',[]);
ip.addParamValue('vid_ts',[]);
ip.addParamValue('plottit',1);
ip.parse(varargin{:});
for i=fields(ip.Results)', eval([i{1} '=ip.Results.' i{1} ';']); end;

%Subplots
if plottit, h = gca; pos = get(h,'Position'); axis off; end;

%Find fields
if isempty(spk_flds_pos)
   try
       if plottit, h = axes('Position',[pos(1) pos(2)+pos(4)/2 pos(3) pos(4)/2]); end
       [~,spk_flds_pos] = find_fields('spk_ts',spk_ts,varargin{:});
   catch err
       throw(MException('analyze_precession:nofields','If spk_flds_pos not provided, input must work for find_fields'));
   end
end

if isempty(vid_ts)
    vid_ts = ((1:numel(signal_filtered+1))-1)/lfp_fs+startt;
end

%Spike theta phase
if isempty(spk_theta_phase)
    if isempty(lfp_phase)
        [btheta,atheta]=butter(3,theta_range/(lfp_fs/2));
        signal_filtered = filtfilt(btheta,atheta,lfp_sig);
        lfp_phase = atan2(imag(hilbert(signal_filtered)), signal_filtered);
    end    
    [~,spk_theta_phase]=histc(spk_ts,vid_ts);
    spk_theta_phase = lfp_phase(spk_theta_phase);
end

% Fits and correlations
b = NaN(size(spk_flds_pos,2),1);
s = NaN(size(spk_flds_pos,2),1);
p = NaN(size(spk_flds_pos,2),1);
rho = NaN(size(spk_flds_pos,2),1);

for i=1:size(spk_flds_pos,2)
    goodspks = ~isnan(spk_flds_pos(:,i))&spk_flds_pos(:,i)>=0&spk_flds_pos(:,i)<1;
    if ~isempty(goodspks)
   [rho(i),p(i),s(i),b(i)] = kempter_lincirc(spk_flds_pos(goodspks,i),spk_theta_phase(goodspks));
    end
end

% Plotting
if plottit
    
    if isempty(ip.Results.spk_flds_pos)
    hi = pos(4)/2;
    else
    hi = pos(4);
    end

    for i=1:size(spk_flds_pos,2)
       axes('Position',[pos(1)+pos(3)/size(spk_flds_pos,2)*(i-1) pos(2) pos(3)/size(spk_flds_pos,2) hi]); 
       goodspks = ~isnan(spk_flds_pos(:,i))&spk_flds_pos(:,i)>=0&spk_flds_pos(:,i)<1;
       scatter(repmat(spk_flds_pos(goodspks,i),[2 1]),rad2deg([mod(spk_theta_phase(goodspks),2*pi);mod(spk_theta_phase(goodspks),2*pi)+2*pi]),'.');
       
       x = linspace(-1,1,750);
       phi = mod(2*pi*s(i)*x+b(i),2*pi);
       k = find(abs(diff(phi))>pi);
       phi(k) = NaN;
       
       hold on;
       text(x(floor(750*0.75)),rad2deg(phi(floor(750*0.75))+2*pi)+25,['rho=' sprintf('%2.2g',rho(i))],'Color',[1 0 0],'FontWeight','bold','BackgroundColor','w');
       text(x(floor(750*0.75)),rad2deg(phi(floor(750*0.75))+2*pi)-25,['p=' sprintf('%2.2g',p(i))],'Color',[1 0 0],'FontWeight','bold','BackgroundColor','w');
       plot(x,rad2deg(phi+2*pi),'r','LineWidth',2);
       plot(x,rad2deg(phi),'r','LineWidth',2);
       
       xlim([0 1]);ylim([0 720]);xlabel('% Through Field');ylabel('\theta phase');
       axis square;
    end
end

axes(h);

end

