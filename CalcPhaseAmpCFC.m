function [comod,phasedistout,phasecenters,p,nulldata,nullsummary] = CalcPhaseAmpCFC(phases,amps,varargin)
% calculates the comodulation index of phase freq and amp freq
% Inputs:
%   phases: a vector of phases
%   amp:    a vector of amplitudes
%   bins:   A number of bins to use from -pi to pi (default is 50)


p=inputParser;
addOptional(p,'LoAmps',[]);
addOptional(p,'Cutoff',0);
addOptional(p,'bins',50,@(a) isnumeric(a) && isscalar(a));
addOptional(p,'runstat',0,@(a) isnumeric(a) && isscalar(a));
parse(p,varargin{:});
bins=p.Results.bins;
runstat=p.Results.runstat;
Cutoff=p.Results.Cutoff; % as a zscore?
LoAmps=p.Results.LoAmps;

% get a real value for the cutoff if we dont designate it
if Cutoff==0 || isempty(Cutoff)
    LoMu=nanmean(LoAmps(:)); LoSig=nanstd(LoAmps(:));
    Cutoff=LoMu+LoSig;
end

% convert LoAmps to logial array
if isempty(LoAmps)
    LoAmps=true(size(phases));
else
    LoAmps=LoAmps>Cutoff;
end

% get these into linearized vectors of the right dimension
if ndims(phases)==1
    phasescat=phases(:); 
else
    phasescat=linearize(phases(LoAmps)'); % linearize, but transpose so it tacks consistently
end
    
if ndims(amps)==1
    ampscat=amps(:);
else
    ampscat=linearize(amps(LoAmps)');
end

% generate our phases
phasebins = linspace(-pi,pi,bins+1);
phasecenters = mean([phasebins(1:end-1); phasebins(2:end)]);
[phasects,~,phaseall] = histcounts(phasescat,phasebins);
phaseAmp = zeros(bins,1);
%this will be a 50 bins of phases where in each is the mean
%amplitude of the amplitude channel across time
% for bb = 1:bins
%    phaseAmp(bb) = mean(amps(phaseall==bb));
% end
phaseAmp=accumarray(phaseall,ampscat)./phasects'; % the faster way, gup array and then divide
phasedist= phaseAmp./sum(phaseAmp); % normalize to sum of it all
phasedistout=gather(phasedist);
%comod will be sum across inds (amp * log(amp/50)) all over
%log numbins...which i think is the sum over
%i think this should be a magnitude squared coherence
%function, which should be available by mscohere(phase signal,
%amplitude signal)
comod = gather(sum(phasedist.*log(phasedist./(ones(bins,1)/bins)))/log(bins));

% from tort et al PNAS 2008 (and from
% github.com/tortlab/phase-amplitude-coupling
%MI=(log(bins)-(-sum(phasedist).*log(phasedist)))/log(bins);
% bootstrap
guardrail=1000; % msec
verbose=0;
phasedistB=nan(runstat,bins);
if runstat>0
    for j=1:runstat
        
        if ndims(amps)==1
            cutpos=randi(length(ampscat)-guardrail*2)+guardrail;
            ampsB=circshift(ampscat,cutpos);
            % now cut to size for lo freq amplitude
        else
            % shuffle trialwise
            trialmix=randperm(size(amps,1));
            mixamps=amps(trialmix,:); % mix amps by trial
            ampsB=linearize(mixamps(LoAmps(trialmix,:))'); % now keep that mix for the loamp cutoff and linearize
        end       
        phaseAmp=accumarray(phaseall,ampsB)./phasects'; % the faster way
        phasedistB(j,:)= phaseAmp./sum(phaseAmp); % normalize to sum of it all
        comodB(j) = gather(sum(phasedistB(j,:).*log(phasedistB(j,:)./(ones(1,bins)/bins)))/log(bins));
    end
    if verbose
    fprintf('finished boot\n');
    end
    nulldata=phasedistB;
    nullsummary=comodB;
    p=1-normcdf(comod,nanmean(comodB),nanstd(comodB));
else
    p=nan; nulldata=nan; nullsummary=nan;
end




end

%{
Sanity check to see if the distribution really is above chance
errsig=2*nanstd(phasedistB,1,1);
errmu=nanmean(phasedistB,1);

figure; patch([1:bins bins:-1:1], [errmu+errsig fliplr(errmu-errsig)],[.7 .7 .7],...
    'FaceAlpha',.5,'LineStyle','none');
hold on;
plot(1:bins,phasedist);
%}


% this code was adapted from the below bit used by antonio fernandez ruiz

% [comod] = bz_CrossFreqMod(lfp,phaserange,amprange,flagPlot)
%
%
%This function calculates the modulation index of phase-amplitude between
%phaserange (lower frequencies) to amplitude range (higher frequencies).
%It can really take a long time if you do very small steps of frequency,
%due to wavelet processing each frequency at a time.
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    phaserange     [min:steps:max] array of frequencies to filter for the
%                   phase signal
%    amprange       [min:stepsmax] array of frequencies range for wavelets
%                   for the power signal
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     phaseCh      channel to compute phase. If empty take first channel
%     ampChans     channels to compute amplitude. If empty take first channel
%     makePlot      default true
%    =========================================================================
%
%OUTPUT
%   comod               phase-frequency versus amplitude-frequency
%                       comodulogram matrix in modulation index units
%
%Dependencies
%   bz_Filter
%   bz_WaveSpec
% Implemented by Eliezyer de Oliveira, 2018
% Update: 08/06/2018%

%   AntonioFR, 2/2019

