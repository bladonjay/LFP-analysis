function [modindex, thetarange, gammarange, powPhsDists, bincenters,pctDataUsed,thetaamps_M,gammaamps_M,stdVals] = thetaModGamma2(root,varargin)
% function [modindex, thetarange, gammarange, powPhsDists, bincenters] = thetaModGamma(root,[optional inputs])
%
% Computes the phase modulation of gamma by theta using the method
% described by Tort et al (2010) J Neurophys.  In short, it computes the
% entropy in the histogram of powers over phases.
%
%
% INPUT ARGS
%  root - (required) cmb object with active_lfp set and lfp loaded
%  optional inputs:
%    filtType - default = [], otherwise 'theta', 'thetaDeltaRatio' or 'vel' to process data only
%               with high theta power or fast running speeds.
%    filtParams - default = [], lower and upper bounds to use in the
%                 filtering of the data as mentioned above
%    shuffles - default 0, use if you want to compute null dist of modulation
%    thetaRange - default = [4:0.25:12] frequencies to look at phases of 
%    gammaRange - defualt = [30:1:120] frequencies to look at the power of
%    stdGamma - default = 0 for no standardization. Give 1 to standardize 
%                 current data, or provide and N by 1 vector with the scale 
%                 factor for each of the N gamma bands (usually the average 
%                 over bins).
%
%
% OUTPUT ARGS
%  modindex - matrix of size [nG x nT] where nT is the number of theta frequencies 
%             analyzed and nG is the number of gamma frequencies analysed
%  thetarange - same as inputed, usable as labels for the matrix when plotting
%  gammarange - same as inputed, usable as labels for the matrix when plotting
%  powPhsDists - normalized power histograms over phases for theta band that 
%                showed the max modulation. has shape [nG x nT x nPh] where
%                nT is the number of theta frequencies analyzed and nG is
%                the number of gamma frequencies analysed and nPh is the
%                number of phase bins used.
%  bincenters - phase bin labels
%  stdValues  - pre-standardization sum for each of the N gamma bands used in
%                standardization.  0 if no zscoring was used.
%
% v2 - modified to use single broadband theta for phase estimates 6-10Hz

% eln 120512
% eln 130308 - updated to parse inputs and use reference theta channel
% eln 130315 - updated to z-score gamma bands

  p = inputParser;

  p.addRequired('root')
  p.addParamValue('filtType', [], @ischar);
  p.addParamValue('filtParams', [], @isnumeric);
  p.addParamValue('shuffles',0,@isnumeric);
  p.addParamValue('thetarange',4:0.25:12,@isnumeric);
  p.addParamValue('gammarange',30:1:120,@isnumeric);
  p.addParamValue('broadband',0,@isnumeric);
  p.addParamValue('thetaRefChan',root.active_lfp,@isnumeric);
  p.addParamValue('stdGamma',0,@isnumeric);
  p.addParamValue('phaseType','wavelet',@ischar);
  p.addParamValue('useHDeeg',1,@isnumeric);

  p.parse(root, varargin{:});

  root = p.Results.root;
  filtType = p.Results.filtType;
  filtParams = p.Results.filtParams;
  shuffles = p.Results.shuffles;
  thetarange = p.Results.thetarange;
  gammarange = p.Results.gammarange;
  broadband = p.Results.broadband;
  thetaRefChan = p.Results.thetaRefChan;
  stdGamma = p.Results.stdGamma;
  phaseType = p.Results.phaseType;
  useHDeeg = p.Results.useHDeeg;
  
  % check that everything gets defined
  if isempty(filtParams) && ~isempty(filtType), 
    error('filtParams must be defined if filtType is specified');
  elseif isempty(filtType),
    filtParams = [];
	end
  
  % select out high quality data (no saturations and no flat lines + apply requested filtering
  [pctDataUsed,inds2cut] = cleanData2(root,filtType,filtParams);
  
  % quit now if pctDataUsed is too low (there should be at least 200 cycles of theta (7Hz), the min mentioned by Tort et al (2010).)
  if pctDataUsed*sum(diff(root.epoch))<(200/7) | isnan(gammarange), 
    [modindex, thetarange, gammarange, powPhsDists, bincenters,thetaamps_M,gammaamps_M,stdVals] = deal(nan);
    return;
  end

  % get gamma amplitudes
  % this grabs the power at each gama band (bands are above)
  [~,gammaamps] = multiphasevec(gammarange,root.lfp.signal,root.lfp.fs,6); 
  if broadband, gammaamps = mean(gammaamps); end; 
  clear signal
  
  % get theta filtered signals and phase (basically a butter fitler aroiund
  % theta, or look for which ones)
  root_ = root;
  if thetaRefChan~=root.active_lfp
    root_ = root_.LoadLFP(thetaRefChan);
    root_.active_lfp = thetaRefChan;
    [~,~,root_] = scaledEEG(root_,useHDeeg,1,250);
  end
  
  % get the phase
  [thetaangles,thetaamps,thetarange] = extractThetaPhase(root_.lfp.signal,root_.lfp.fs,phaseType,thetarange);
  clear root_ 
  
  % apply data clean up filter now
  gammaamps(:,inds2cut) = [];
  thetaangles(:,inds2cut) = [];

  % return mean theta and gamma amps
  thetaamps_M = mean(thetaamps,2);
  gammaamps_M = mean(gammaamps,2);
  
  % standardize gammaamps as needed
  if any(stdGamma~=0) && ~any(isnan(stdGamma)), [gammaamps,stdVals] = standardizeGamma(gammaamps,stdGamma); end
  
  % double check that there are at least 200 cycles of theta (7Hz), the min
  % mentioned by Tort et al (2010).
  if size(thetaangles,2)<200*(root.lfp.fs/7) || any(isnan(stdGamma))
    modindex = nan;
    powPhsDists = nan;
    bincenters = nan;
  else
    
    if ~shuffles
      % compute the values of interest
      [modindex, powPhsDists, bincenters] = computeModIndex(gammaamps,thetaangles,thetarange,root.lfp.fs);  
    else
      % compute significance thresholds for each frequency pairing
      [modindex, powPhsDists, bincenters] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange,root.lfp.fs);
    end

    % filter out phs values for theta band with peak amp
    powPhsDists = permute(powPhsDists,[1 3 2 4]);
    %[~,highAmpInd] = max(mean(thetaamps,2));
    %powPhsDists = squeeze(powPhsDists(:,highAmpInd,:));
  end
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the JUICY stuff
function [modindex, powPhsDists, bincenters] = computeModIndex(gammaamps,thetaangles,thetarange,fs)
  % bin gamma amplitudes by theta phases
  meanamps = nan(36,size(gammaamps,1),size(thetaangles,1));
  bincenters = nan(36,size(thetaangles,1));
  for i=1:size(thetaangles,1)
    [meanamps(:,:,i),bincenters]=Findmeanamps(thetaangles(i,:)',gammaamps',thetarange(i),fs);
  end
  
  % get back to the way things used to be
  meanamps = permute(meanamps,[2 3 1]);

  % normalize the area to one
  meanampsNorm = meanamps ./ repmat(sum(meanamps,3),[1,1,size(meanamps,3)]);  

  % compute the mod indexes
  shannonent = sum(meanampsNorm .* log10(meanampsNorm),3);
  modindex = (log10(36)+shannonent) ./ log10(36);

  % power distrobution over theta phases for the theta band with the
  powPhsDists = meanamps;
end


function [modindex_out, powPhsDists_out, bincenters] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange,fs)
  fprintf('Shuffling (%i): ',shuffles);
  modindex = nan(size(gammaamps,1),size(thetaangles,1),shuffles);
  powPhsDists = nan(size(gammaamps,1),size(thetaangles,1),36,shuffles);

  for i = 1:shuffles
    % keep user feel confident that all is well
    fprintf('%i ',i);
    
    % find random offset to shift all theta angles by
    offsetOkay = false;
    while ~offsetOkay
      offset = round(rand*size(thetaangles,2));
      offsetOkay = offset>1 && offset<size(thetaangles,2);
    end

    % shift theta angles
    ta = [thetaangles(:,offset+1:end), thetaangles(:,1:offset)];
    % compute mod indices for this shift
    [modindex(:,:,i), powPhsDists(:,:,:,i), bincenters] = computeModIndex(gammaamps,ta,thetarange,fs);
  end

  % find distro values
  modindex_out = modindex;
  powPhsDists_out = powPhsDists;

  %   modindex_out(:,:,1) = mean(modindex,3);
%   modindex_out(:,:,2) = std(modindex,[],3);
%   powPhsDists_out(:,:,1) = mean(powPhsDists,3);
%   powPhsDists_out(:,:,2) = std(powPhsDists,[],3);
   fprintf('\n');
end


function [meanamps bincenters]=Findmeanamps(slowphase,fastamp,thetaFrq,fs)
  %this function finds the normalized mean amplitude of the fast wave for each phase of the
  %slow wave
  %it recieves the phase of the slow wave (in radians from 0-2pi) and the
  %corresponding amplitude of the fast wave (in mV)

  %it returns a matrix with the normalized mean amplitude (row 2) corresponding to the
  %phases(row 1-- max of corresponding phase), and the middle of the phase
  %for a bar graph(row 3)

  
  
  %%make into column vectors if not
  if size(slowphase,1)==1
    slowphase=slowphase';
  end

  % double check data matches
  if size(fastamp,1)~=size(slowphase,1)
    error('fastamp must have time in rows');
  end

  % define phase bins
  baseSorting = linspace(-pi, pi, 37);
  sorting = linspace(-pi, pi, floor((fs/thetaFrq)/2)+1);
  
  % perform sorting
  meanamps = nan(length(sorting)-1,size(fastamp,2));
	binN = nan(length(sorting),1);
  for bin = 1:length(sorting)-1
    currInds = slowphase > sorting(bin) & slowphase <= sorting(bin+1);
    meanamps(bin,:) = nanmean(fastamp(currInds,:),1);
		binN(bin) = length(currInds);
  end  
  
  % add in bin labsls
  bc = mean([sorting(1:end-1);sorting(2:end)])'+pi;
  bincenters = mean([baseSorting(1:end-1);baseSorting(2:end)])'+pi;

  % remap phases back onto 36 bins
  meanamps = interp1([bc(1)-mean(diff(bc));bc;bc(end)+mean(diff(bc))],[meanamps(end,:);meanamps;meanamps(1,:)],bincenters);
  
%    bc = mean([sorting(1:end-1);sorting(2:end)])'+pi; bincenters = mean([baseSorting(1:end-1);baseSorting(2:end)])'+pi;meanamps2 = interp1([bc(1)-mean(diff(bc));bc;bc(end)+mean(diff(bc))],[meanamps(end,:);meanamps;meanamps(1,:)],bincenters);

end