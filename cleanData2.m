function [pctDataUsed,inds2cut,varargout] = cleanData2(root,filtType,filtParams,varargin) %gammaamps,thetaangles)

% root fields:
%   root.lfp.signal
%   root.lfp.fs
%   root.epoch


  varargout = []; if ~exist('varargin','var'), varargin = []; end;

  % compute specific filters if needed
  if exist('filtType','var') && ~isempty(filtType)
   
    % NOTE: if you want scaled EEG, use scaledEEG.m to update the root object before submitting it to cleanData2;

    
    % find good amps
    switch(filtType)
      case 'theta' % relative amplitude of theta in percentage of it max value
          % command in would be 'theta',.5 (50% of max theta)
        th_Amp = abs(hilbert(buttfilt(root.lfp.signal,[6 12],root.lfp.fs,'bandpass',4)));
        excInds = th_Amp < max(th_Amp)*filtParams(1) | th_Amp > max(th_Amp)*filtParams(2);
      case 'thetaMag' % absolute magnitude of theta in mV
        th_Amp = abs(hilbert(buttfilt(root.lfp.signal,[6 12],root.lfp.fs,'bandpass',4)));
        excInds = th_Amp < filtParams(1) | th_Amp > filtParams(2);
      case 'thetaDeltaRatio'
        th_Amp = abs(hilbert(buttfilt(root.lfp.signal,[6 10],root.lfp.fs,'bandpass',4))); 
        dl_Amp = abs(hilbert(buttfilt(root.lfp.signal,[1 3],root.lfp.fs,'bandpass',4)));
        rat = th_Amp ./ dl_Amp;
        excInds = true(size(rat));        
        incEpchs = CMBHOME.Utils.OverThresholdDetect(rat,filtParams(1),ceil(root.lfp.fs/4),ceil(root.lfp.fs/2));
        if size(incEpchs,1)
          root.b_lfp(root.active_lfp).b_myvar = [1:length(root.b_lfp(root.active_lfp).signal)]';
          OE = root.epoch; root.epoch = root.lfp.ts(incEpchs); incInds = CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar); 
          root.epoch = OE; incInds = incInds - root.lfp.myvar(1) + 1;
          excInds(incInds) = false;
        end
        %         [~,th_Amp] = multiphasevec(6:12,root.lfp.signal,root.lfp.fs,6);  th_Amp = mean(th_Amp);% Find theta power
        %         [~,dl_Amp] = multiphasevec(1: 3,root.lfp.signal,root.lfp.fs,6);  dl_Amp = mean(dl_Amp);% Find delta power
        %        excInds = th_Amp ./ dl_Amp <= filtParams(1);
      case 'vel'
        hdvel = interp1(root.ts,root.vel*root.spatial_scale,root.lfp.ts);
        excInds = hdvel < filtParams(1) | hdvel > filtParams(2);
      otherwise
        error('Unknown filtType');
    end
    excInds = reshape(excInds,1,[]);
  else
    excInds = zeros(1,length(root.lfp.signal));
  end

  % record original data length
  ODL = length(root.lfp.signal);
  
  % drop saturated or flatlined data
  t0 = [0-1e-4 0+1e-4]; % tresholds for no dV/dt
  badInds = CMBHOME.Utils.ThresholdBandDetect(diff(root.lfp.signal),t0(1),t0(2),1,round(0.030 * root.lfp.fs)); % flatline data
  if ~isempty(badInds)
    OE = root.epoch;
    % add 500ms buffer around each bit of bad data to remove potential edge effects on the pahse estimates
    badInds = root.lfp.ts(min(max([badInds(:,1)-ceil(root.lfp.fs/2) badInds(:,2)+ceil(root.lfp.fs/2)],1),length(root.lfp.ts)));
    if isvector(badInds), badInds = reshape(badInds,1,[]); end
    % define epochs as only the bad data with the surrounding buffers
    root.epoch = badInds;
    root = CMBHOME.Utils.MergeEpochs(root);
    % build an 'inds' vector in lfp object to find inds of the bad data
    root.b_lfp(root.active_lfp).b_myvar = [1:length(root.b_lfp(root.active_lfp).signal)]';
    % build bad ind boolean
    badInds = zeros(length(root.b_lfp(root.active_lfp).signal),1);
    % use myvar to mark bad data epochs to set the boolean flags
    badInds(CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar)) = 1;
    % finally, repass this boolean matrix through to get appropriately windowed vector
    root.b_lfp(root.active_lfp).b_myvar = badInds;
    % restore epoch to normal so we can still refer to root.lfp
    root.epoch = OE;
    % grab correctly epoched data out
    badInds = CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar)';
  else
    badInds = zeros(1,length(root.lfp.signal));
  end
  % perform filter
  inds2cut = excInds|badInds;
  for i = 1:length(varargin)
    if isvector(varargin{i}) 
      varargin{i}(excInds|badInds) = []; 
    else
      varargin{i}(:,excInds|badInds) = []; 
    end
  end
  varargout = varargin;

  % document how much data remains
  pctDataUsed = sum(~inds2cut)/ODL;
end


