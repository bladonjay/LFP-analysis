function [thetaPhs, thetaAmp, band] = extractThetaPhase(signal,fs,method,band)

if isempty(signal)
    error('Please load lfp with active_lfp before calling');
end

switch(lower(method))
    
    case 'waveform'
        
        min_sep = fs / 20; % 50ms separation between peak and next trough
        min_length = min_sep; % same for the separation between trough and next peak
        
        %  method described by
        broadBand = [1 80];
        signal = buttfilt(signal,broadBand,fs,'bandpass',4);
        
        % find filtered theta peaks and troughs
        if exist('band','var') & ~isempty(band),	thetarange = minmax(band);
        else thetarange = [6 10]; end
        thetaPhs = angle(hilbert(buttfilt(signal,thetarange,fs,'bandpass',4))) + pi;
        
        % find bounds of peaks as pi/2 < peaks < 3pi/2
        peakBounds = CMBHOME.Utils.ThresholdBandDetect(thetaPhs, pi/2, 3*pi/2, min_sep, min_length);
        
        % find bounds of troughs as  pi/2 < peaks < 3pi/2 after adding pi
        trphBounds = CMBHOME.Utils.ThresholdBandDetect(mod(thetaPhs+pi,2*pi), pi/2, 3*pi/2, min_sep, min_length);
        
        % find mean peak duration
        meanPeakDur = ceil(mean(diff(peakBounds,[],2)));
        peak_iM = ones(size(peakBounds,1),meanPeakDur); % prepare to pull 150 samples of data around each event
        peak_iM(:,1) = peakBounds(:,1);
        peak_iM = cumsum(peak_iM,2); % build indice matrix
        peakBounds(max(peak_iM,[],2)>length(signal),:) = [];% drop cycle which extend past the end of signal
        peak_iM(max(peak_iM,[],2)>length(signal),:) = []; % drop cycle which extend past the end of signal
        peakSig = signal(peak_iM);
        peakInd = peakBounds(:,1) + maxInd(peakSig,[],2) - 1;
        
        % find mean trough duration
        meanTrphDur = ceil(mean(diff(trphBounds,[],2)));
        trph_iM = ones(size(trphBounds,1),meanTrphDur); % prepare to pull 150 samples of data around each event
        trph_iM(:,1) = trphBounds(:,1);
        trph_iM = cumsum(trph_iM,2); % build indice matrix
        trphBounds(max(trph_iM,[],2)>length(signal),:) = [];% drop cycle which extend past the end of signal
        trph_iM(max(trph_iM,[],2)>length(signal),:) = []; % drop cycle which extend past the end of signal
        trphSig = signal(trph_iM);
        trphInd = trphBounds(:,1) + minInd(trphSig,[],2) - 1;
        
        % assemble ordered list of peaks and troughs
        peaksAndTrphs = [peakInd ones(size(peakInd)); trphInd -ones(size(trphInd))];
        [~,srtInds] = sort(peaksAndTrphs(:,1));
        peaksAndTrphs = peaksAndTrphs(srtInds,:);
        % remove double-troughs and double-peaks
        badInds = diff(peaksAndTrphs(:,2))==0;
        peaksAndTrphs(badInds,:) = [];
        
        % interpolate the intermediate points
        thetaPhs = nan(size(signal'));
        for i = 1:length(peaksAndTrphs)-1
            if peaksAndTrphs(i,2) == -1 % starting at trough
                thetaPhs(peaksAndTrphs(i,1):peaksAndTrphs(i+1,1)) = linspace(-pi,0,diff(peaksAndTrphs(i:i+1,1))+1);
            else %peaksAndTrphs(i,2) == +1 % starting at peak
                thetaPhs(peaksAndTrphs(i,1):peaksAndTrphs(i+1,1)) = linspace(0,pi,diff(peaksAndTrphs(i:i+1,1))+1);
            end
        end
        thetaAmp = ones(size(thetaPhs));
        band = mean(thetarange);
        %keyboard
        
    case 'wavelet'
        % returns a vector of values, spanning theta range, per time point
        if exist('band','var') & ~isempty(band)
            thetarange = band;
        else
            thetarange = 4:0.1:12;
        end
        [thetaPhs,thetaAmp] = multiphasevec(thetarange,signal,fs,6);
        band = thetarange;
        
    case 'hilbert'
        % returns a single phase estimate by time point using bandpass filtering
        if exist('band','var') & ~isempty(band)
            thetarange = minmax(band);
        else
            thetarange = [4 12];
        end
        
        hilbSig = hilbert(buttfilt(signal,thetarange,fs,'bandpass',4));
        thetaPhs = angle(hilbSig)';
        thetaAmp = abs(hilbSig)';
        band = mean(thetarange);
        
    otherwise
        error('Unknown method')
end