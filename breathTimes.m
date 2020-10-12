function [beginBreath,endBreath] = breathTimes(respirationSignal,cleanUp,Fs,winSz,prominence,plotResults)

% This function identifies local maxima (probable end-of-inhalation events)
% in the input vector, after optionally smoothing, detrending, and
% rescaling. Note that it's optimised for use with speech breathing, which 
% has a distinct and asymmetric profile.

% The peak-finding parameter here is topographical prominence (how much
% each peak stands out from its neighbours). The value of 0.25 is the
% default and it's what worked best for us based on piloting an in-house
% dataset, so go ahead and adjust for better results.

% If you supply a sample rate (Fs), the smoothing is based on a default 50 ms
% window. Otherwise, the window is defined by the total number of 
% samples / 500. You can optionally specify a different window length in
% samples, which will override the default ms even if you're given an Fs.

% Once the maxima are determined, a signal change-based algorithm
% searches for the likely beginnings of inhalation events. The true
% beginning of the breath is difficult to verify with most breath-belt type
% equipment, which are by nature a little shakey, but is at least more
% objective and consistent than human annotation. This part of the function
% is made with thanks to Jonas from the MATLAB Answers community.

% If there are any missing beginning points (i.e., smaller count than count
% of peaks), they are interpolated from successful detections.

% Finally, the time points for the inhalation initiations and maxima are
% returned. You can optionally plot your vector alongside the time points.

% If you provide the Fs, arguments out and plots are in milliseconds.

% Examples:

% 1. Default options are no clean up, window is total samples/500, 
% no plotting, sample rate is unknown, default prominence = 0.25

% [beginBreath,endBreath] = breathTimes(breathData);

% 2. Plot, but no clean up, sample rate is unknown and default prominence,
% window for smoothing is 1000 samples

% [beginBreath,endBreath] = breathTimes(breathData,0,[],1000,[],1);

% 3. Clean up and plot with a 16 kHz sample rate, prominence of 0.5,
% default smoothing of 50 ms

% [beginBreath,endBreath] = breathTimes(breathData,1,16000,[],0.5,1);

% Troubleshooting:

% You can adjust the prominence setting manually if too many or too few
% peaks are found

% You can also smooth more or less by setting a larger or smaller window
% size, which can help with noisy or erroneous initiations in particular

% Feel free to contact me with suggestions!
% Alexis Deighton MacIntyre a.macintyre.17@ucl.ac.uk

if exist('cleanUp','var')
    if isempty(cleanUp) || cleanUp ~= 1
        cleanUp = 0;
    end
else
    cleanUp = 0;
end

if exist('winSz','var')
    if isempty(winSz)
        winSz = 0;
    end
else
    winSz = 0;
end

if exist('plotResults','var')
    if isempty(plotResults) || plotResults ~= 1
        plotResults = 0;
    end
else
    plotResults = 0;
end

if ~exist('Fs','var')
    Fs = 0;
elseif isempty(Fs)
    Fs = 0;
else % Resample to 1 kHz
    [N,D] = rat(1000/Fs);
    respirationSignal = resample(respirationSignal,N,D); 
end

if exist('prominence','var')
    if isempty(prominence)
        prominence = 0.25;
    end
else
    prominence = 0.25;
end


% Optional Cleaning Up

if cleanUp == 1
    respirationSignal = detrend(respirationSignal); % detrend
    respirationSignal = rescale(respirationSignal,-1,1); % rescale to -1,1
    if Fs > 0 && winSz == 0 % smooth using a moving mean at 50 ms
        win = round(.05*Fs);
        respirationSignal = movmean(respirationSignal,win);
    elseif Fs == 0 && winSz == 0 % smooth using a moving mean of count samples/500
        respirationSignal = movmean(respirationSignal,round(numel(respirationSignal)/500));
    elseif winSz ~= 0 % smooth using a moving mean based on directed samples
        respirationSignal = movmean(respirationSignal,winSz);        
    end
end

% Detect peaks (end of breath events/begin of exhalation)

peaks = islocalmax(respirationSignal,'MinProminence',prominence,...
    'FlatSelection','first');
peaks = find(peaks);

%% Detect slope bases (begin of breath events/end of exhalation)

minThresh = 100;
candidates = [];

while numel(candidates)<numel(peaks) % approximate count of bases according to count of peaks
    minThresh = minThresh-1;
    if minThresh < 0
        idx = findchangepts(respirationSignal,'Statistic','linear','MinThreshold',0);
        dx = diff(idx);
        dy = diff(respirationSignal(idx));
        slope = dy./dx;
        z = find(slope > 0.0001);
        candidates = idx(z);
        break
    end
    idx = findchangepts(respirationSignal,'Statistic','linear','MinThreshold',minThresh);
    dx = diff(idx);
    dy = diff(respirationSignal(idx));
    slope = dy./dx;
    z = find(slope > 0.0001);
    candidates = idx(z);
end

slopeBases = []; % Select only the slope bases that immediately precede a peak
for i = 1:numel(peaks)
    currentCandidates = candidates((peaks(i)-candidates)>0);
    [~,idx] = min(peaks(i)-currentCandidates);
    
    if ~isempty(idx) & sum(slopeBases==currentCandidates(idx))==0
        slopeBases(i) = peaks(i) - ...
            (peaks(i)-currentCandidates(idx));
    else % no slopeBase found, or is already taken
        slopeBases(i) = NaN; % infer this value later by taking median value of other slope bases
    end
end

idx = isnan(slopeBases); % Infer any missing slopeBases
if ~isempty(idx)
    slopeBases(idx) = (peaks(idx)-nanmean(peaks-slopeBases'));
end

for i = 1:numel(slopeBases)-1 % Remove any duplicate slopeBases
    if slopeBases(i+1)-slopeBases(i) < peaks(i)-slopeBases(i)
        slopeBases(i) = NaN;
    end
end

slopeBases(isnan(slopeBases)) = [];

beginBreath = slopeBases';
endBreath = peaks;

if isempty(beginBreath)
    error('No breathing events found - try adjusting prominence downward?');
end

%% Plot, if desired

if plotResults == 1
    figure;
    plot(respirationSignal,'LineWidth',0.75,'Color',[0.3010 0.7450 0.9330]);
    %grid on
    hold on; 
    
    plot(beginBreath,respirationSignal(round(beginBreath)),'go','MarkerSize',8)
    plot(endBreath,respirationSignal(round(endBreath)),'ro','MarkerSize',8)
    
    if Fs > 0
        xlabel('Time (seconds)');
        xticks(1:10000:numel(respirationSignal));
        xticklabels(0:10:(numel(respirationSignal)-1)/1000)
    else
        xlabel('Samples');
    end
    
    grid on
    set(gcf, 'Position', [20, 200, 800, 350])
    legend({'Breath Signal','Initiation','Maxima'},'Location','northeastoutside') 
    title('Respiratory Initiation and Maxima Events')
    set(gcf,'color','w')
end
end