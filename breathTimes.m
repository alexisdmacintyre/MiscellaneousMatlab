% This function identifies respiratory events (inhalation beginnings and
% endings). Requires MATLAB 2016 or later.

% Required arguments: vector (Nx1 breath signal)

% Optional arguments: Fs (sample rate) 
% win (length of window in ms for smoothing the breath belt signal, 
% defaults to 250 ms/250 samples if no Fs given) 
% plotResults (0 or 1)

% Example usage: [begins,ends] = breathtimes(vector,Fs,win,plotResults) Fs
% = 1000; win = 250; [begins,ends] = breathTimes(breathBelt,1000,250,1)

% About: Likely inhalation endings are detected as local maxima in the 
% input breath belt vector.
% Likely beginnings of inhalation events are then determined based on
% peaks in the 2nd derivative of the breath belt signal. A few
% adjustments are performed to fix spuriously early points, etc. The true
% beginning of the breath is difficult to verify with most breath-belt type
% equipment, which are by nature a little shakey, but is at least more
% objective and consistent than human annotation. % If no Fs is given, 1 kHz
% is assumed for smoothing purposes

% If there are any missing beginning points (i.e., smaller count than count
% of peaks), they are removed, but you can also simply edit the script (see
% comments) to interpolate from successful detections.

% Troubleshooting: if the inhalation beginning are occuring too late in the
% breath belt positive going slope, try increasing the smoothing window

% Alexis Deighton MacIntyre
% a.macintyre.17@ucl.ac.uk


function [begins,ends] = breathTimes(vector,Fs,win,plotResults)

% Parse Inputs
if exist('win','var')
    if isempty(win)
        win = 250;
    end
else
    win = 250;
end

if exist('plotResults','var')
    if isempty(plotResults)
        plotResults = 0;
    end
else
    plotResults = 0;
end

if exist('Fs','var')
    if isempty(Fs)
        flagFs = 1;
        Fs = 1000;
    else
        flagFs = 0;
    end
else
    flagFs = 1;
    Fs = 1000;
end

% Smooth with moving window
vector_smooth = movmean(vector,Fs*(win/1000));
vector_smooth = rescale(vector_smooth,-1,1);

peaks = islocalmax(vector_smooth,'MinProminence',0.1,... % Detect peaks (inhale endings) in breath signal
    'FlatSelection','first');
peaks = find(peaks);

sigDeriv = [NaN(2,1) ; rescale(movmean(diff(diff(vector_smooth)),Fs/2)-0.25,-1,1)]; % use smoothed 2nd derivative to detect inhale beginnings
peakBases = islocalmax(sigDeriv,'MinProminence',0.1,'FlatSelection','first');
peakBases = find(peakBases);
peakHeights = sigDeriv(peakBases);
peaks = [peaks NaN(numel(peaks),1)];

for ii = 1:size(peaks,1) % Join inhale candidate beginning with the endings
    
    if ii ~= 1
        if ~isnan(peaks(ii-1)) % Don't include candidate beginnings prior to previous ending
            first = peaks(ii-1,1);
        else
            first = 0;
        end
    else
        first = 0;
    end
    c = peakBases(peakBases<peaks(ii,1) & peakBases>first);
    
    if isempty(c) % If no beginning, delete current ending. Alternatively, 
        % you can infer missing points using the median breath duration
        % value by commenting out this part and uncommenting below.
        peaks(ii,1) = NaN; 
    else
        if numel(c) == 1 % only one candidate beginning
            peaks(ii,2) = c; 
        else  % If multiple options, choose the beginning associated with highest value in 2nd derivative
            % d = peaks(ii,1)-c;
            % c = c(d>minBr); Uncomment and set minBr = to X if Fs known,
            % e.g., minBr = Fs*0.075;
            
            ch = peakHeights(peakBases<peaks(ii,1) & peakBases>first);
            [~,id] = max(ch);
            peaks(ii,2) = c(id);
        end
        
        % Adjust placement of beginning to make sure it's located at the
        % base of the breath belt slope (i.e., not spuriously earlier)
        checkLinear = vector_smooth(peaks(ii,2):peaks(ii,1));
        linDeriv = diff(checkLinear);
        linDeriv(linDeriv<=0) = NaN;
        linDeriv = isnan(linDeriv);
        linDeriv([1 numel(linDeriv)]) = 1;
        
        a = strfind(linDeriv',[0 1]);
        b = strfind(linDeriv',[1 0]);
        [~,id] = max(a-b);
        peaks(ii,2) = peaks(ii,2) + ...
            b(id) + 1;
    end
end

peaks(isnan(peaks(:,1)),:) = [];
    
% Uncomment to infer any missing beginning values using median breath duration
% brDur = peaks(:,1)-peaks(:,2);
% brDur = median(brDur(~isnan(brDur)));

for ii = find(isnan(peaks(:,2)))'
    peaks(ii,2) = peaks(ii,1)-brDur;
end

peaks(peaks(:,2)<0,2) = 1;

begins = peaks(:,2);
ends = peaks(:,1);

if plotResults == 1
    plot(vector_smooth,'LineWidth',0.75,'Color',[0.3010 0.7450 0.9330]);
    hold on
    plot(peaks(:,2),vector_smooth(peaks(:,2)),'gx','MarkerSize',8)
    plot(peaks(:,1),vector_smooth(peaks(:,1)),'rx','MarkerSize',8)
    
    yticks(-1:0.5:1)
    
    if flagFs == 0
        xlabel('Time (seconds)');
        xticks(1:round(numel(vector)/10):numel(vector)+1);
        t = numel(vector)/Fs;
        xticklabels(0:round(t/11):t);
    else
        xlabel('Samples');
        xticks(1:round(numel(vector)/11):numel(vector));
    end
    
    grid on
    set(gcf, 'Position', [500, 500, 900, 350])
    legend({'Breath Signal','Initiation','Maxima'},'Location','northeastoutside')
    title('Inhalation Events')
    set(gcf,'color','w')
end


end
