% This function identifies inhalation events (inhalation beginnings and
% endings) from the breath belt signal. Requires MATLAB 2016 or later. 
%
% Usage: [onsets,offsets] =
% breathTimes(vector,Fs,'WinSz',200,'MinDur',100,'MinHeight',0.1,'Plot',1)
%
% Required arguments: vector (Nx1 breath signal); Fs (sample rate)
%
% Optional name pair arguments: 'WinSz' (for smoothing, default is 250 ms);
% 'MinDur' (minimum duration between onset and offset, default is 150 ms);
% 'MinHeight' (minimum vertical distance between onset and offset), default
% is 0.05 A.U.); 'Plot' (set to 1 to view results, default is 0)
%
% Suggestions: Some breath belt patterns look like inhalations but are
% associated with vocal exhalation, so make sure you corroborate this
% function's output with known speech onsets, etc.
%
% Alexis Deighton MacIntyre
% a.macintyre.17@ucl.ac.uk


function [onsets,offsets] = breathTimes(vector,Fs,varargin)

defaultWin = 250;
defaultDuration = 150;
defaultHeight = 0.5;
defaultPlot = 0;

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validPlot = @(x) (x == 0) || (x == 1);
addRequired(p,'vector');
addRequired(p,'Fs',validScalarPosNum);
addParameter(p,'WinSz',defaultWin,validScalarPosNum);
addParameter(p,'MinDur',defaultDuration,validScalarPosNum);
addParameter(p,'MinHeight',defaultHeight,validScalarPosNum);
addParameter(p,'Plot',defaultPlot,validPlot);

parse(p,vector,Fs,varargin{:});

win = p.Results.WinSz;
minDur = (p.Results.MinDur/1000)*Fs;
minHt = p.Results.MinHeight;
plotResults = p.Results.Plot;

timeChunk = round(Fs/20); % Checks linearity of postive-going slope (i.e., inhalations) every 50 ms

% Smooth
vector_smooth = movmean(vector,Fs*(win/1000));
vector_smooth = movmean(vector_smooth,round(Fs/50));
vector_smooth = rescale(vector_smooth,-1,1);

peaks = islocalmax(vector_smooth,'MinProminence',0.1,... % Detect peaks (inhale endings) in breath signal
    'FlatSelection','first');
peaks = find(peaks);

onsets = [];
offsets = [];

for ii = 1:size(peaks,1)
    
    % Find peak base
    if ii == 1
        t1 = 1;
    else
        t1 = peaks(ii-1)-minDur;
    end
    t2 = peaks(ii);
    vec = vector_smooth(t1:t2);
    
    % Check for positive linearity
    meanSlope = [];
    c = 1;
    
    for iii = 1:round((numel(vec)/timeChunk))-1
        t1_loop = timeChunk*iii;
        if iii == round((numel(vec)/timeChunk))-1
            t2_loop = numel(vec);
        else
            t2_loop = t1_loop+timeChunk;
        end
        
        vOut = mean(diff(vec(t1_loop:t2_loop)))*1e+4;
        
        if vOut < 0.5
            vOut(2) = 0;
            c = 1;
        else
            vOut(2) = c;
            c = c + 1;
        end
        
        meanSlope = [meanSlope ; vOut];
    end
    
    [~,maxSlope] = max(meanSlope(:,2));
    peakBase = find(meanSlope(1:maxSlope,2)==1,1,'last');
    
    % Adjust inhale end
    if (timeChunk*maxSlope)+timeChunk < numel(vec)
        v = vec((timeChunk*peakBase):(timeChunk*maxSlope)+timeChunk);
    else
        v = vec((timeChunk*peakBase):numel(vec));
    end
    p = quantile(v,0.99);
    p = find(v>=p,1,'first');
    offsets = [offsets ; t1+(timeChunk*peakBase)+p];
    
    % Adjust inhale begin
    if (timeChunk*peakBase)-timeChunk > 1
        v = vec((timeChunk*peakBase)-timeChunk:(timeChunk*peakBase)+timeChunk);
    else
        v = vec(1:(timeChunk*peakBase)+timeChunk);
    end
    p = quantile(v,0.01);
    p = find(v<=p,1,'last');
    onsets = [onsets ; t1+(timeChunk*peakBase)-timeChunk+p];

end

brDur = offsets-onsets;
brHt = vector_smooth(offsets)-vector_smooth(onsets);
idx = brDur < minDur | brHt < minHt;
onsets(idx) = [];
offsets(idx) = [];

if plotResults == 1
    plot(vector_smooth,'LineWidth',0.75,'Color',[0.3010 0.7450 0.9330]);
    hold on
    plot(onsets,vector_smooth(onsets),'gx','MarkerSize',8)
    plot(offsets,vector_smooth(offsets),'rx','MarkerSize',8)
    
    yticks(-1:0.5:1)
    
    xlabel('Time (seconds)');
    xticks(1:round(numel(vector)/10):numel(vector)+1);
    t = numel(vector)/Fs;
    xticklabels(0:round(t/11):t);

    grid on
    set(gcf, 'Position', [500, 500, 900, 350])
    legend({'Breath Signal','Initiation','Maxima'},'Location','northeastoutside')
    title('Inhalation Events')
    set(gcf,'color','w')
end


end