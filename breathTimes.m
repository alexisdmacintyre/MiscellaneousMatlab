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

timeChunk = round(Fs/25); % Interval to check linearity of postive-going slope (i.e., inhalations) every ~40 ms
minPeakDist = 300; % Separation in ms between peaks

% Smooth
vector_smooth = movmean(vector,Fs*(win/1000));
vector_smooth = movmean(vector_smooth,round(Fs/50));
vector_smooth = rescale(vector_smooth,-1,1);

% Initial search for peaks (inhalation ends)
meanSlope = [];
c = 1;

for iii = 1:round((numel(vector_smooth)/timeChunk))-1
    
    t1_loop = timeChunk*iii;
    
    if iii == round((numel(vector_smooth)/timeChunk))-1
        t2_loop = numel(vector_smooth);
    else
        t2_loop = t1_loop+timeChunk;
    end
    
    vOut = mean(diff(vector_smooth(t1_loop:t2_loop)))*1e+4;
    
    if vOut < 1
        vOut(2) = 0;
        vOut(3) = 0;
        c = 1;
    else % Continuous positive-going slope
        vOut(2) = c;
        if c == 1
            vOut(3) = vOut(1);
        elseif c == 2
            vOut(3) = meanSlope(iii-1,1)+vOut(1);
        elseif c > 2
            vOut(3) = meanSlope(iii-1,3)+vOut(1);
        end
        c = c + 1;
    end
    meanSlope = [meanSlope ; vOut];
end

idx = zeros(size(meanSlope,1),1);
idx(meanSlope(:,2)>=floor(minDur/timeChunk)) = 1;
idx = strfind(idx',[1 0]);

peaks = [];

for ii = 1:numel(idx) % Localised peaks/ends of positive slopes
    if ii == 1
        t1 = 1;
        t2 = (timeChunk*idx(ii))+timeChunk;
    elseif ii == numel(idx)
        t1 = (timeChunk*idx(ii))-timeChunk;
        t2 = numel(vector_smooth);
    else
        t1 = (timeChunk*idx(ii))-timeChunk;
        t2 = (timeChunk*idx(ii))+timeChunk;
    end
    
    [~,m] = max(vector_smooth(t1:t2));
    peaks = [peaks ; t1+m];
end

peaks(peaks>=numel(vector_smooth)) = numel(vector_smooth);

% Within slopes, determine inhalation beginning/end points

onsets = [];
offsets = [];

for ii = 1:numel(peaks)
    
    if ii == 1
        t1 = 1;
        t2 = peaks(ii);
    else
        t1 = peaks(ii-1)+1;
        t2 = peaks(ii);
    end

    vec = vector_smooth(t1:t2);
    
    if numel(vec) > minDur

        % Analyse slope and trim flat segments preceding true onset/after
        % true peak
        
        meanSlope = [];
        c = 1;
        
        for iii = 1:floor((numel(vec)/timeChunk))
            
            t1_loop = timeChunk*iii;
            if iii == floor((numel(vec)/timeChunk))
                t2_loop = numel(vec);
            else
                t2_loop = t1_loop+timeChunk;
            end
            
            vOut = mean(diff(vec(t1_loop:t2_loop)))*1e+4;
            
            if vOut < 0
                vOut(2) = 0;
                vOut(3) = 0;
                c = 1;
            else
                vOut(2) = c;
                if c == 1
                    vOut(3) = vOut(1);
                elseif c == 2
                    vOut(3) = meanSlope(iii-1,1)+vOut(1);
                elseif c > 2
                    vOut(3) = meanSlope(iii-1,3)+vOut(1);
                end
                c = c + 1;
            end
            meanSlope = [meanSlope ; vOut];
        end

        % Narrow down end point (peak)
        
        maxTot = find(meanSlope(:,2)==max(meanSlope(:,2)),1,'last');
        pb = find(meanSlope(1:maxTot,2)==1,1,'last');
        
        if (pb+1)*timeChunk >= numel(vec)
            pb = 1;
        end
        
        mp = ceil((maxTot-pb)/2);
        maxR = quantile(meanSlope(pb+mp:maxTot,1),0.5);
        maxR = find(meanSlope(pb+mp:maxTot,1)>=maxR,1,'last');
        maxR = (pb+mp+maxR)-1;

        if maxTot-maxR > 2 & ...
                meanSlope(maxR,1)-meanSlope(maxTot,1) > ...
                0.5 * meanSlope(maxR,1)
            m = meanSlope(maxR,1);
            m = find(meanSlope(maxR:maxTot,1)<=0.5*m,1,'first');
            if isempty(m)
                m = round((maxTot-maxR)/2);
            end
            maxSlope = maxR+m;
        else
            maxSlope = maxTot;
        end

        % Specify inhale begin
        v = vec(timeChunk*pb:timeChunk*(pb+mp));
        p = quantile(v,0.1);
        p = find(v<=p,1,'last');
        onsets = [onsets ; t1+(timeChunk*pb)+p];
        
        % Specify inhale end
        if (timeChunk*maxSlope)+timeChunk < numel(vec)
            v = vec((timeChunk*maxSlope):(timeChunk*maxSlope)+timeChunk);
        else
            v = vec((timeChunk*maxSlope):end);
        end
        [~,p] = max(v);
        if isempty(p)
            p = 1;
        end
        offsets = [offsets ; t1+(timeChunk*maxSlope)+p];
        
    end

end

onsets(offsets>=numel(vector_smooth)) = [];
offsets(offsets>=numel(vector_smooth)) = [];

% If two peaks are close together, choose the first one (based on piloting,
% it seems that closely following peaks, even if larger, tend to be speech
% and not inhalation)

for i = 2:numel(offsets)
    t1 = offsets(i-1);
    if ~isnan(t1)
        if (offsets(i) - t1) <= minPeakDist
            t2 = offsets(i);
            m = max([vector_smooth(t1) vector_smooth(t2)]);
            thresh = 0.5*m;
            
            if vector_smooth(t1) >= thresh
                offsets(i) = NaN;
                onsets(i) = NaN;
            else
                offsets(i-1) = NaN;
                onsets(i) = NaN;
            end
        end
    end
end

onsets(isnan(onsets)) = [];
offsets(isnan(offsets)) = [];

% Join split breaths

x = [1 ; offsets];
y = [onsets ; numel(vector_smooth)];
z_dist = y-x;
z_ht = vector_smooth(y)-vector_smooth(x);

for i = 2:numel(x)-1
    if (z_dist(i) < 2*minDur) || (z_ht(i) > -0.05 & z_dist(i) < 1.5*minDur)

        old_t1 = x(i);
        old_t2 = y(i);
        new_t1 = onsets(i-1);
        new_t2 = offsets(i);
        
        if isempty(new_t1) || isnan(new_t1)
            j = i - 1;
            while isnan(new_t1) || isnan(new_t1)
                new_t1 = onsets(j);
                j = j-1;
            end
        end
        
        newSlope = vector_smooth(new_t1:new_t2);
        
        if ~(any(diff(find(newSlope<quantile(newSlope,0.45)))>timeChunk))
            
            % Check that the new longer breath does not contain a distinct
            % peak/valley within
            
            [~,m] = max(newSlope);
            idx = find(offsets==new_t2);
            offsets(offsets==old_t1) = NaN;
            if (new_t1 + m) - 1 ~= new_t2
                offsets(idx) = (new_t1 + m) - 1;
                x(x==new_t2) = (new_t1 + m) - 1;
            end

            [~,m] = min(newSlope);
            idx = find(onsets==new_t1);
            onsets(onsets==old_t2) = NaN;
            if (new_t1 + m) - 1 ~= new_t1
                onsets(idx) = (new_t1 + m) - 1;
                y(y==new_t1) = (new_t1 + m) - 1;
            end

        end
    end
end

onsets(isnan(onsets)) = [];
offsets(isnan(offsets)) = [];

brDur = offsets-onsets;
brHt = vector_smooth(offsets)-vector_smooth(onsets);
idx = brDur < minDur | brHt < minHt;
onsets(idx) = [];
offsets(idx) = [];

brDur = offsets-onsets;
brHt = vector_smooth(offsets)-vector_smooth(onsets);

% Remove breaths < 50% of median duration/height
if numel(onsets) > 4
    idx = brDur<0.5*median(brDur) | brHt<0.5*median(brHt);
    onsets(idx) = [];
    offsets(idx) = [];
end

if plotResults == 1
    figure;
    plot(vector_smooth,'LineWidth',0.75,'Color',[0.3010 0.7450 0.9330]);
    hold on
    plot(onsets,vector_smooth(onsets),'gx','MarkerSize',8)
    plot(offsets,vector_smooth(offsets),'rx','MarkerSize',8)
    yticks(-1:0.5:1)
    
    xlabel('Time (seconds)');
    xticks(1:round(numel(vector)/20):numel(vector)+1);
    t = numel(vector)/Fs;
    xticklabels(0:round(t/19):t);
    axis tight

    grid on
    set(gcf, 'Position', [500, 500, 900, 350])
    legend({'Breath Signal','Initiation','Maxima'},'Location','northeastoutside')
    title('Inhalation Events')
    set(gcf,'color','w')
end


end
