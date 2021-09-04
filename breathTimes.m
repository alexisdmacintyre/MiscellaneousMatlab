% This function identifies inhalation events (inhalation beginnings and
% endings) from the breath belt signal. Requires MATLAB 2016 or later. 
%
% Usage: [onsets,offsets] =
% breathTimes(vector,Fs,'WinSz',20,'MinDur',100,'MinHeight',0.1,'Plot',1)
%
% Required arguments: vector (Nx1 breath signal); Fs (sample rate)
%
% Optional name pair arguments: 'WinSz' (for moving mean smoothing, default is 20 ms);
% 'MinDur' (minimum duration between onset and offset, default is 100 ms);
% 'MinHeight' (minimum vertical distance between onset and offset), default
% is 0.075 A.U.); 'Plot' (set to 1 to view results, default is 0)
%
% Suggestions: Some breath belt patterns look like inhalations but are
% associated with vocal exhalation, so make sure you corroborate this
% function's output with known speech onsets, etc.
%
% Alexis Deighton MacIntyre
% a.macintyre.17@ucl.ac.uk


function [onsets,offsets] = breathTimes(vector,Fs,varargin)

defaultWin = 20;
defaultDuration = 100;
defaultHeight = 0.075;
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

timeChunk = round(Fs/30); % Interval to check linearity of postive-going slope (i.e., inhalations) every ~33 ms
minPeakDist = 0.3*Fs; % Separation in ms between peaks
splitPeakDist = 0.15*Fs; % Threshold at which to join rather than delete breaths

% Smooth
vector_smooth = movmean(vector,Fs*(win/1000));
vector_smooth = rescale(vector_smooth,0,1);

% Initial round-up of peaks
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
    
    if vOut < 0.1
        vOut(2) = 0;
        vOut(3) = 0;
        c = 1;
    else % Is a continuous, positive-going slope
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
idx(meanSlope(:,2)>0) = 1;
j = strfind(idx',[1 0 1]);
idx(j+1) = 1;
j = strfind(idx',[0 1 0]);
idx(j+1) = 0;

idx_st = strfind(idx',[0 1]); 
idx_en = strfind(idx',[1 0]);

idx = idx_en'-idx_st';
ht = vector_smooth(idx_en);

% Cut some of the obvious noise-related peaks
for i = 2:numel(idx)
    if idx(i)*timeChunk < minDur & ...
            ht(i) < ht(i-1)
    idx_en(i) = NaN;
    end
end

idx_en(isnan(idx_en)) = [];

peaks = [];

for ii = 1:numel(idx_en) % Localised peaks/ends of positive slopes
    if ii == 1
        t1 = 1;
        t2 = (timeChunk*idx_en(ii))+timeChunk;
    elseif ii == numel(idx_en)
        t1 = (timeChunk*idx_en(ii))-timeChunk;
        t2 = numel(vector_smooth);
    else
        t1 = (timeChunk*idx_en(ii))-timeChunk;
        if (timeChunk*idx_en(ii))+timeChunk > numel(vector_smooth)
            t2 = numel(vector_smooth);
        else
            t2 = (timeChunk*idx_en(ii))+timeChunk;
        end
    end
    
    [~,m] = max(vector_smooth(t1:t2));
    peaks = [peaks ; t1+m];
    
end

peaks(peaks>=numel(vector_smooth)) = numel(vector_smooth);

peaks(isnan(peaks)) = [];


% Join split slopes

for i = 1:numel(peaks)-1
    if ~isnan(peaks(i))
        t1 = peaks(i);
        t2 = peaks(i+1);
        if t2-t1 < splitPeakDist
            if vector_smooth(t1) >= 0.75*vector_smooth(t2)
                peaks(i+1) = NaN;
            else
                peaks(i) = NaN;
            end
        end
    end
end

peaks(isnan(peaks)) = [];

% Within slopes, determine more precise inhalation beginning/end points

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

    t1_new = t1;
    t2_new = t2;
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
        
        % Join splits
        for j = 2:size(meanSlope,1)-1
            if meanSlope(j,2) == 0 & ...
                    meanSlope(j-1,2) ~= 0 & ...
                    meanSlope(j+1,2) ~= 0

                meanSlope(j,2) = meanSlope(j-1,2) + 1;                
                meanSlope(j,3) = meanSlope(j-1,3);
                jj = 1;
                while meanSlope(j+jj,2) ~= 0
                    
                    meanSlope(j+jj,2) = meanSlope(j+jj,2) + meanSlope(j,2);
                    meanSlope(j+jj,3) = meanSlope(j+jj,3) + meanSlope(j-1,3);
                    
                    jj = jj + 1;
                    if j+jj > size(meanSlope,1)
                        break
                    end
                end
            end
        end
        
        % Narrow down begin/end point windows
        
        maxTot = find(meanSlope(:,2)==max(meanSlope(:,2)),1,'last'); % end of slope
        pb = find(meanSlope(1:maxTot,2)==1,1,'last'); % beginning of slope
        
        if pb ~= 1
            t1_new = t1 + (pb*timeChunk);
            vec = vector_smooth(t1_new:t2);
        end
         meanSlope = meanSlope(pb:maxTot,:);

        rt = meanSlope(:,3)/max(meanSlope(:,3));
        
        % Check grade of slope
        
        if find(rt<0.06,1,'last') > 1 % trim any preceding flat section
            pb = find(rt<0.06,1,'last') - 1;
            meanSlope = meanSlope(pb:end,:);
            if pb ~= 1
                t1_new = t1_new + (pb*timeChunk);
                vec = vector_smooth(t1_new:t2);
            end
        end
        
        
        if find(rt>=0.95,1,'first') < numel(rt) % trim any subsequent flat section
            slopeMax = find(rt>=0.95,1,'first') + 1;
            if slopeMax ~= numel(rt)
                t2_new = t1_new + (slopeMax*timeChunk);
                vec = vector_smooth(t1_new:t2_new);
            end
        end

        % Specify inhale begin
        p = quantile(vec,0.01);
        p = find(vec<=p,1,'last');
        
        new_beg = t1_new + p;

        % Specify inhale end
        
        vec = vector_smooth(new_beg:t2_new);
        
        [~,p] = max(vec);
        new_en = new_beg + p - 1;
        
        if new_en - new_beg >= minDur
            onsets = [onsets ; new_beg];
            offsets = [offsets ; new_en];
        end
        
    end
    
end

onsets(offsets>=numel(vector_smooth)) = [];
offsets(offsets>=numel(vector_smooth)) = [];

% If two distinct peaks are close together, choose the first one (based on
% piloting, it seems that distinct, closely following peaks, even if
% larger, tend to be speech and not inhalation)

for i = 1:numel(offsets)-1
    t1 = offsets(i);
    t2 = offsets(i+1);
    if ~isnan(t1)
        if t2-t1 <= minPeakDist
            m = max([vector_smooth(t1) vector_smooth(t2)]);
            thresh = 0.5*m;
            
            if vector_smooth(t1) >= thresh
                onsets(i+1) = NaN;
                offsets(i+1) = NaN;
            else
                offsets(i) = NaN;
                onsets(i) = NaN;
            end
        end
    end
end

onsets(isnan(onsets)) = [];
offsets(isnan(offsets)) = [];

% Join any remaining split breaths
for ii = 1:numel(offsets)-1
    t0 = onsets(ii);
    if isnan(t0)
        j = ii;
        while isnan(t0)
            t0 = onsets(j);
            j = j - 1;
        end
    end
    t1 = offsets(ii);
    t2 = onsets(ii+1);
    rt = (vector_smooth(t1)-vector_smooth(t2))/...
        (vector_smooth(t1)-vector_smooth(t0));

    if (t2 - t1) < minDur & rt < 0.25
        onsets(ii+1) = NaN;
        offsets(ii) = NaN;
    end
end

onsets(isnan(onsets)) = [];
offsets(isnan(offsets)) = [];

% Remove sub-threshold breaths
z_ht = vector_smooth(offsets)-vector_smooth(onsets);
onsets(z_ht<minHt) = [];
offsets(z_ht<minHt) = [];

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
