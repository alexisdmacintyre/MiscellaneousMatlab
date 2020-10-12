function [asynchronies,unmatched1,unmatched2,matched1,matched2] = calculateAsynchrony(vector1,vector2,plotMatches)

% Calculates asynchrony (absolute distance) between paired datapoints from
% two vectors. This function was developed for rhythm behavioural 
% research, but could be applied for other closest two-way match problems.

% For each observation, the closest two-way matched timepoints
% are found and their absolute asynchrony is calculated using euclidean
% distance. Remaining unpaired points are ignored, and can be specified
% by using optional output arguments that return their indices.

% Optional 2nd and 3rd outputs return a vector containing unpaired indices
% from vectors 1 and 2, respectively; 4th and 5th outputs are the paired
% equivalents.

% The method attempts to account for ordinality, in case a series of
% best pairings are consistently a little out of time (but not out of step)
% from each other.

% To visualise the pairings and make sure it's working right, use the 
% optional 3rd input argument (set flag to 1).

% Example usage including a visualisation:
% asynchronies = calculateAsynchrony(participant1Drum,participant2Drum,1)

% Example usage without visualisation, only returning paired and unpaired
% datapoints from participant 1:
% [asynchronies,unmatched1,~,matched1,~] = calculateAsynchrony(participant1Drum,participant2Drum)

% comments/suggations: alexisdharrison@gmail.com or a.macintyre.17@ucl.ac.uk

vector1 = sort(vector1);
vector2 = sort(vector2);

if ~exist('plotMatches','var') || plotMatches~=1
    plotMatches = 0;
end

if numel(vector1)>0 && numel(vector2)>0
    
    if size(vector1,1) == 1
        vector1 = vector1';
    end
    
    if size(vector2,1) == 1
        vector2 = vector2';
    end
    
    if numel(vector1) > numel(vector2) % Longer vector becomes vectorB
        vectorA = vector2;
        vectorB = vector1;
        
        longerVector = 1;
        
    elseif numel(vector1) < numel(vector2)
        vectorB = vector2;
        vectorA = vector1;
        
        longerVector = 2;
        
    else
        vectorA = vector1;
        vectorB = vector2;
        
        longerVector = 0;
    end
    
    if plotMatches == 1 % plot both vectors together
        close all
        scatter(vectorA,repmat(1,numel(vectorA),1))
        hold on
        scatter(vectorB,repmat(1.1,numel(vectorB),1))
        ylim([.95 1.15])
    end
    
    % Pair up observations to determine asynchronies between speakers
    counter = 1;
    asynchronies = [];
    numObs = numel(vectorB);
    
    for i = 1:numObs
        if sum(~isinf(vectorA))~=0 && sum(~isinf(vectorB))~=0
            
            maxLen = max(numel(vectorB),numel(vectorA));
            padding = inf(maxLen,1);
            vectorA = [vectorA ; padding(1:abs(numel(vectorA)-maxLen))];
            vectorB = [vectorB ; padding(1:abs(numel(vectorB)-maxLen))];
            
            [closestA,idxA] = min(abs(vectorA(i)-vectorB));
            [closestB,idxB] = min(abs(vectorB(i)-vectorA));
            [~,ind] = max([idxA idxB]);
            
            if idxA == idxB % vectorA(i)/vectorB(i) are each other's closest match
                
                asynchronies(counter) = closestA;
                            
                if plotMatches == 1
                    plot([vectorA(i),vectorB(idxA)], [1,1.1], '*-')
                end
                
                hits(1,counter) = vectorA(i);
                hits(2,counter) = vectorB(i);
                
                vectorA(i) = inf;
                vectorB(i) = inf;
                
            else % there are potentially unpaired points
                switch ind
                    case 1
                        lower = vectorB(vectorB<vectorB(idxA));
                        if numel(lower) > 1
                            vectorB(vectorB<lower(end)) = inf;
                        end
                    case 2
                        lower = vectorA(vectorA<vectorA(idxB));
                        if numel(lower) > 1
                            vectorA(vectorA<lower(end)) = inf;
                            
                        end
                end
                
                compA = vectorA(~isinf(vectorA));
                compB = vectorB(~isinf(vectorB));
                maxElement = min(numel(compA),numel(compB));
                
                switch ind
                    case 1
                        [~,diffAB] = min(abs(compA(1)-compB));
                    case 2
                        [~,diffAB] = min(abs(compB(1)-compA));
                end
                
                diffAB = diffAB - 1;
                
                if maxElement == 1
                    compDiff = 1;
                elseif maxElement < 8 + diffAB
                    compA = compA(1:maxElement);
                    compB = compB(1:maxElement);
                    compDiff = maxElement-diffAB;
                else
                    compA = compA(1:8+diffAB);
                    compB = compB(1:8+diffAB);
                    compDiff = 8;
                end
                
                switch ind
                    
                    case 1 % vectorB(i) may be unpaired
                        
                        if sum(abs(compA(1:compDiff)-compB(1:compDiff))) ...
                                > sum(abs(compA(1:compDiff) ...
                                - compB(1+diffAB:compDiff+diffAB)))
                            
                            asynchronies(counter) = closestA;
                            
                            if plotMatches == 1
                                plot([vectorA(i),vectorB(idxA)], [1,1.1], '*-')
                            end
                            
                            hits(1,counter) = vectorA(i);
                            hits(2,counter) = vectorB(idxA);
                            
                            vectorA(i) = inf;
                            vectorB(idxA) = inf;
                            vectorB(i:idxA-1) = [];
                            
                        else
                            
                            if vectorB(i) ~= inf
                                
                                asynchronies(counter) = abs(vectorA(i)-vectorB(i));

                                if plotMatches == 1
                                    plot([vectorA(i),vectorB(i)], [1,1.1], '*-')
                                end
                                
                                hits(1,counter) = vectorA(i);
                                hits(2,counter) = vectorB(i);
                                
                                vectorA(i) = inf;
                                vectorB(i) = inf;
                                
                            else
                                % vectorB(i) is inf now!
                                counter = counter - 1;
                                vectorA = [vectorA(1:i) ; vectorA(i:end)];
                                vectorA(i) = inf;
                            end
                            
                        end
                        
                    case 2 % vectorA(i) may be unpaired
                        
                        if sum(abs(compB(1:compDiff)-compA(1:compDiff))) ...
                                > sum(abs(compB(1:compDiff) ...
                                - compA(1+diffAB:compDiff+diffAB)))
                            
                            asynchronies(counter) = closestB;
                            
                            if plotMatches == 1
                                plot([vectorA(idxB),vectorB(i)], [1,1.1], '*-')
                            end
                            
                            hits(1,counter) = vectorA(idxB);
                            hits(2,counter) = vectorB(i);
                            
                            vectorB(i) = inf;
                            vectorA(idxB) = inf;
                            vectorA(i:idxB-1) = [];
                            
                        else
                            
                            if vectorA(i) ~= inf
                                asynchronies(counter) = abs(vectorB(i)-vectorA(i));
                                
                                if plotMatches == 1
                                    plot([vectorA(i),vectorB(i)], [1,1.1], '*-')
                                end
                                
                                hits(1,counter) = vectorA(i);
                                hits(2,counter) = vectorB(i);
                                
                                vectorA(i) = inf;
                                vectorB(i) = inf;
                                
                            else
                                counter = counter - 1;
                                vectorB = [vectorB(1:i) ; vectorB(i:end)];
                                vectorB(i) = inf;
                                
                            end
                            
                        end
                        
                end
            end
            
        end
        counter = counter + 1;
    end
    
    if longerVector == 1
        
        matched1 = find(ismember(vector1,hits(2,:)));
        matched2 = find(ismember(vector2,hits(1,:)));
        unmatched1 = find(~ismember(vector1,hits(2,:)));
        unmatched2 = find(~ismember(vector2,hits(1,:)));
        
    else
        
        matched1 = find(ismember(vector1,hits(1,:)));
        matched2 = find(ismember(vector2,hits(2,:)));
        unmatched1 = find(~ismember(vector1,hits(1,:)));
        unmatched2 = find(~ismember(vector2,hits(2,:)));
        
    end
    
else
    asynchronies = NaN;
    unmatched1 = NaN;
    unmatched2 = NaN;
    matched1 = NaN; 
    matched2 = NaN;
end

asynchronies = asynchronies';

if plotMatches == 1
    xlabel('Observations');
    xticks(1:10:max([vector1 ; vector2]));
    xticklabels([]);
    yticks([1 1.1]);
    if longerVector == 1
        yticklabels([2 1]);
    else
        yticklabels([1 2]);
    end
    grid on
    set(gcf, 'Position', [20, 200, 800, 350])
    title('Mutual Closest Matches')
    set(gcf,'color','w')
end


end
