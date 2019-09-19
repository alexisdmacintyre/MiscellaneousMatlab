function asynchronies = calculateAsynchrony(vector1,vector2,plotMatches)

% Calculates asynchrony (absolute distance) between paired datapoints from
% two vectors. This function was developed for rhythm behavioural 
% research, but could be applied for other closest two-way match problems.

% For each observation, the closest two-way matched timepoints
% are found and their absolute asynchrony is calculated using euclidean
% distance. Single unpaired points are ignored, i.e., there is no 'penalty'
% or note for unpaired points, so adjust accordingly.

% The method attempts to account for ordinality, in case a series of
% best pairings are consistently a little out of time (but not out of step)
% from each other.

% To visualise the pairings and make sure it's working right, use the 
% optional 3rd argument (1).

% Example usage including a visualisation:
% asynchronies = calculateAsynchrony(participant1Drum,participant2Drum,1)

if ~exist('plotMatches','var') || plotMatches~=1
    plotMatches = 0;
end

if numel(vector1) > numel(vector2) % Longer vector becomes vectorB
    vectorA = vector2;
    vectorB = vector1;
elseif numel(vector1) < numel(vector2)
    vectorB = vector2;
    vectorA = vector1;
else
    vectorA = vector1;
    vectorB = vector2;
end

if plotMatches == 1 % plot both vectors together
    close all
    scatter(vectorA,repmat(1,numel(vectorA),1))
    hold on
    scatter(vectorB,repmat(1.1,numel(vectorB),1))
    ylim([.95 1.15])
end

% Pair up observations
counter = 1;
asynchronies = [];
numObs = numel(vectorB);

for i = 1:numObs
    if sum(~isinf(vectorA))~=0 && sum(~isinf(vectorB))~=0 % stop when one vector runs out
        
        maxLen = max(numel(vectorB),numel(vectorA));
        padding = inf(maxLen,1);
        vectorA = [vectorA ; padding(1:abs(numel(vectorA)-maxLen))];
        vectorB = [vectorB ; padding(1:abs(numel(vectorB)-maxLen))];
        
        [closestA,idxA] = min(abs(vectorA(i)-vectorB));
        [closestB,idxB] = min(abs(vectorB(i)-vectorA));
        [~,idx] = max([idxA idxB]);
        
        if idxA == idxB % vectorA(i)/vectorB(i) are each other's closest match
            
            asynchronies(counter) = closestA;
            
            if plotMatches == 1
                plot([vectorA(i),vectorB(i)], [1,1.1], '*-')
            end
            
            vectorA(i) = inf;
            vectorB(i) = inf;
            
        else % there are potentially unpaired points
            switch idx
                case 1 % vectorA(i) is closer to vectorB(i+X)
                    lower = vectorB(vectorB<vectorB(idxA));
                    if numel(lower) > 1
                        vectorB(vectorB<lower(end)) = inf;
                    end
                case 2 % vectorB(i) is closer to vectorA(i+X)
                    lower = vectorA(vectorA<vectorA(idxB));
                    if numel(lower) > 1
                        vectorA(vectorA<lower(end)) = inf;
                    end
            end
            
            compA = vectorA(~isinf(vectorA));
            compB = vectorB(~isinf(vectorB));
            maxElement = min(numel(compA),numel(compB));
            
            switch idx
                case 1
                    [~,diffAB] = min(abs(compA(1)-compB));
                case 2
                    [~,diffAB] = min(abs(compB(1)-compA));
            end
            
            diffAB = diffAB - 1;
            
            % Check if future pairings (i+7) are associated with a smaller
            % error before accepting skip and leaving point(s) unpaired
            
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
            
            switch idx
                
                case 1 % check if vectorB(i) should be unpaired
                    
                    if sum(abs(compA(1:compDiff)-compB(1:compDiff))) ...
                            > sum(abs(compA(1:compDiff) ...
                            - compB(1+diffAB:compDiff+diffAB)))
                        
                        asynchronies(counter) = closestA;
                        
                        if plotMatches == 1
                            plot([vectorA(i),vectorB(idxA)], [1,1.1], '*-')
                        end
                        
                        vectorA(i) = inf;
                        vectorB(idxA) = inf;
                        vectorB(i:idxA-1) = [];
                        
                    else % skipping doesn't reduce future error, only
                        % skip if there are more than one unpaired points
                        % before the next match
                        
                        if vectorB(i) ~= inf
                            asynchronies(counter) = abs(vectorA(i)-vectorB(i));
                            
                            if plotMatches == 1
                                plot([vectorA(i),vectorB(i)], [1,1.1], '*-')
                            end
                            
                            vectorA(i) = inf;
                            vectorB(i) = inf;
                        else
                            counter = counter - 1;
                            vectorA = [vectorA(1:i) ; vectorA(i:end)];
                            vectorA(i) = inf;
                        end
                        
                    end
                    
                case 2 % check if vectorA(i) should be unpaired
                    
                    if sum(abs(compB(1:compDiff)-compA(1:compDiff))) ...
                            > sum(abs(compB(1:compDiff) ...
                            - compA(1+diffAB:compDiff+diffAB)))
                        
                        asynchronies(counter) = closestB;
                        
                        if plotMatches == 1
                            plot([vectorA(idxB),vectorB(i)], [1,1.1], '*-')
                        end
                        
                        vectorB(i) = inf;
                        vectorA(idxB) = inf;
                        vectorA(i:idxB-1) = [];
                        
                    else % skipping doesn't reduce future error, only
                        % skip if there are more than one unpaired points
                        % before the next match
                        
                        if vectorA(i) ~= inf
                            asynchronies(counter) = abs(vectorB(i)-vectorA(i));
                            
                            if plotMatches == 1
                                plot([vectorA(i),vectorB(i)], [1,1.1], '*-')
                            end
                            
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

set(gca,'YTick', [])

end