function display(obj,tag)

% Displays the results of the different timers

if ~isempty(obj.times)
    
    % Find the stopped timers with the correct tag 
    isRunning = cellfun(@isempty,obj.times(:,2));
    if nargin == 1
        maxTagLen = max(cellfun(@numel,obj.times(~isRunning,3)));
        isTarget = ~isRunning;
    else
        isTag = cellfun(@(a) strcmp(a,tag),obj.times(:,3));
        isTarget = ~isRunning & isTag;
        maxTagLen = numel(tag);
    end
    
    maxTagLen = max(maxTagLen,5) + 2;
    maxTimeLen = 13;
    
    % Display Title
    tagTitleStr = blanks(maxTagLen); tagTitleStr(1:3) = 'Tag';
    timeTitleStr = blanks(maxTimeLen); timeTitleStr(1:4) = 'Time';
    disp(repmat('-',1,maxTagLen+maxTimeLen));
    disp([tagTitleStr timeTitleStr])
    disp(repmat('-',1,maxTagLen+maxTimeLen));
    
    % Create strings to display
    cumTDelta = 0;
    for l=find(isTarget)'
        
        tDeltaStrDisp = blanks(maxTimeLen);
        tDelta = obj.times{l,4};
        tDeltaStr = secs2hms(tDelta);
        tDeltaStrDisp(1:numel(tDeltaStr)) = tDeltaStr;
        
        tagStrDisp = blanks(maxTagLen);
        tagStr = obj.times{l,3};
        tagStrDisp(1:numel(tagStr)) = tagStr;
        
        % Display
        if size(obj.times,2) == 5
            dataStrDisp = blanks(maxTimeLen);
            if isscalar(obj.times{l,5}) % The data is scalar
                dataStr = num2str(obj.times{l,5},'%.2f');
                dataStrDisp(1:numel(dataStr)) = dataStr;
            elseif isempty(obj.times{l,5}) % There is no data
                dataStrDisp(1) = '-';
            else % The data is an array or a cell array
                dataStrDisp(1) = '+';
            end
            disp([tagStrDisp tDeltaStrDisp dataStrDisp])
        else
            disp([tagStrDisp tDeltaStrDisp])
        end
        
        cumTDelta = cumTDelta + tDelta;
        
    end
    
    % Display the total
    if nargin > 1
        totalTitleStr = blanks(maxTagLen); totalTitleStr(1:5) = 'TOTAL';
        totalTimeStrDisp = blanks(maxTimeLen);
        totalTimeStr = secs2hms(cumTDelta);
        totalTimeStrDisp(1:numel(totalTimeStr)) = totalTimeStr;
        
        disp(repmat('-',1,maxTagLen+maxTimeLen));
        disp([totalTitleStr totalTimeStrDisp])
        disp(repmat('=',1,maxTagLen+maxTimeLen));
    end
    
    % Show how many timers are still running
    if nnz(isRunning) > 0
        fprintf('And %d timers are still running!\n',nnz(isRunning))
    end
else
    disp('Timing: Nothing to display!')
end

end
