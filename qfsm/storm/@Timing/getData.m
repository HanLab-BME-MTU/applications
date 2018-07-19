function data = getData(obj,tag)

% Returns the data that was attached to a specified timer

if nargin == 1
    isTag = true(numel(obj.times(:,3)),1);
else
    isTag = cellfun(@(a) strcmp(a,tag),obj.times(:,3));
end

isRunning = cellfun(@isempty,obj.times(:,2));
isTarget = ~isRunning & isTag;

if all(~isTarget)
    fprintf('Timing: No stopped timers with tag ''%s'' found!\n',tag);
    data = {};
else
    if nnz(isTarget) > 1
        if nargin == 1
            fprintf('Timing: Multiple stopped timers found!\n');
        else
            fprintf('Timing: Multiple stopped timers with tag ''%s'' found!\n',tag);
        end
    end
    
    data = obj.times(isTarget,5);
end

end
