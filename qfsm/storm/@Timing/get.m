function tDelta = get(obj,tag)

% Returns the elapsed time of a specified timer in seconds

if nargin == 1
    isTag = true(numel(obj.times(:,3)),1);
else
    isTag = cellfun(@(a) strcmp(a,tag),obj.times(:,3));
end

isRunning = cellfun(@isempty,obj.times(:,2));

isTarget = ~isRunning & isTag;

if all(~isTarget)
    fprintf('Timing: No stopped timers with tag ''%s'' found!\n',tag);
    tDelta = [];
else
    if nnz(isTarget) > 1
        if nargin == 1
            fprintf('Timing: Multiple stopped timers found!\n');
        else
            fprintf('Timing: Multiple stopped timers with tag ''%s'' found!\n',tag);
        end
    end
    
    tDelta = vertcat(obj.times{isTarget,4});
end

end
