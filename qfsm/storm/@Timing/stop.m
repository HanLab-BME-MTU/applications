function tDelta = stop(obj,tag,data)

% Stop all the timers with the specified tag (optionally) and append data (optionally)

if nargin == 1 
    tag = '';
end

tStop = toc(obj.ticID);

isRunning = cellfun(@isempty,obj.times(:,2));
isTag = cellfun(@(a) strcmp(a,tag),obj.times(:,3));

isTarget = isRunning & isTag;

if all(~isTarget)
    fprintf('Timing: No timers with tag ''%s'' are running!\n',tag);
    tDelta = [];
else
    if nnz(isTarget) > 1
        fprintf('Timing: Multiple timers with tag ''%s'' stopped!\n',tag);
    end
    obj.times(isTarget,2) = {tStop};
    tStart = vertcat(obj.times{isTarget,1});
    tDelta = tStop-tStart;
    obj.times(isTarget,4) = num2cell(tDelta);
    if nargin == 3
        if ~iscell(data)
            data = {data};
        end
        obj.times(isTarget,5) = data;
    end
end

end
