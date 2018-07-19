function isRunningTag = running(obj)

% Find running timers

isRunning = cellfun(@isempty,obj.times(:,2));

if all(~isRunning)
    fprintf('Timing: No timers are running!\n');
    isRunningTag = {};
else
    % Return the tags of the running timers
    isRunningTag = obj.times(isRunning,3);
end

end
