function tDeltaStr = getFormatted(obj,tag)

% Returns the elapsed time of a specified timer in human readable form

if nargin == 1
    tDelta = obj.get();
else
    tDelta = obj.get(tag);
end

tDeltaStr = arrayfun(@secs2hms,tDelta,'UniformOutput',false);

end
