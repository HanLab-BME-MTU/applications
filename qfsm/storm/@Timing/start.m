function start(obj,tag)

% Start a timer with a specified tag (optionally)

if nargin == 1
    tag = '';
end

tStart = toc(obj.ticID);

timeIdx = size(obj.times,1)+1;
obj.times(timeIdx,1) = {tStart};
obj.times(timeIdx,3) = {tag};

end
