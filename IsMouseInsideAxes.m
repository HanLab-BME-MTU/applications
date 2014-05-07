function [isInside] = IsMouseInsideAxes(hAxes)

    mousePos = get(hAxes,'CurrentPoint');  % The current point w.r.t the figure.
    xylim = [get(hAxes, 'xlim'); get(hAxes, 'ylim')];
    
    if any( mousePos(1,1:2) < xylim(1:2) ) || any( mousePos(1,1:2) > xylim(3:4) )
        isInside = false;
    else
        isInside = true;
    end
    
end