function [relMousePos] = getRelativeMousePos(hAxes)

    mousePos = get(hAxes,'CurrentPoint');  % The current point w.r.t the figure.
    xylim = [get(hAxes, 'xlim'); get(hAxes, 'ylim')];
    
    relMousePos =(mousePos(1,1:2) - xylim(1:2)) ./ (xylim(3:4) - xylim(1:2));

end