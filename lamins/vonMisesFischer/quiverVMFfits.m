function [ handles ] = quiverVMFfits( idx, fits )
%quiverVMFfits Draws arrows at each idx for each fit found

holdState = ishold;

hold on;

% p.x = xlim;
% p.y = ylim;
% I_size = [p.y(2) - p.y(1) p.x(2) - p.x(1) ];
I_size = size(get(imhandles(imgca),'CData'));
I_size = I_size([1 2]);

[Y,X] = ind2sub(I_size,idx);

fit_lengths = cellfun('length',fits)/2 - 1;

max_fits = max(fit_lengths);

handles(max_fits) = 0;

for i=1:max_fits
    screen = fit_lengths > i - 1;
    [U,V] = cellfun(@(x) deal(x(i*2)*cos(x(1+(i-1)*2)),x(i*2)*sin(x(1+(i-1)*2))),fits(screen));
    handles(i) = quiver(X(screen)-U,Y(screen)-V,U*2,V*2,0,'ShowArrowHead','off');
end

hold off;

if(holdState)
    hold on;
end


end

