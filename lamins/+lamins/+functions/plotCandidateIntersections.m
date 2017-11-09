function [ handles ] = plotCandidateIntersections(candidates )
%plotCandidateIntersections plot output of findCandidateIntersections
% candidates is a 1D cell array with the lengths corresponding to angles
% in the range [0,pi)
N = length(candidates);

angles = 0:pi/N:pi-pi/N;
U = cos(angles);
V = sin(angles);

holdState = ishold;

hold on;

% p.x = xlim;
% p.y = ylim;
% I_size = [p.y(2) - p.y(1) p.x(2) - p.x(1) ];
I_size = size(get(imhandles(imgca),'CData'));
I_size = I_size([1 2]);

handles(N) = 0;

for i=1:N
    [Y,X] = ind2sub(I_size,candidates{i});
%     handles(i) = quiver(X-U(i),Y-V(i),repmat(U(i)*2,size(X)),repmat(V(i)*2,size(Y)),0,'ShowArrowHead','off');
    handles(i) = quiver(X-U(i),Y-V(i),repmat(U(i)*2,size(X)),repmat(V(i)*2,size(Y)),0,'Color','w');
end

hold off;

if(holdState)
    hold on;
end

keyboard;

end

