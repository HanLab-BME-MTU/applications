function [ConnectivityMatrix,ConnectivityMatrixRussianForm,PositionVectices] = GetConnectivityPointbyPoint(matrix)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

[y,x] = find(matrix == 1);
PositionVectices     = [x,y];
Con   = zeros(length(x),length(x));
for it1 = 1:length(x)-1
    for it2 = it1+1:length(x)
        if (abs(x(it1)-x(it2)) == 1  & abs(y(it1)-y(it2)) == 0) |(abs(x(it1)-x(it2)) == 0  & abs(y(it1)-y(it2)) == 1)
            Con(it1,it2) = 1;
        elseif (abs(x(it1)-x(it2)) == 1  & abs(y(it1)-y(it2)) == 1)
            Con(it1,it2) = sqrt(2);
        end
    end
end
ConnectivityMatrix = Con + Con';
ConnectivityMatrixRussianForm=getRussianForm(ConnectivityMatrix);