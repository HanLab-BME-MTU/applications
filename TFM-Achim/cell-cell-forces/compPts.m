function [result]=compPts(pt1,pt2,myEps)
if nargin < 3 || isempty(myEps)
    myEps=eps;
end
    if abs(pt1(1)-pt2(1))<=myEps && abs(pt1(2)-pt2(2))<=myEps
        result=true;
    else
        result=false;
    end
end