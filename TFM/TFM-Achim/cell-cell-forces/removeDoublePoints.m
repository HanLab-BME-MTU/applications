function [curve]=removeDoublePoints(curve)
    n=size(curve,1);
    i=1;
    while i<n
        if compPts(curve(i,:),curve(i+1,:))
            curve(i+1,:)=[];
            n=n-1;
        else
            i=i+1;
        end
    end
end