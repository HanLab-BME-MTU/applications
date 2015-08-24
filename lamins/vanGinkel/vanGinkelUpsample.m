function [ y ] = vanGinkelUpsample( orientationMatrix, target )
%vanGinkelUpsample upsample orientation information
%
% orientationMatrix : 
% target            :

    % 
    s = size(orientationMatrix);
    M = reshape(orientationMatrix,s(1)*s(2),s(3));
    % todo, deal with even case
    n = (s(3)-1)/2;

    
    x = [0:n -n:-1];
    xx = wraparound(bsxfun(@minus,x,(0:2*n)'),[-n n+1]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});

    A = 2*exp(-xx.^2);

    if(isscalar(target))
        target = 0:target:pi-target;
        target = target./(pi/s(3));
    end

    T = bsxfun(@minus,x,target')';
    T = wraparound(T,[-n n+1]);
    T = 2*exp(-T.^2);
    w = A\T;

    y = M*w;
    y = reshape(y,s(1),s(2),[]);

end