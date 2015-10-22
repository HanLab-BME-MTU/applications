function [ y, yy, varargout ] = vanGinkelDerivative( orientationMatrix, target )
%vanGinkelUpsample upsample orientation information
%
% INPUT
% orientationMatrix : YxXxOrientation
% target            : Orientation angles to estimate
%
% OUTPUT
% y                 : YxXxTarget

    import vanGinkel.*;

    % 
    s = size(orientationMatrix);
    M = reshape(orientationMatrix,s(1)*s(2),s(3));
    % todo, deal with even case
    n = (s(3))/2;
%     n = n+0.5;

    if(nargin < 2)
        target = pi/s(3);
    end

    
%     x = [0:n -n:-1];
    x = wraparoundN(0:s(3)-1,-n,n);
    % toeplitz?
    xx = wraparoundN(bsxfun(@minus,x,(0:(2*n-1))'),[-n n]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});

    A = exp(-xx.^2/2);

    if(isscalar(target))
        target = 0:target:pi-target;
    end
    target = target./(pi/s(3));
    

    T = bsxfun(@minus,x,target')';
    T = wraparoundN(T,[-n n]);
    eT = exp(-T.^2/2);
    % derivative with respect to target
    w = (eT.*T);
    % normalize the weights such that the columns sum to 1
%     w = w./repmat(sum(w),s(3),1);
    M = M/A;
    y = M*w;
    
    % new size 3
    ns3 = size(y,2);
    
    if(nargout > 1)
        w2 = w.*T - eT;
        yy = M*(w2);
        yy = reshape(yy,s(1),s(2),ns3);
    end
    if(nargout > 2)
        % calculate higher derivatives recursively
%         w = {w w2};
%         w{nargout} = [];
        varargout{nargout-2} = [];
        for i = 3:nargout
            wn = -(i-1)*w + T.*w2;
            varargout{i-2} = reshape(M*wn,s(1),s(2),ns3);
            w = w2;
            w2 = wn;
        end
    end
    
    y = reshape(y,s(1),s(2),ns3);


end