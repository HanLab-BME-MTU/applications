function [ response ] = interpolate( orientationMatrix, theta )
%orientationSpace.interpolate Evaluate each orientation repsonse at arbitrary angle
%at each position
%
% INPUT
% orientationMatrix : YxXxOrientation
% theta            : Orientation angles to interpolate
%
% OUTPUT
% y                 : YxXxTarget
    import orientationSpace.*;
    % 
    s = size(orientationMatrix);
    M = reshape(orientationMatrix,s(1)*s(2),s(3));
    % todo, deal with even case
    n = (s(3))/2;
%     n = n+0.5;

    
%     x = [0:n -n:-1];
    x = wraparoundN(0:s(3)-1,-n,n);
    % toeplitz?
    xx = wraparoundN(bsxfun(@minus,x,(0:(2*n-1))'),[-n n]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});

    A = exp(-xx.^2/2);
    M = M/A;
    
%     x = shiftdim(x,-1);
  
    % note we factor out the exp(1/2) constant factor
    
    theta = theta*s(3)/pi;
    theta = reshape(theta,s(1)*s(2),1,[]);

    T = bsxfun(@minus,x,theta);
    T = wraparoundN(T,[-n n]);
    T = exp(-T.^2/2);
%     T = reshape(T,s(1)*s(2),1,[]);
    
    response = sum(bsxfun(@times,M,T),2);

    % normalize the weights such that the columns sum to 1
%     w = w./repmat(sum(w),s(3),1);


    response = reshape(response,s(1),s(2),[]);



end

