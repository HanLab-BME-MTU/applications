function [ response ] = rescaleAngles( orientationMatrix, scaleFactor )
%orientationSpace.scaleDown Evaluate each orientation repsonse at different
%scale
%
% INPUT
% orientationMatrix : YxXxOrientation
% theta            : Orientation angles to interpolate
%
% OUTPUT
% y                 : YxXxTarget
    import orientationSpace.*;
    % 
    
    if(~isreal(orientationMatrix))
        imagOrientationMatrix = imag(orientationMatrix);
        imagOrientationMatrix = cat(3,imagOrientationMatrix,-imagOrientationMatrix);
        imagResponse = rescaleAngles(imagOrientationMatrix,scaleFactor);
        response = rescaleAngles(real(orientationMatrix),scaleFactor) + 1i*imagResponse(:,:,1:end/2);
        return;
    end
    
    s = size(orientationMatrix);
    M = reshape(orientationMatrix,s(1)*s(2),s(3));
    % todo, deal with even case
    n = (s(3))/2;
%     n = n+0.5;

    
%     x = [0:n -n:-1];
    x = wraparoundN(0:s(3)-1,-n,n);
    % toeplitz?
    xx = wraparoundN(bsxfun(@minus,x,(0:(s(3)-1))'),[-n n]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});
    tt = wraparoundN(bsxfun(@minus,wraparoundN((0:s(3)/scaleFactor-1)*scaleFactor,-n,n),(0:(s(3)-1))'),[-n n]);
    A = exp(-xx.^2/2);  
    T = exp(-tt.^2/2/(scaleFactor*scaleFactor));
    response = M*(A\T);
    
%     T = reshape(T,s(1)*s(2),1,[]);
    
%     response = sum(bsxfun(@times,M,T),2);

    % normalize the weights such that the columns sum to 1
%     w = w./repmat(sum(w),s(3),1);


    response = reshape(response,s(1),s(2),[]);



end

