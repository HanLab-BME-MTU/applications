function responseAtOrder = getResponseAtOrderVec(response,Korg,Kg,freq)
%getResponseAtOrderVec Applies fft before running getResponeAtOrderVecHat
%
% See also orientationSpace.getResponseAtOrderVecHat

if(nargin < 4)
    freq = false;
end

if(~freq)
    response_hat = fft(response);
end
responseAtOrder = orientationSpace.getResponseAtOrderVecHat(response_hat,Korg,Kg);
if(~freq)
    responseAtOrder = ifft(responseAtOrder);
end


end

