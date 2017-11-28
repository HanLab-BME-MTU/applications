function responseAtOrder = getResponseAtOrderVec(response,Korg,Kg)
%getResponseAtOrderVec Applies fft before running getResponeAtOrderVecHat
%
% See also orientationSpace.getResponseAtOrderVecHat



response_hat = fft(response);
responseAtOrder = orientationSpace.getResponseAtOrderVecHat(response_hat,Korg,Kg);
responseAtOrder = ifft(responseAtOrder);


end

