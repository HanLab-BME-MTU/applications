function [ Ihat_skel ] = showSkelOnIntensity( I, skel, c )
%showSkelOnIntensity Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    c = [1 0 0];
end

Ia = imadjust(I,[0.2 0.8]);
Ihat = imsubtract(imadd(Ia,imtophat(Ia,strel('disk',5))),imbothat(Ia,strel('disk',5)));
Ihat_skel = repmat(Ihat,[1 1 3]);
Ihat_skel = shiftdim (Ihat_skel,2);
Ihat_skel(1,skel) = c(1);
Ihat_skel(2,skel) = c(2);
Ihat_skel(3,skel) = c(3);
Ihat_skel = shiftdim (Ihat_skel,1);

imshow(Ihat_skel);


end

