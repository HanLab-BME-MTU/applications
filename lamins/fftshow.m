function [ h ] = fftshow( F , maglog, cm)
%fftshow shows the amplitude of a signal using intensity and the phase
%using color

if(nargin < 2 || isempty(maglog))
    maglog = false;
end
if(nargin < 3 || isempty(cm))
    cm = hsv(360);
end


N_colors = size(cm,1);

magnitude = abs(F);

if(islogical(maglog) && maglog)
    magnitude = log(magnitude+1);
elseif(isnumeric(maglog))
    magnitude = magnitude.^maglog;
end

magnitude = magnitude - min(magnitude(:));
magnitude = magnitude / max(magnitude(:));

phase = angle(F) + pi;
phase = round(phase / (2*pi) * (N_colors-1))+1;

I = reshape( repmat(magnitude(:),1,3) .* cm(phase(:),:) ,[size(magnitude) 3]);

h = imshow(I);


end

