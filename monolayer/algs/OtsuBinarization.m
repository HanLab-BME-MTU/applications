% Thresholding-based binarization
% http://en.wikipedia.org/wiki/Otsu's_method
% http://www.mathworks.com/access/helpdesk/help/toolbox/images/index.html?/
% access/helpdesk/help/toolbox/images/graythresh.html&http://www.google.co.il/search?sourceid=navclient&ie=UTF-8&rlz=1T4SKPB_en&q=Otsu+matlab
% A - image
% factor - factor of the Otsu level
%
function [Abin, EM] = OtsuBinarization(A,factor)

if nargin < 2
    factor = 1;
end

debug = 0;

[level EM] = graythresh(A);
Abin = im2bw(A,level* factor);


if (debug)
    I = A;
    perim = dilate(bwperim(Abin),3);
    I(perim) = 255;
    figure; imagesc(I); colormap(gray); title('Otsu binarization');
end