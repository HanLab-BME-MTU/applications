function [ shift, tformshift ] = getChromaticShiftCrop( MD, maxShift )
%getChromaticShiftCrop Obtains the chromatic shift based upon a cropped
%region selected by the user

if(nargin < 2)
    maxShift = 3;
end

hmv = movieViewer(MD);
userData = get(hmv,'UserData');
movieFig = userData.getFigure('Movie');
movieAxes = findobj(movieFig,'Type','Axes');
hr = imrect(movieAxes);
wait(hr);
position = round(getPosition(hr));
delete(hr);
close(hmv);
delete(hmv);
drawnow;

largerPosition = position;
largerPosition(1:2) = largerPosition(1:2) - maxShift;
largerPosition(3:4) = largerPosition(3:4) + maxShift*2;

templateReader = CellReader(CropReader(MD.getReader(),position));
largerReader = CellReader(CropReader(MD.getReader(),largerPosition));

template = templateReader(1,1,maxShift:end-maxShift).to3D;
% control.im = templateReader(2,1,maxShift:end-maxShift).to3D;
reference = largerReader(2,1,:).to3D;

template = double(template);
reference = double(reference);

out = pdollar.images.normxcorrn(template,reference,'valid');
[~,outmaxi] = max(out(:));
% note row = y , col = x
[y,x,z] = ind2sub(size(out),outmaxi);


% control.out = pdollar.images.normxcorrn(double(control.im),double(reference),'valid');
% [~,control.outmaxi] = max(control.out(:));
% [control.y,control.x,control.z] = ind2sub(size(control.out),control.outmaxi);

% disp(control);
shift = [x y z] - [maxShift+1 maxShift+1 maxShift];
% shift = [x y z] - [control.x control.y control.z]

% keyboard;

if(nargout > 1)
    [optimizer, metric] = imregconfig('multimodal');
    tform = imregtform(template,reference,'translation',optimizer,metric);
    tformshift = tform.T(end,1:end-1) - [maxShift maxShift maxShift-1];
    
    % control.tform = imregtform(control.im,reference,'translation',optimizer,metric)
    % tformshift = tform.T(end,1:end-1) - control.tform.T(end,1:end-1);
    
end



end

