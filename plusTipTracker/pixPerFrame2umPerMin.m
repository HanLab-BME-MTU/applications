function [umPerMin] = pixPerFrame2umPerMin(pixPerFrame,secPerFrame,pixSizeNm)
% converts pixels/frame -> microns/min given sec/frame and nm pixel size
%
% [umPerMin] = pixPerFrame2umPerMin(pixPerFrame,secPerFrame,pixSizeNm)
umPerMin = (60/1000) * pixPerFrame * (pixSizeNm/secPerFrame);
