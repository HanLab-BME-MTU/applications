function [pixPerFrame] = umPerMin2pixPerFrame(umPerMin,secPerFrame,pixSizeNm)
% converts microns/min -> pixels/frame given sec/frame and nm pixel size
%
% [pixPerFrame] = umPerMin2pixPerFrame(umPerMin,secPerFrame,pixSizeNm)
pixPerFrame = (1000/60) * umPerMin .* (secPerFrame/pixSizeNm);
