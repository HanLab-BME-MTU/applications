function [pixPerFrame] = umPerMin2pixPerFrame(umPerMin,secPerFrame,pixSizeNm)

pixPerFrame = (1000/60) * umPerMin .* (secPerFrame/pixSizeNm);
