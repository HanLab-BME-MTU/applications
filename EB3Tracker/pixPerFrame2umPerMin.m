function [umPerMin] = pixPerFrame2umPerMin(pixPerFrame,secPerFrame,pixSizeNm)

umPerMin = (60/1000) * pixPerFrame * (pixSizeNm/secPerFrame);
