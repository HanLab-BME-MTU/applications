function compressed = compressSteerable(steerableS)
    compressed.res = steerableS.res;
    compressed.nms = steerableS.nms ~= 0;
    compressed.theta = steerableS.theta(compressed.nms);
end
