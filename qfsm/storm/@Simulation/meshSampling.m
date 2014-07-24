function meshSampling(obj,height,length,nHorzLines,densityStart,densityEnd,sigmaNoise)
nVertLines = round(nHorzLines/(2*height)*length);
for i=1:nHorzLines
    xStart = 0;
    yStart = height - (i-1)*2*height/(nHorzLines-1);
    xEnd = length;
    yEnd = yStart;
    obj.lineWithDensityGradient([xStart yStart 0],[xEnd yEnd 0],densityStart,densityEnd);
end
for i=1:nVertLines
    xStart = (i-1)*length/(nVertLines-1);
    yStart = -height;
    xEnd = xStart;
    yEnd = height;
    density = densityStart+(densityEnd-densityStart)*(i-1)/(nVertLines-1);
    obj.lineWithDensity([xStart yStart 0],[xEnd yEnd 0],density);
end
obj.addGaussianNoise(sigmaNoise);
end