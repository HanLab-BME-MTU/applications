function meshNoise(obj,height,length,nHorzLines,density,sigmaNoiseStart,sigmaNoiseEnd)
nVertLines = round(nHorzLines/(2*height)*length);
for i=1:nHorzLines
    xStart = 0;
    yStart = height - (i-1)*2*height/(nHorzLines-1);
    xEnd = length;
    yEnd = yStart;
    obj.lineWithDensity([xStart yStart 0],[xEnd yEnd 0],density);
end
for i=1:nVertLines
    xStart = (i-1)*length/(nVertLines-1);
    yStart = -height;
    xEnd = xStart;
    yEnd = height;
    obj.lineWithDensity([xStart yStart 0],[xEnd yEnd 0],density);
end
obj.addGaussianNoiseGradientX(sigmaNoiseStart,sigmaNoiseEnd,0,length)
end