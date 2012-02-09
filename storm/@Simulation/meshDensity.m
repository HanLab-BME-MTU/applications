function meshDensity(obj,height,length,nHorzLines,squeezeFactor,density,sigmaNoise)
slope = (2*height/(nHorzLines-1)/squeezeFactor-2*height/(nHorzLines-1))/length;
for i=1:nHorzLines
    xStart = 0;
    yStart = height - (i-1)*2*height/(nHorzLines-1);
    xEnd = length;
    yEnd = yStart/squeezeFactor;
    obj.lineWithDensity([xStart yStart 0],[xEnd yEnd 0],density);
end
x = 0;
while x < length
    obj.lineWithDensity([x height+(nHorzLines-1)/2*slope*x 0],[x -(height+(nHorzLines-1)/2*slope*x) 0],density);
    x = x + (2*height/(nHorzLines-1) + slope*x);
end
obj.addGaussianNoise(sigmaNoise);
end