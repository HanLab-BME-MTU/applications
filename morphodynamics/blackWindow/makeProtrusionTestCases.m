function makeProtrusionTestCases(outFolder,caseName)
%This function is used to generate images of idealized examples for testing
%protrusion calculations and window propagation methods. 
%
%Hunter Elliott, 1/2011

imSize = [500,500];
nFrames = 100;

xVel = 1;

circRad = 50;

if nargin < 1 || isempty(outFolder)
    outFolder = uigetdir(pwd,'Select folder for output images:');
end

mkClrDir(outFolder)

%Create format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];

for j = 1:nFrames

    switch caseName



        case 'MovingCircle'
            
                if j == 1
                    circStrel = strel('disk',circRad,0);
                    currPos = [round(imSize(1)/2) circRad + 100];
                else
                    currPos = [currPos(1) currPos(2)+xVel];
                end                
                
                currMask = false(imSize);
                currMask(currPos(1),currPos(2)) = true;
                currMask = imdilate(currMask,circStrel);

                
        case 'Protrusion'
            
            if j == 1
               
                coefs = [linspace(0,-1,nFrames/2) linspace(-1,0,nFrames/2)];
                
            end
            
            x = -100:1:100;
            y = coefs(j)*100*(exp(-(x/25).^2))+200;
            
            xP = [200 400 300+x(end:-1:1)];
            yP = [400 400 y(end:-1:1)];
                
            currMask = poly2mask(xP,yP,imSize(1),imSize(2));
            
        case 'ProtrusionRetraction'
            
            if j == 1
               
                coefs = [linspace(0,-1,nFrames/4) ...
                         linspace(-1,1,nFrames/2) ...
                         linspace(1,0,nFrames/4)];
                
            end
            
            x = -100:1:100;
            y = coefs(j)*50*(exp(-(x/25).^2))+200;
            
            xP = [200 400 300+x(end:-1:1)];
            yP = [400 400 y(end:-1:1)];
                
            currMask = poly2mask(xP,yP,imSize(1),imSize(2));
            
        case 'Wave'
            
            if j == 1
                currPos = -50;
                coefs = [linspace(0,-1,nFrames/2) linspace(-1,0,nFrames/2)];                
            end
            
            x = -100:1:100;
            y = coefs(j)*50*(exp(-((x+currPos)/25).^2))+200;
            
            currPos = currPos + 1;
            
            xP = [200 400 300+x(end:-1:1)];
            yP = [400 400 y(end:-1:1)];
                
            currMask = poly2mask(xP,yP,imSize(1),imSize(2));
            
        case 'WaveConstantAmplitude'
            
            if j == 1
                currPos = -50;
                coefs = -1*ones(1,nFrames);
            end
            
            x = -100:1:100;
            y = coefs(j)*50*(exp(-((x+currPos)/25).^2))+200;
            
            currPos = currPos + 1;
            
            xP = [200 400 300+x(end:-1:1)];
            yP = [400 400 y(end:-1:1)];
                
            currMask = poly2mask(xP,yP,imSize(1),imSize(2));
            
        case 'StandingWave'
            
             if j == 1                
                currPos = -50;
                coefs = [linspace(0,-1,nFrames/2) linspace(-1,0,nFrames/2)];                
            end
            
            x = -100:1:100;
            y = coefs(j)*50*cos((x+currPos)*4*pi ./ 200)+200;                        
            
            xP = [200 400 300+x(end:-1:1)];
            yP = [400 400 y(end:-1:1)];
                
            currMask = poly2mask(xP,yP,imSize(1),imSize(2));
            %maskArea(j) = nnz(currMask);
            
        otherwise 

            error('Unrecognized method!')


    end
    
    fName = [outFolder filesep caseName '_' num2str(j,fString) '.tif'];
    currMask = uint16(bwdist(~currMask)*1000);
    imwrite(currMask,fName);
    
    
end

