function [backLevelFirst, backLevelLast, nucLevelFirst, nucLevelLast, haloLevelFirst, haloLevelLast] = ptGetThresholds (ptJob, jobNumber)
% ptGetThresholds allows the user to specify the intensity values of his images interactively. Values of interest are: background, nuclei, halos.
% each value get's specified for the first and the last frame of the current analysis, so that for every frame the values can be interpolated
%
% SYNOPSIS       ptGetThresholds (ptJob, jobNumber)
%
% INPUT          ptJob : job structure containing all the info on the current job
%                jobNumber : current job number
%
% OUTPUT         backLevelFirst : level of the background in the first image
%                backLevekLast : level of the background in the last image
%                nucLevelFirst : level of the nuclei in the first image
%                nucLevelLast : level of the nuclei in the last image
%                haloLevelFirst : level of the halos in the first image
%                haloLevelLast : level of the halos in the last image
%
% DEPENDENCIES   ptGetThresholds uses { nothing }
%                                  
%                ptGetThresholds is used by { PolyTrack }
%
% REMARK         Details on the ptJob structure (partly):
% 		    imagedirectory : where the images are
% 		    imagenameslist : list of images within imagedirectory with imagename
% 		    firstimage     : which images shall we start with (refers to imagenameslist)
% 		    lastimage      : which images shall be the last one (refers to imagenameslist)
% 		    intensityMax   : highest value image can have (calc from bitdepth)
% 					
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Apr 04          Cleaned up source, make function independent of GUI handle

imageDirectory = ptJob.imagedirectory;
firstImageNr = ptJob.firstimage;
lastImageNr = ptJob.lastimage;
imageNameList = ptJob.imagenameslist;
intensityMax = ptJob.intensityMax;

% This is the 7x7 area (in pixels) where we will average the thresholds over
intensityArea = 49;

% Goto the image directory and fetch the first image (subtract the background and
% do normalization as well)
cd (imageDirectory);
name = char (imageNameList (firstImageNr));
firstTempImage = imreadnd2 (name, 0, intensityMax);
[firstImage, backgroundFirst] = ptGetProcessedImage (firstTempImage, 20);
clear firstTempImage;

% And fetch the last image as well
name = char (imageNameList (lastImageNr));
lastTempImage = imreadnd2 (name, 0, intensityMax);
[lastImage, backgroundLast] = ptGetProcessedImage (lastTempImage, 20);
clear lastTempImage;

% Get the image size
[img_h,img_w] = size (firstImage);

%--------------------------------------------------------------------------------

% Get the intensity values of the background of the first image
firstBackgroundImage = figure, imshow (firstImage, []);
title('Click on the background (at least 8 times) on an evenly spread out area. Right-click or press enter to finish.');

% Get the coordinates the user clicks on with the mouse
[X,Y] = getpts (firstBackgroundImage);

intensity = [];
for dots = 1 : size(X,1)
   if round(X(dots)) > 3 & round(X(dots)) < img_w - 3 & ...
      round(Y(dots)) > 3 & round(Y(dots)) < img_h - 3  
      intensity(end+1) = sum (sum (firstImage (round(Y(dots)) - 3 : round(Y(dots)) + 3, round(X(dots)) - 3 : round(X(dots)) + 3))) / intensityArea;
   end
end

backLevelFirst = sum (intensity) / length (intensity);

% Close the figure
close

% And clear the memory space for X and Y
clear X; clear Y;

%--------------------------------------------------------------------------------

% Get the intensity values of the nuclei of the first image
firstNucleusImage = figure, imshow (firstImage, []);
title('Click on the nuclei (at least 8 times) on an evenly spread out area. Right-click or press enter to finish.');

% Get the coordinates the user clicks on with the mouse
[X,Y] = getpts (firstNucleusImage);

intensity = [];
for dots = 1 : size(X,1)
    
	if round(X(dots)) > 3 & round(X(dots)) < img_w - 3 & ...
           round(Y(dots)) > 3 & round(Y(dots)) < img_h - 3  
      
        intensity(end+1) = sum( sum( firstImage(round(Y(dots)) - 3 : round(Y(dots)) + 3, round(X(dots)) - 3 : round(X(dots)) + 3))) / intensityArea;
    end
end
nucLevelFirst = sum(intensity) / length(intensity);

% Close the figure
close

% And clear the memory space for X and Y
clear X; clear Y;

%--------------------------------------------------------------------------------

% Get the intensity values of the halos of the first image
firstHaloImage = figure, imshow(firstImage,[]);
title('Click on the halos (at least 8 times) on an evenly spread out area. Right-click or press enter to finish.');

% Get the coordinates the user clicks on with the mouse
[X,Y] = getpts (firstHaloImage);
intensity = [];
for dots = 1 : size(X,1)
    
	if round(X(dots)) > 3 & round(X(dots)) < img_w - 3 & ...
           round(Y(dots)) > 3 & round(Y(dots)) < img_h - 3  
      
        intensity(end+1) = sum (sum (firstImage (round(Y(dots)) - 3 : round(Y(dots)) + 3, round(X(dots)) - 3 : round(X(dots)) + 3))) / intensityArea;
    end
end
haloLevelFirst = sum(intensity) / length(intensity);

% Close the figure
close

% And clear the memory space for X and Y
clear X; clear Y;

%--------------------------------------------------------------------------------

% Get the intensity values of the background of the last image
lastBackgroundImage = figure, imshow(lastImage,[]);
title('Click on the background (at least 8 times) on an evenly spread out area. Right-click or press enter to finish.');

% Get the coordinates the user clicks on with the mouse
[X,Y] = getpts (lastBackgroundImage);
intensity = [];
for dots = 1 : size(X,1)
    
	if round(X(dots)) > 3 & round(X(dots)) < img_w - 3 & ...
           round(Y(dots)) > 3 & round(Y(dots)) < img_h - 3  
      
        intensity(end+1) = sum (sum (lastImage (round(Y(dots)) - 3 : round(Y(dots)) + 3, round(X(dots)) - 3 : round(X(dots)) + 3))) / intensityArea;
    end
end
backLevelLast = sum(intensity) / length(intensity);

% Close the figure
close

% And clear the memory space for X and Y
clear X; clear Y;

%--------------------------------------------------------------------------------

% Get the intensity values of the nuclei of the last image
lastNucleusImage = figure, imshow(lastImage,[]);
title('Click on the nuclei (at least 8 times) on an evenly spread out area. Right-click or press enter to finish.');

% Get the coordinates the user clicks on with the mouse
[X,Y] = getpts (lastNucleusImage);
intensity = [];
for dots = 1 : size(X,1)
    
	if round(X(dots)) > 3 & round(X(dots)) < img_w - 3 & ...
           round(Y(dots)) > 3 & round(Y(dots)) < img_h - 3  
      
        intensity(end+1) = sum (sum (lastImage (round(Y(dots)) - 3 : round(Y(dots)) + 3, round(X(dots)) - 3 : round(X(dots)) + 3))) / intensityArea;
    end
end
nucLevelLast = sum(intensity) / length(intensity);
% Close the figure
close

% And clear the memory space for X and Y
clear X; clear Y;

%--------------------------------------------------------------------------------

% Get the intensity values of the halos of the last image
lastHaloImage = figure, imshow(lastImage,[]);
title('Click on the halos (at least 8 times) on an evenly spread out area. Right-click or press enter to finish.');

% Get the coordinates the user clicks on with the mouse
[X,Y] = getpts (lastHaloImage);
intensity = [];
for dots = 1 : size(X,1)
    
	if round(X(dots)) > 3 & round(X(dots)) < img_w - 3 & ...
           round(Y(dots)) > 3 & round(Y(dots)) < img_h - 3  
      
        intensity(end+1) = sum (sum (lastImage (round(Y(dots)) - 3 : round(Y(dots)) + 3, round(X(dots)) - 3 : round(X(dots)) + 3) ) ) / intensityArea;
    end
end
haloLevelLast = sum(intensity) / length(intensity);

% Close the figure
close

% And clear the memory space for X and Y
clear X; clear Y;

return;
