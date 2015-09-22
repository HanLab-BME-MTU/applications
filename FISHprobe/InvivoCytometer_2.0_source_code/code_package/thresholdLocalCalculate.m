function [level_img, level_whole] = thresholdLocalCalculate(imageIn,choice_of_threshold,level_local_radius, pace, varargin)
% local thresholding based on global thresholding level using one of the
% three methods "Otsu","Rosin",or "FluorescenceImage"
%
% level_img = thresholdOtsu_local(imageIn,level_local_radius, pace, showPlots)
% 
% This function selects a threshold for the input fluorescence image using
% the thresholding method the user askes for as in choice_of_threshold,
% then find the local thresholds also
% decided by thresholding method the user askes to get a local threshold
% map. After some smoothing on this threshold map(not included here) can used to segment the input
% image.
% 
% Input:
% 
%   imageIn:            2D input image to be thresholded.
%   choice_of_threshold: string input "Otsu","Rosin",or "FluorescenceImage"
%   level_local_radius: the radius of local patch
%   pace:               the pace to calculate local threshold, mostly to
%                       speed up the process which is computational expensive
%
% Output:
% 
%   level_img - The intensity value image selected for local thresholding,
%               before further smoothing over it
%
% Liya Ding, 1/2012, cleaned 12/2012

ip=inputParser;
ip.addRequired('imageIn',@isnumeric);
ip.addRequired('choice_of_threshold', ...
               @(x) ((ischar(x) && ismember(x, {'Otsu', 'Rosin', 'FluorescenceImage'})) || ...
                      isa(x, 'function_handle')));
ip.addRequired('level_local_radius',@isnumeric);
ip.addRequired('pace',@isnumeric);
ip.addOptional('showPlots',0,@isnumeric)
ip.parse(imageIn, choice_of_threshold, level_local_radius, pace, varargin{:});
showPlots=ip.Results.showPlots;

% Define the half pace 
half_pace = round(pace-1)/2;

% Convert to double if necessary
imageIn = double(imageIn);

% Find nonzero values (if due to masking)
nzInd = find(imageIn);

if isempty(nzInd) % all pixels have zero intensity (have been masked)
    level_img = zeros(size(imageIn)) + eps;
    level_whole = eps;
    return;
end

% Get minumum and maximum pixel values in image
minSignal = min(imageIn(nzInd));
maxSignal = max(imageIn(nzInd));

% Normalize nonzero value between 0 and 1
imageInNorm = zeros(size(imageIn));
imageInNorm(nzInd) = (imageIn(nzInd)- minSignal) / (maxSignal - minSignal);

try
    % Obtain the global threshold using the method user askes for 
    if isa(choice_of_threshold, 'function_handle')
        level = choice_of_threshold(imageInNorm);                
    else        
        switch choice_of_threshold
            case 'Otsu'
                level = thresholdOtsu(imageInNorm);
            case 'Rosin'
                level = thresholdRosin(imageInNorm);
            case 'FluorescenceImage'
                level = thresholdFluorescenceImage(imageInNorm);
        end
    end    
    level_whole = level*(maxSignal - minSignal)+minSignal;
catch
    level = nan;
    level_whole = nan;
end

% Initialized the threshold map(level_img) as setting everywhere as the global level
level_img = imageIn*0+level_whole;

% If the radius defined here is larger than the size of image, simply
% return with the global threshold
if level_local_radius> min(size(imageIn))
    return;
end

x_grid = level_local_radius + 1 : pace : size(level_img,1) - level_local_radius;
y_grid = level_local_radius + 1 : pace : size(level_img,2) - level_local_radius;


level_local_radius = round(level_local_radius);
half_pace = round(half_pace);
pace=round(pace);

% For every patch according to the pace(grid)
for img_x = level_local_radius + 1 : pace : size(level_img,1) - level_local_radius
    for img_y = level_local_radius + 1 : pace : size(level_img,2) - level_local_radius
% 
%         img_x = round(img_x);
%         img_y = round(img_y);
        
        % Cropping a local image path
        local_img = imageIn(img_x - level_local_radius : img_x + level_local_radius ,img_y - level_local_radius:img_y + level_local_radius );
        
        % Find nonzero values (due to masking)
        nzInd = find(local_img);
        
        % Get minumum and maximum pixel values in image
        minSignal = min(local_img(nzInd));
        maxSignal = max(local_img(nzInd));
        
        % Normalize nonzero value between 0 and 1
        imageInNorm = zeros(size(local_img));
        imageInNorm(nzInd) = (local_img(nzInd)- minSignal) / (maxSignal - minSignal);
        
        % Find the threshold for this patch
               
        try
            
            if isa(choice_of_threshold, 'function_handle')
                level = choice_of_threshold(imageInNorm);                
            else                    
                switch choice_of_threshold
                    case 'Otsu'
                        level = thresholdOtsu(imageInNorm);
                    case 'Rosin'
                        level = thresholdRosin(imageInNorm);
                    case 'FluorescenceImage'
                        level = thresholdFluorescenceImage(imageInNorm);
                end                
            end
            
        catch
            level = [];
        end
        
        if isnan(level)
            level = [];
        end 
        
        % level could be empty if input patch is all background zeros
        if isempty(level)
            level_img(img_x-half_pace:img_x+half_pace,img_y-half_pace:img_y+half_pace) = level_whole;
            % For the boundary
            if img_x == x_grid(1)
                level_img(1:img_x+half_pace,img_y-half_pace:img_y+half_pace) = level_whole;
            end
            if img_x == x_grid(end)
                level_img(img_x-half_pace:end,img_y-half_pace:img_y+half_pace) = level_whole;
            end
            if img_y == y_grid(1)
                level_img(img_x-half_pace:img_x+half_pace,1:img_y+half_pace) = level_whole;
            end
            if img_y == y_grid(end)
                level_img(img_x-half_pace:img_x+half_pace,img_y-half_pace:end) = level_whole;
            end
        else
            
             level = level*(maxSignal - minSignal)+minSignal;       
            
            level_img(img_x-half_pace:img_x+half_pace,img_y-half_pace:img_y+half_pace) = level;
            if img_x == x_grid(1)
                level_img(1:img_x+half_pace,img_y-half_pace:img_y+half_pace) = level;
            end
            if img_x == x_grid(end)
                level_img(img_x-half_pace:end,img_y-half_pace:img_y+half_pace) = level;
            end
            if img_y == y_grid(1)
                level_img(img_x-half_pace:img_x+half_pace,1:img_y+half_pace) = level;
            end
            if img_y == y_grid(end)
                level_img(img_x-half_pace:img_x+half_pace,img_y-half_pace:end) = level;
            end
        end
    end
end


% For the corners

CL = level_local_radius + pace;

level_img(1:CL, 1:CL)= ...
    mean2(level_img(CL+1:CL+10, 1:CL))/2 +...
    mean2(level_img(1:CL,CL+1:CL+10))/2;

level_img(1:CL, end-CL:end)= ...
    mean2(level_img(CL+1:CL+10, end-CL:end))/2 +...
    mean2(level_img(1:CL,end-CL-10:end-CL-1))/2;
 
level_img(end-CL:end, 1:CL)= ...
    mean2(level_img(end-CL-10:end-CL-1, 1:CL))/2 +...
    mean2(level_img(end-CL:end,CL+1:CL+10))/2;

level_img(end-CL:end, end-CL:end)= ...
     mean2(level_img(end-CL-10:end-CL-1,end-CL:end))/2 +...
     mean2(level_img(end-CL:end,end-CL-10:end-CL))/2;
