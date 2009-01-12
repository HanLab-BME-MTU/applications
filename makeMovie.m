function makeMovie(param,tracks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeMovie takes param and tracks as imputs to create sequential images of vesicle position 
% out of the position data from simulations.
%
% INPUT         param             : simulation parameter structure
%               tracks            : simulation output structutre
%
% OUTPUT        Images are written to [HOME filesep 'Data' filesep
%               'vesicle_simulation' filesep param.result_id]
%
% DEPENDENCES   makeMovie makes no use of other functions
%               makeMovie is used by 
%
% REMARKS
%
% Created December 6, 2007 by DA Nunez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imaging parameters
frameHeight = 128;
frameWidth = 256;
pixelSize = 0.126;
minOffset = 5;

%make directory for images
HOME = getenv('HOME');
userSpecifiedDirectory = uigetdir([HOME filesep 'Data' filesep 'vesicle_simulation']);
mkdir(userSpecifiedDirectory,'images');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This loop creates the individual images
for i = 1 : param.time_max/param.sampling_rate
    currentFrame = zeros(frameHeight, frameWidth);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This loop puts each vesicle in a given image
    for j = 1 : param.vesicle_number
    currentMotorPosH = ceil(tracks(j).points(i,1)) + minOffset;
    currentFrame( round(0.5*frameHeight) , max(1, currentMotorPosH - 1) : min(currentMotorPosH + 1, frameWidth)) = 225;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This loop names the images
    if (i < 10)
        imageName = [param.result_id,'frame_00', num2str(i), '.tif'];
    elseif (i < 100)
        imageName = [param.result_id,'frame_0', num2str(i), '.tif'];
    else
        imageName = [param.result_id,'frame_', num2str(i), '.tif'];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %filter image with gaussian
    filtered_image = Gauss2D(currentFrame, 1);
    
    %save image
    imwrite(uint8(filtered_image), [userSpecifiedDirectory filesep 'images' filesep imageName], 'tiff');
    
    %print progress
    fprintf('Making Movie Progress:  %f\n', i/param.time_max/param.sampling_rate);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%