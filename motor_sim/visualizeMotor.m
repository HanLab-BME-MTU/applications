function [] = visualizeMotor(param,motorInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualizeMotor plots the positions vesicles with each motor and
% colorcoats each motor according to the state of the motor
%
% This function
%
% INPUT         motorInfo: structure containg vesicle positions,
%                           motor positions, and motor states
%               param: structure containing all
%                       parameters used in the simulation
%
% OUTPUT
%
% DEPENDENCES
%
% REMARKS
%
% Created March 8, 2008 by DA Nunez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imaging parameters
frameHeight = 20;
frameWidth = 100;
pixelSize = 0.126; %position data is in microns; pixelSize is in microns
minOffset = 5;

%make directory for images
HOME = getenv('HOME');
userSpecifiedDirectory = uigetdir([HOME filesep 'Data' filesep 'vesicle_simulation']);
mkdir(userSpecifiedDirectory,'imagesWithMotors');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This loop creates the individual images
for iTime = 1 : param.time_max/param.time_step
    currentFrame = zeros(frameHeight, frameWidth,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This loop puts each vesicle in a given image
    for iVesicle = 1 : param.vesicle_number
        
        figure;
        hold on;
        
        plot (motorInfo.vesiclePosition(1,iVesicle,iTime),frameHeight/2,'bx');
        
        %currentVesiclePosH = ceil(motorInfo.vesiclePosition(1,iVesicle,iTime)/0.126) + minOffset;
        %currentFrame( round(0.5*frameHeight) , max(1, currentVesiclePosH - 1) : min(currentVesiclePosH, frameWidth),:) = 200;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This loop puts each kinesin in a given image. Kinesin will be
        %represented by a triton/pawlike point composed of a cross of pixels
        %with an extra two pixels in the direction of forward (anterograde)
        %motion fo the kinesin motor. By deviding the pixel size by three, one
        %pixel of the image allows for 9 pixels to control shape of kinesin
        %marker.
        for iKinesin = 1 : param.total_kinesin
            %color motor blue if dettached, red if walking, green if
            %paused
            %currentKinesinPosH = ceil(motorInfo.kinesinPosition(iKinesin,iVesicle,iTime)/pixelSize) + minOffset;
            %currentFrame( round(0.5*frameHeight) , max(1, currentKinesinPosH) : min(currentKinesinPosH + 3, frameWidth),motorInfo.kinesinState(iKinesin,iVesicle,iTime)) = 400;
            %currentFrame( round(0.5*frameHeight)-1:round(0.5*frameHeight)+1, max(1, currentKinesinPosH) : min(currentKinesinPosH, frameWidth),motorInfo.kinesinState(iKinesin,iVesicle,iTime)) = 400;
            
            if motorInfo.kinesinState(iKinesin,iVesicle,iTime) ==1
                color = 'dr-';
            elseif motorInfo.kinesinState(iKinesin,iVesicle,iTime) ==2
                color = 'db-';
            else
                color = 'dg-';
            end
            
            plot (motorInfo.kinesinPosition(iKinesin,iVesicle,iTime),frameHeight/2,color);
        end %of kinesin loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This loop puts each kinesin in a given image. Kinesin will be
        %represented by a triton/pawlike point composed of a cross of pixels
        %with an extra two pixels in the direction of forward (anterograde)
        %motion fo the kinesin motor. By deviding the pixel size by three, one
        %pixel of the image allows for 9 pixels to control shape of kinesin
        %marker.
        for iDynein = 1 : param.total_dynein
            %color motor blue if dettached, red if walking, green if paused
            %currentdyneinPosH = ceil(motorInfo.dyneinPosition(iDynein,iVesicle,iTime)/pixelSize) + minOffset;
            %currentFrame( round(0.5*frameHeight)-2:round(0.5*frameHeight)+2 , max(1, currentdyneinPosH)-2 : min(currentdyneinPosH+2, frameWidth),motorInfo.dyneinState(iDynein,iVesicle,iTime)) = 400;
            
            if motorInfo.dyneinState(iDynein,iVesicle,iTime) ==1
                color = 'sr-';
            elseif motorInfo.dyneinState(iDynein,iVesicle,iTime) ==2
                color = 'sb-';
            else
                color = 'sg-';
            end
            
            plot (motorInfo.dyneinPosition(iDynein,iVesicle,iTime),frameHeight/2,color);
        end %of dynein loop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end %of vesicle for loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This loop names the images
    if (iTime < 10)
        imageName = [param.result_id 'frame_00', num2str(iTime), '.tif'];
    elseif (iTime < 100)
        imageName = [param.result_id 'frame_0', num2str(iTime), '.tif'];
    else
        imageName = [param.result_id 'frame_', num2str(iTime), '.tif'];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %filter image with gaussian
    %filtered_image = Gauss2D(currentFrame, 1);
    
    %save image
    %imwrite([uint8(currentFrame)],[userSpecifiedDirectory filesep 'imagesWithMotors' filesep imageName], 'tiff');
    
    saveas(gcf,[userSpecifiedDirectory filesep 'imagesWithMotors' filesep imageName], 'tiff');
    
    %print progress
    fprintf('Making Movie Progress:  %f\n', iTime/param.time_max);
    
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
