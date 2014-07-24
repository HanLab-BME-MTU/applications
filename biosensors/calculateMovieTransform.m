function xForm = calculateMovieTransform(baseImage,inImage,varargin)
%CALCULATEMOVIETRANSFORM calculates an alignment transform using reference images
% 
% xForm = calculateMovieTransform
% xForm = calculateMovieTransform(baseImage,inputImage)
% xForm = calculateMovieTransform(baseImage,inputImage,'OptionName1',optionValue1,...) 
% 
% This function can be used to create a transform which can then be used to
% align images taken from two different channels/cameras. This is
% accomplished using reference images of the same object(s) imaged in both
% channels. Usually these are images of multi-spectral fluorescent beads
% taken in both the channels, or a grid-mircometer. Images of biological
% objects can also be used to generate a transform, but this is not
% recommended as there may be other differences between the images in the
% two channels besides their alignment.
% 
% The transformation is generated in two steps. Both steps are optional,
% but at least one step must be enabled, or there's nothing to do! The
% first step specifies an initial transform in one of three ways:
% 
%       1 - For bead images with sparse beads, by detecting beads in both
%       channels and determining their misalignment.
% 
%       2 - For dense bead images, or other images, by asking the user to
%       click on the same object in both images several times.
% 
%       3 - An "initial guess" transform can be directly input. See options
%       below.
% 
% The second step attempts to refine the initial transformation by
% minimizing the difference (RMSD) between the two images via non-linear
% optimization. This step is computationally intensive and may take several
% minutes to complete - please be patient. Additionally, it is not
% guaranteed to work and, depending on the quality of the initial
% transformation and the images, this step may fail in some cases. If this
% step fails for you, try a better initial transform, different images or
% pre-processing of the images. As a last result, you may try using masks
% generated from the images rather than the images themselved.
%  
%
% Input:
% 
%   baseImage - The "reference" image. This is the image that the second
%   channel will be alligned to. Optional. If no input, the user will be
%   asked to click on the image .tif file.
% 
%   inputImage - The image which will be aligned with the baseImage.
%   That is, when the resulting transform is applied to this image it will
%   align with the bsdr image. Optional. If no input, the user will be
%   asked to click on the image file.
% 
%       Note: In general, it doesn't really matter which image is the
%       base and which is the input image. What's important is
%       to remember which is which and apply the transform appropriately!
% 
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('InitialTransform'->image transform). If input, the first step
%       will be skipped and this transform will be used as the intial
%       guess. 
%   
%       ('DoInitial'->true/false) If true, the inital transform will be
%       generated. If false, the refinement will be attempted directly from
%       the input images. This requires that the mis-alignment between the
%       images be fairly small (<5 pixels)
%
%       ('DoRefinment'->true/false) If true, the function will attempt to
%       refine the initial transform. Optional, default is True. This step
%       is highly recommended, only skip it if it is failing!
% 
%       ('TransformType->'projective' or 'polynomial') The
%       type of transform to use to align the images. Default is
%       projective, which should handle most dual-camera image alignment
%       problems. Polynomial allows higher-order distortion of the image,
%       which can attempt to correct for other effects such as chromatic
%       aberration, but is more difficult to refine.
%
%       ('SaveName'->Character array) If input, this specifies a path and
%       file name to save the resulting transform to. If not specified, the
%       user will be asked where to save the transform, but you aren't
%       forced to.
%
%
% Output:
% 
%   xForm - The resulting alignment transform. This transform, when applied
%   to the input image, will align the two images. This transform is in the
%   format used by the matlab image processing toolbox function
%   imtransform.m This can then be used with the function imtransform or
%   transformMovie.m to align images.
%
% 
% Hunter Elliott
% 11/2010
%


%% ------------- Input -------------- %%

if nargin < 1 || isempty(baseImage)
    
    %Ask user to specify an image file
    [fName,fPath] = uigetfile({'*.tif';'*.TIF';'*.tiff';'*.TIFF'},...
                               'Select the base image file:');
    if fName == 0
        error('You must specify a base image file to continue!')
    end
    %Attempt to load it
    try
        baseImage = imread([fPath filesep fName]);
    catch em
        error(['Could not load base image! Error: ' em.message])
    end
    
end
    
if nargin < 2 || isempty(inImage)
    
    %Ask user to specify an image file
    [fName,fPath] = uigetfile({'*.tif';'*.TIF';'*.tiff';'*.TIFF'},...
                               'Select the input image file:');
    if fName == 0
        error('You must specify an input image file to continue!')
    end
    %Attempt to load it
    try
        inImage = imread([fPath filesep fName]);
    catch em
        error(['Could not load base image! Error: ' em.message])
    end
    
end

if ~isequal(size(baseImage),size(inImage))
    error('The images to be aligned must be the same size!')
end

[initXform,doInit,doRefine,xFormType,saveName] = parseInput(varargin);

% ----- Defaults ------ %

if ~isempty(initXform) 
    if ~istransform(initXform)
        error('The specified initial transform is not in the format used by imtransform.m! Check transform!')
    else        
        doInit = false;
    end      
end

if isempty(doInit)
    doInit = true;
end

if isempty(doRefine)
    doRefine = true;
end

if isempty(xFormType)
    xFormType = 'projective';
end%TEMP - CHECK THAT ITS A VALID TYPE!

if ~doInit && ~doRefine
    error('You must enable either the initial or refined transformation or both!')
end

%% ------------- Initial Transformation ----------- %%
    

if doInit
       
    
    %Ask the user how they want to specify the initial transformation
    bPressed = questdlg('How would you like to create the initial transformation?',...
                'Initial Transformation Creation',...
                'Manually','Spot Detection','Cancel Initial xForm','Manually');

            
    disp('Determining intial transform...')
    
    switch bPressed
        
        
        %Let the user manually click on both pictures to produce initial
        %transformation
        case 'Manually'
            
            %Let the user know whats going on
            waitHan = msgbox('After you click "Ok", you will be shown both of the images you want to align. You must click on several points (Called control points) in both images, so that the two images can be aligned based on these points. The more points you click, the better the resulting transform will be! The minimum # of pairs of points for projective transforms is 4 and for polynomial it is 10. Try to spread the points out evenly over the image area, as this will also improve the resulting transformation. Simply close the control-point selection window when you are finished to continue generating the transform.');
            uiwait(waitHan);            
            
            %Call this but scale the images first because it displays to
            %the whole range            
            [cpIn,cpBase]= cpselect(mat2gray(inImage),mat2gray(baseImage),'Wait',true);
               
            if ~isempty(cpIn) && ~isempty(cpBase)
                if strcmp(xFormType,'polynomial')
                    %Use second-order polynomial
                    initXform = cp2tform(cpIn,cpBase,xFormType,2);
                else
                    initXform = cp2tform(cpIn,cpBase,xFormType);
                end
            end
        case 'Spot Detection'
                    
            beadRad = str2double(inputdlg('What is the approximate bead diameter, in pixels?'));
            if isnan(beadRad) || beadRad < 1
                error('Invalid bead radius! Need bead radius to perform detection!');
            end
            
            try
                %Call the bead-alignment routine
                initXform = getTransformFromBeadImages(baseImage,inImage,xFormType,beadRad,1);
            catch em
                error(['The bead-detection based alignment failed! Try manual-alignment or different images! Error : ' em.message])
            end
            
            
        otherwise
            disp('No initial transform used!')
            if ~doRefine
                error('If refinement is disabled, an initial transform MUST be generated!');
            end                    
    
    end
    
    
    %Show the pre- and post-initial transform alignment, if one was created
    fsFigure(.75);
    if ~isempty(initXform)
        subplot(1,2,1)
    end
    image(cat(3,mat2gray(baseImage),mat2gray(inImage),zeros(size(baseImage))));
    hold on,axis image,axis off
    title('Overlay, before any transformation. Red: base, Green: Input')    
    if ~isempty(initXform)
        subplot(1,2,2)        
        xIn = imtransform(inImage,initXform,'XData',[1 size(baseImage,2)],...
                                            'YData',[1 size(baseImage,1)]);
        image(cat(3,mat2gray(baseImage),mat2gray(xIn),zeros(size(baseImage))));
        hold on,axis image,axis off
        title('Overlay, after initial transformation');
    end
    
end


%% ------------- Transformation Refinement ------------ %%

if doRefine
                
    if ~isempty(initXform)
        %If possible, use the initial transformation as initial guess for
        %refinement
        switch xFormType
            
            case 'projective'
                xForm = findOptimalXform(baseImage,inImage,0,xFormType,initXform.tdata.T);    
                
            case 'polynomial'
                xForm = findOptimalXform(baseImage,inImage,0,xFormType,initXform.tdata);    
                
            otherwise
                error(['"' xFormType '" is not a supported transform type!'])
        end
                
    else
        xForm = findOptimalXform(baseImage,inImage,0,xFormType);            
    end
    
else    
    xForm = initXform;
end

if isempty(xForm)
    error('The transform refinement failed! Try: a different initial guess, different images, different transform type. Sorry!')
end

%Show the pre- and post-final transform alignment, if one was created
fsFigure(.75);    
xIn = imtransform(inImage,xForm,'XData',[1 size(baseImage,2)],...
                                'YData',[1 size(baseImage,1)]);
image(cat(3,mat2gray(baseImage),mat2gray(xIn),zeros(size(baseImage))));
hold on,axis image,axis off
title('Overlay, after refined transformation');




%In case it didn't work, we wait until the end to ask to save
if isempty(saveName)
    
    [saveName,saveDir] = uiputfile('*.mat','Save your transform:');   
    
    if isequal(saveDir,0) || isequal(saveName,0)
        %In case they clicked "cancel"
        saveName = [];
    else
        saveName = [saveDir saveName];
    end
    
end
%Check again in case they cancelled the dialogue
if ~isempty(saveName)
    save(saveName,'xForm');
end

function [initXform,doInit,doRefine,xFormType,saveName] = parseInput(argArray)



%Init output
initXform = [];
doInit = [];
doRefine = [];
xFormType = [];
saveName = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg

    switch argArray{i}                     

        case 'InitialTransform'
            
           initXform = argArray{i+1};

        case 'DoRefinement'
           doRefine = argArray{i+1};
           
        case 'DoInitial'
           doInit = argArray{i+1};

        case 'TransformType'
            xFormType = argArray{i+1};        
    
        case 'SaveName'

           saveName = argArray{i+1};
       otherwise

           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
    end    
end



