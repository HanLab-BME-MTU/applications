function [objectSizes]=determineDetectedObjectSizes(data);
% determine size of clathrin pits in all images for a given condition
%
% SYNOPSIS:
% [pitResults]=determineAllPitsizes(Data);
%
% INPUT     :   Data   = field containign all data for the condition
%
% OUTPUT    :   pitResults = field containing for every movie
%                   .bg = [mean std] - background
%                   .amp = [mean std] - amplitude
%                   .sigma = [mean std] - sigma auf Gaussian fit
%                         
%
% REMARKS   :   
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% DATE: 05-Oct-2007
% last modified
%
%

% loop over all movies

for i=1:length(data)
    
    % locate detection data
    orDirectory = cd;
    currDirectory = data(i).source;
    cd(currDirectory);
    
    readDetectionFilename = 'detection.mat';
    
    
    % load the first image of the movie stack
    [IFileName,IPathName] = uigetfile('*.tif','Select first original image');
    cd(IPathName);
    currentImage = imread(IFileName);
    cd(currDirectory);
    ImageData(i).image = currentImage;
    
    % now load the detection data if necessary, using the appropriate
    % name
        
    cd('DetectionStructures');
    loadfile = load(readDetectionFilename);
    if isfield(loadfile,'detection')
        movieInfo = loadfile.detection;
    elseif isfield(loadfile,'cdet')
        movieInfo = loadfile.cdet;
    else
        error('no detection data file of specified format found');
    end
    ImageData(i).detection = movieInfo;
    
end

fprintf('first image in movie ');

for i=1:length(data)
    
    fprintf('#%02d',i);
    
    currentImage    = ImageData(i).image;
    movieInfo       = ImageData(i).detection;
    
    % current image is currentImage, current pit positions are movieInfo,
    % which has to be converted into the x,y format
    xcoord = movieInfo(1).xCoord(:,1);
    ycoord = movieInfo(1).yCoord(:,1);
    xypos = [xcoord ycoord];
        
    % determine pit sizes in this image
    [CurrMovResults] = determinePitsizeFromImage_2d(currentImage, xypos);
    
    objectSizes(i).bg    = CurrMovResults(:,1);
    objectSizes(i).amp   = CurrMovResults(:,2);
    objectSizes(i).sig   = CurrMovResults(:,3);
    
    cd(orDirectory);
    
    fprintf('\b\b\b');
    
end % of loop

allSig = vertcat(objectSizes(:).sig);
xvec = [0:0.02:3];
bar(xvec+0.01,histc(allSig,xvec));

end % of function




%% =======================================================================
%                           Main Subfunction
%=========================================================================

function [results]=determinePitsizeFromImage_2d(image, xypos);
% determine size of clathrin pits in image (in pixels)
%
% SYNOPSIS:
% [results]=determinePitsizeFromImage(image, xypos);
%
% INPUT     :   image   = tiff image containing pits
%               xypos   = nx2 vector containing x,y positions of pits,
%                         which can be extracted from detection results
%
% OUTPUT    :   results = vector containing the results of the Gaussian fit
%                         to the pits, with
%                         [Background Amplitude Gausswidth]; 
%
%
% REMARKS   :   Note that for the Gaussian fit, the x,y positions are
%               considered known and kept fixed in this version of the
%               function. Thus, lack of precise position determination
%               will result in a widening of the Gauss distribution and
%               larger values of sigma, so that for unprecise positions,
%               this function needs to be modified to include a fit to x,y
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% DATE: 05-Oct-2007
% last modified
%
%
[isx,isy] = size(image');
radius = 5;


% positions of mini-window in grid form
[miniImX, miniImY] = ndgrid(-radius:radius,-radius:radius);
miniDist = sqrt(miniImX.^2 + miniImY.^2);
% extract from the square grid positions those below the specified radius
% (convert to circular area)
[XinPix, YinPix] = find(miniDist<=radius);
XinPix = XinPix-(radius+1);
YinPix = YinPix-(radius+1);
% now XinPix and YinPix contain the positions of the pixels in the circular
% area around the origin with radius radius


% f1 = figure;
% subplot(1,2,1);
% imshow(image,[]);
% hold on

IntMin = min(image(:));
IntMax = max(image(:));

fprintf(' progress');

% loop over know positions
for i=1:length(xypos)
    
    pr = round(100*i/length(xypos));
    fprintf(' %03d',pr);
    
    % cut out appropriate region in image
    cpos = xypos(i,:);
    x0 = cpos(1);
    y0 = cpos(2);
            
    x0i = round(x0);
    y0i = round(y0);
    
    % crop limits x
    x1 = max( x0i-radius,1 );   
    x2 = min( x0i+radius,isx ); 
    % distance of objects from x-edge
    xl1 = x0i - x1 + 1;
    xl2 = isx - x0i + 1;
    
    % crop limits x
    y1 = max( y0i-radius,1 );   
    y2 = min( y0i+radius,isy ); 
    % distance of objects from y-edge
    yl1 = y0i - y1 + 1;
    yl2 = isy - y0i + 1;
    
    cutOutIm = image(y1:y2,x1:x2);
    
        
    if min([xl1,xl2,yl1,yl2])>4
        
        % for the fitting, yl and xl have to be switched just because of the
        % way the convention for x,y is inside the function 
        imax = double(max(cutOutIm(:)));
        imean = double(mean(cutOutIm(:)));
        
        % estimates and startvector have the form
        % [amplitude sigma x-position y-position offset]
        
        % first round: fix sigma and allow positions and amplitudes to
        % converge freely
        startvec1 = [imax-imean 1.5     yl1     xl1     imean];
        fixvec1   = [0          1       0       0       0    ];
        [estimates1] = fit2DGaussToImage(cutOutIm,startvec1,fixvec1);
        
        % second round: fix background and positions and now only allow the
        % height above background and sigma to converge
        startvec2 = estimates1;
        fixvec2   = [0          0       1       1       1    ];
        [estimates2] = fit2DGaussToImage(cutOutIm,startvec2,fixvec2);
        
                
        % enter sigma into list
        Gausswidth(i) = estimates2(2);
    
        % amplitude above background
        Amplitude(i) = estimates2(1);
    
        % background
        Background(i) = estimates2(5);
    
    else
        % enter sigma into list
        Gausswidth(i) = nan;
    
        % amplitude above background
        Amplitude(i) = nan;
    
        % background
        Background(i) = nan;
    end
    
    fprintf('\b\b\b\b');
    
    
end % of loop

fprintf('\b\b\b\b\b\b\b\b\b');


% average sigma
results = [Background' Amplitude' Gausswidth']; 

end % of function



%% ========================================================================
%                       Secondary Subfunction
%==========================================================================



function [estimates] = fit2DGaussToImage(Image, ParameterGuess, ParameterFix)
% function fits 2-D-Gaussian function to an image, based on provided
% estimates of position
%
% INPUT:    Image = 2-D image (2-D matlab matrix)
%           ParameterGuess (optional) 
%               = vector containing guess for all parameters; it
%               should be in the form
%               [amplitude  sigma   x-position  y-position  offset]
%               DEFAULT (if no value is entered explicitly) is
%               [imax-imin  1.2     xpos-of-max ypos-of-max imin]
%           ParameterFix (optional)
%               = vector fixing the value of given parameters
%
% OUTPUT:   estimates = vector of fitted estimate of the parameters, has
%                   the form:
%                   [amplitude sigma x-position y-position offset]
% 
% last modified: 
% Dinah Loerke, 07/10/08

%% =======================================================================
%  set fixvector

global fixvec

fixvec = [0 0 0 0 0];
if nargin>2
    fixvec = ParameterFix;
end


%========================================================================
%% STEP 1:
% make x- and y-vector to correspond to image positions
%========================================================================

% figure;
% imshow(Image,[]); 

[xs,ys]=size(Image);

xymat(:,:,1) = 0*double(Image);
xymat(:,:,2) = 0*double(Image);
[xymat(:,:,1), xymat(:,:,2)] = ndgrid(1:xs,1:ys);


%========================================================================
%% STEP 2:
% set startpoint for fitting
%========================================================================

% image mean
imean = mean(double(Image(:)));
% image max
imax = max(double(Image(:)));
% pos of max
[pmaxx,pmaxy] = find(Image==imax);



%===========================
% NOTE: if you want to keep the value of sigma fixed (e.g. if you know the
% size of your PSF, have very noisy data, and want to improve the
% convergence of the other parameters), then add a line here: 
% start_point(2) = fixed_value (in pixels);
%===========================


% if guess is entered explicitly, use the entered guess values
% else enter default values
if nargin>1
    start_point = ParameterGuess;
           
else
    % the convergence of the fit depends on a halfway decent estimate of
    % the background, therefore set the background value as well
    % background = mean
    start_point(5) = imean;
    % amplitude = max-mean
    start_point(1) = double(imax-imean);
    % sigma = 1/5 of shorter side length, at minimum 3 pixels
    start_point(2) = max((min([xs,ys])/6),3);

    % x-position = x-position of max
    start_point(3) = pmaxx;
    % y-position = y-position of max
    start_point(4) = pmaxy;
end

% display the choice of startpoint on the image, to allow the user to
% check whether the values make sense (note the x- and y-axis in this
% form of image display)
%     hold on;
%     plot(start_point(4),start_point(3),'r.'); 
%     xlabel('y-pos');
%     ylabel('x-pos');
%     title('position guess');
%     pause(0.2);
%     figure;


%========================================================================
%% STEP 3:
% perform fit
%========================================================================

CurrOptions = optimset('Display', 'off');
estimates = fminsearch(@Gauss2Dfun, start_point, CurrOptions);
% expfun accepts curve parameters as inputs and outputs sse,
% the sum of squares error for FittedFunction - Data.
    function sse = Gauss2Dfun(params)
        
        % if any parameters are fixed, do it here
        if max(fixvec)>0
            for p=1:5
                if fixvec(p)==1
                    params(p) = start_point(p);
                end
            end
        end
        
        % build up current fitted image and error        
        for n=1:xs
            for m=1:ys
                FittedCurve(n,m) = func2DGauss(xymat(n,m,1),xymat(n,m,2),params);
                %FittedCurveBin(n,m) = func2DGaussBin(xymat(n,m,1),xymat(n,m,2),params);
                ErrorVector(n,m) = FittedCurve(n,m) - double(Image(n,m));
            end
        end

%         % update figure with fit (display as crosshair)
%         x0 = params(3); xlo = params(3)-params(2); xhi = params(3)+params(2);
%         y0 = params(4); ylo = params(4)-params(2); yhi = params(4)+params(2);
%         imshow(Image,[]); hold on; 
%         plot([y0 y0],[xlo xhi],'c-'); 
%         plot([ylo yhi],[x0 x0],'c-'); 
%         hold off; pause(0.01);
        % comment/uncomment up to here
        
        
         % if any parameters are fixed, re-set their values to input here
        if max(fixvec)>0
            for p=1:5
                if fixvec(p)==1
                    params(p) = start_point(p);
                end
            end
        end
        
        sse = sum( ErrorVector(:).^ 2);
    end


% if any parameters are fixed, re-set their values to input here to ensure
% the output equals the input
if max(fixvec)>0
    for p=1:5
        if fixvec(p)==1
            estimates(p) = start_point(p);
        end
    end
end
                

% FinalCurve(n,m) = func2DGauss(xymat(n,m,1),xymat(n,m,2),estimates);
% x0 = estimates(3); xlo = estimates(3)-estimates(2); xhi = estimates(3)+estimates(2);
% y0 = estimates(4); ylo = estimates(4)-estimates(2); yhi = estimates(4)+estimates(2);
% % imshow(Image,[]); hold on; 
% % plot([y0 y0],[xlo xhi],'r-'); 
% % plot([ylo yhi],[x0 x0],'r-'); 
% 


end % of function




%========================================================================
% subfunction
%========================================================================



function [ivalue]=func2DGauss(x,y,parameters)

% amplitude over background
amp = abs(parameters(1));
% sigma (symmetric in x and y!!)
sig = parameters(2);
% location of Gauss in x
x0 = parameters(3);
% location of Gauss in y
y0 = parameters(4);
% background offset
offset = parameters(5);

% zz = square distance of point x,y from center x0,y0
zz = ( (x-x0).^2 + (y-y0).^2 );
ivalue = offset + amp .* exp( -zz/(2*sig^2) );


end % of subfunction


