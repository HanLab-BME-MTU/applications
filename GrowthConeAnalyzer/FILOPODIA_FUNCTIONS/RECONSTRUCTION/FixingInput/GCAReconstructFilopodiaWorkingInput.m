function [analInfoC,TSFigs] = GCAReconstructFilopodia(img,analInfoC,varargin)
%% GCAReconstructFilopodia: Adds Filopodia and Branch Structures to a Veil/Stem Estimation
% This function rebuilds and records the filopodia network around a
% veil/stem mask (in the case of the neurite) or any binary cell mask 
% where small scale (often low fidelity) protrusion detections have been 
% removed. 
% 
% STEP I: It applies a small scale steerable ridge filter followed 
%         by NMS (non-maximum suppression) to detect filopodia. 
%         It then reconstructs the filopodia network from this 
% 
% 

%INPUT: 
%       img: (REQUIRED) RxC double array of image to analyze where R is the height 
%       (ny) and C is the width (nx) of the input image
%
%       analInfoC: (REQUIRED) an 1x1 structure array containting the current frame
%       information corresponding to the veil/stem estimation : 
%       a filoInfo field will be added to this structure 
% 
%       protrusionC: (OPTIONAL) output from the protrusion process (see getMovieProtrusion.m)
%                                Structure with fields:
%                                   .normals: a rx2 double of array of unit normal                               
%                                             vectors along edge: where r is the number of 
%                                             coordinates along the veil/stem edge
%                                   (see output: output.xy_normal in prSamProtrusion.m)
%                              
%                                   .smoothedEdge: a rx2 double array of edge coordinates 
%                                                  after spline parameterization: 
%                                                  where r is the number of
%                                                  coordinates along the
%                                                  veil/stem edge
%                                   (see output: output.pixel_tm1_output in prSamProtrusion)
%                                Default : [] , NOTE: if empty the field
%                                filoInfo(xFilo).orientation for all filodpodia attached
%                                to veil will be set to NaN (not
%                                calculated)- there will be a warning if
%                                to the user if this is the case
% 
%  PARAMETERS: 
%  STEP I: Detect Thin Ridge Structures : 
%      input into gcaMultiscaleSteerableDetector.m 
% 
%      steerFiltOrder: (PARAM) Scalar Default = 4 
%                     See gcaMultiscaleSteerableDetector.m 
%
%      steerFiltScales: (PARAM) Rx1 vector where R is the number of scale over 
%                               which the response will be integrated : Default = 1.5 
%                     See gcaMultiscaleSteerableDetector.m 
% 
%  STEP II: Clean Ridge Candidates For Reconstruction
%      multSTDNMSResponse: Scalar  (Default 3) 
%  
%      smallCCSizeThreshold: Scalar (Default 3 Pixels)
% 
%      
%       
%     
% 
%       p : parameter structure with fields
%          .steerFiltOrder  = order of steerable ridge filter (Default = 4)
%          .steerFiltScales = scales of filopodia (Default: 1.5) to be
%           intergrated over
%          .plots.internal =  make trouble shooting plots for internal
%          filodoia (Default True) 
% 
%       iFrame : Quick fix until get figure handles 
%
%      

% OUTPUT: adds a field to analInfo called filoInfo
% filoInfo is a N structure x 1 structure with fields providing information
% regarding each filopodia: where N is the number of filopodia pieces
% reconstructed in iFrame.  Currently 'filopodia' can be 
% IDed into select branch groups. 
% groups via the field groupCount. 
% (My version notes addition : Nested fields were avoided here to facilitate the data
% extraction in later steps (ie each filopodia was given an field ID 
% specifying branch order rather than nesting the structure)

%% Check Input% NOTE TO SELF FIX INPUT 
filterOrder = 4; 
scales = 1.5; 


 
if nargin<3 
  p.steerFiltOrder  = 4; 
  p.steerFiltScales = 1.5; 
  p.embedded.on = true; % lifeAct only
  p.embedded.troubleshoot = true; 
  p.candRidgeHists = true; 
end 

if nargin<4
   iFrame  = 1 ;   
end 


%% STEP I: Detect Thin Ridge Structures 
    
[maxRes, maxTh ,maxNMS ,scaleMap]= gcaMultiscaleSteerableDetector(img,filterOrder,scales); 

% save these as we will need them as input for later 
analnfoC.filterInfo.maxTh = maxTh;
analnfoC.filterInfo.maxRes = maxRes;

analnfoC.scaleMap = scaleMap; 

%% Mask out all NMS 
%perform an initial thresholding first to get out main components 
[maskBack,backMu,backSig] = gcaEstimateBackgroundArea(img); 

%% Threshold out maxNMS values less than 3* first mode gaussian fit from the area 
% in the image falling within the forframe
 forValues = maxNMS.*~maskBack; % take out background response based on fluorescence intensity
 valuesFilter = forValues(forValues~=0); 
 

    % try a second step where you filter out the background response from the 
    % estimation around the area of interest: keep thresholding on a pixel
    % by pixel basis. 
    
    [respNMSMean,respNMSSTD]   = fitGaussianModeToPDF(valuesFilter); 
    cutoffTrueResponse = respNMSMean+3*respNMSSTD; % can make this a variable 
    n1 = hist(valuesFilter,100);
    %% Plot the Candidate Ridge Thresholding Histogram If Desired for Troubleshooting
        if p.candRidgeHists == true % plot the histogram with cut-off overlay so can see what losing 
         
          TSFig1H = figure('visible','off'); 
          TSFigs(1).h  = TSFig1H;
          TSFigs(1).name = 'candidateRidges'; 
           
          hist(valuesFilter,100); 
          hold on ed
          line([cutoffTrueResponse cutoffTrueResponse],[0,max(n1)],'color','r','Linewidth',2); 
          title('Red line 3*std of first mode'); 
        end
        
    canRidges = maxNMS.*~maskBack;
    canRidges(canRidges<cutoffTrueResponse) = 0; 
    
    analnfoC.filterInfo.ThreshNMS = canRidges; 
       
%% clean the thresholded NMS response using the body estimation information 
    
skelIn = canRidges; 
skelIn = bwmorph(skelIn,'thin',inf ); 



% Main Function (should rename) that performs the reconstructions
[reconstruct,filoInfo,TSFigs] = gcaAttachFilopodiaStructuresMain(img,skelIn,analInfoC,protrusionC);


analInfoC.filoInfo = filoInfo; % 
analInfoC.reconstructInfo = reconstruct;


end % 


        
        
        




