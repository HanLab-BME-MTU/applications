function [analInfoC,TSFigsFinal] = GCAReconstructFilopodia(img,analInfoC,varargin)
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
%
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

%% Input Parser 
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true; 

ip.addRequired('img'); 
ip.addRequired('analInfoC'); 


% PARAMETERS
 ip.addParameter('steerFiltOrder',2,@(x) ismember(x,[2,4])); 
 ip.addParameter('steerFiltScales',1.5); 
 
 
 ip.addParameter('multSTDNMSResponse',3);
 ip.addParameter('smallCCSizeThreshold',3); 
 %ip.addParmaeter(
 
 ip.addParameter('detectEmbedded',true);
 
 ip.addParameter('TSOverlays',true); 
 
ip.parse(varargin{:});
p = ip.Results;  

%% Initiate 
countFigs = 1; 


%% STEP I: Detect Thin Ridge Structures 
    
[maxRes, maxTh ,maxNMS ,scaleMap]= gcaMultiscaleSteerableDetector(img,filterOrder,scales); 

% save these as we will need them as input for later 
analnfoC.filterInfo.maxTh = maxTh;
analnfoC.filterInfo.maxRes = maxRes;

analnfoC.scaleMap = scaleMap; 

%% Mask out all NMS 
%perform an initial thresholding first to get out main components 
[maskBack,~,backSig] = gcaEstimateBackgroundArea(img); 

%% Threshold out maxNMS values less than 3* first mode gaussian fit from the area 
% in the image falling within the forframe
 forValues = maxNMS.*~maskBack; % take out background response based on fluorescence intensity
 valuesFilter = forValues(forValues~=0); 
 

    % try a second step where you filter out the background response from the 
    % estimation around the area of interest: keep thresholding on a pixel
    % by pixel basis. 
    
    [respNMSMean,respNMSSTD]   = fitGaussianModeToPDF(valuesFilter); 
    cutoffTrueResponse = respNMSMean+ip.Results.multSTDNMSResponse*respNMSSTD; % can make this a variable 
    n1 = hist(valuesFilter,100);
%% OPTIONAL TS PLOT 
        if ip.Results.TSOverlays == true % plot the histogram with cut-off overlay so can see what losing 
         
          TSFigs(countFigs).h = figure('visible','off'); 
       
          TSFigs(countFigs).name = 'candidateRidges'; 
           
          hist(valuesFilter,100); 
          hold on 
          line([cutoffTrueResponse cutoffTrueResponse],[0,max(n1)],'color','r','Linewidth',2); 
          title('Red line 3*std of first mode'); 
          countFigs = countFigs+1; % close figure 
        end
 %%      
    canRidges = maxNMS.*~maskBack;
    canRidges(canRidges<cutoffTrueResponse) = 0; 
    
    analInfoC.filterInfo.ThreshNMS = canRidges; 
        
skelIn = canRidges; 
skelIn = bwmorph(skelIn,'thin',inf ); 



% Main Function (should rename) that performs the reconstructions
[reconstruct,filoInfo,TSFigs2] = gcaAttachFilopodiaStructuresMainFixInput(img,skelIn,analInfoC,protrusionC,p);

TSFigsFinal = [TSFigs; TSFigs2]; 

analInfoC.filoInfo = filoInfo; % 
analInfoC.reconstructInfo = reconstruct;


end % 


        
        
        




