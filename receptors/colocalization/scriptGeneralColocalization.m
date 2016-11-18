function [MD] = scriptGeneralColocalization%(test_img)

%% Create Movie Data object
MD = MovieData('/home2/avega/Documents/test_0001.tif'); %Indicate image file to be analyzed; all channels should be in single tiff file

%% Initialize and add all processes
process = SubResolutionProcess(MD); %Detection
MD.addProcess(process);
process = MultiThreshProcess(MD); %Masking
MD.addProcess(process);
process = ColocalizationProcess(MD); %Colocalization
MD.addProcess(process);


%% Detection Process
% The core function being used here is detectSubResFeatures2D_StandAlone
p = MD.getProcess(1).getParameters();
p.ChannelIndex = 2; %Channel in image to undergo detection process 
p.detectionParam.psfSigma = 1; %Point spread function (approximated by Gaussian) standard deviation (in pixels) 
p.detectionParam.alphaLocMax = 0.05;%alpha-value for initial detection of local maxima
p.detectionParam.testAlpha = struct('alphaR',0.1,'alphaA',0.1,'alphaD',0.1,'alphaF',0);%alpha-values for Gaussian mixture-model fitting
p.detectionParam.doMMF = 0;%1 if adding Gaussians iteratively when mixture-model fitting, 0 otherwise
MD.getProcess(1).setParameters(p); %Save parameters
MD.getProcess(1).run; %Run process


%% Masking Process
% The core function used for masking is calcCellBoundaryImage
p = MD.getProcess(2).getParameters();
p.ChannelIndex = 3;%Channel in image to undergo masking process, usually continuum channel
p.GaussFilterSigma = 2;% Sigma of gaussian filter used to smooth image
p.MaxJump = 1; %If function fails to find a threshold in a stack of images, any value >=1 indicates to use the previous threshold 
MD.getProcess(2).setParameters(p);
MD.getProcess(2).run;


%% Colocalization
% The core function used for colocalization analysis is colocalMeasurePt2Cnt
p = MD.getProcess(3).getParameters();
p.ChannelRef = 2; %Punctate channel which underwent detection process
p.ChannelObs = 3; %Continuum channel
p.ChannelMask = 3; %Channel that was masked
p.SearchRadius = 2; %Radius around detection to use to read out intensity
p.RandomRuns = 1;% Number of times randomized data is analyzed
MD.getProcess(3).setParameters(p);
MD.getProcess(3).run;
        
        
        

end