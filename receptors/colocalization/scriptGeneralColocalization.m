function [] = scriptGeneralColocalization()

     %% Create Movie Data objec
        MD = MovieData('dataSample.tif'); %Indicate image file to be analyzed; all channels should be in single tiff file
     %% Initialize and add all processes
        process = SubResolutionProcess(MD); %Detection
        MD.addProcess(process);
        process = MultiThreshProcess(MD); %Masking
        MD.addProcess(process);
        process = ColocalizationProcess(MD); %Colocalization
        MD.addProcess(process);
        
        
    %% Detection Process
        p = MD.getProcess(1).getParameters();
        p.ChannelIndex = 1; %Channel in image to undergo detection process 
        p.detectionParam.psfSigma = 1; %Point spread function sigma of image (in pixels)
        p.detectionParam.testAlpha = struct('alphaR',0.1,'alphaA',0.1,'alphaD',0.1,'alphaF',0);%alpha-values for detection statistical tests
        p.detectionParam.alphaLocMax = 0.05;%alpha-value for initial detection of local maxima
         p.detectionParam.doMMF = 0;%1 if mixture-model fitting, 0 otherwise
        MD.getProcess(1).setParameters(p); %Save parameters
        MD.getProcess(1).run; %Run process
        

    %% Masking Process
        p = MD.getProcess(2).getParameters();
        p.ChannelIndex = 1;%Channel in image to undergo masking process, usually continuum channel
        p.GaussFilterSigma = 2;% Sigma of gaussian filter used to smooth image
        p.MaxJump = 1; %If function fails to find a threshold in a stack of images, any value > 1 indicates to use the previous threshold 
        MD.getProcess(2).setParameters(p);
        MD.getProcess(2).run;

    %% Colocalization
        p = MD.getProcess(3).getParameters();
        p.ChannelRef = 2; %Punctate channel which underwent detection process
        p.ChannelObs = 1; %Continuum channel
        p.ChannelMask = 1; %Channel that was masked
        p.SearchRadius = 2; %Radius around detection used for analysis
        p.RandomRuns = 1;% Number of times randomized data is analyzed
        MD.getProcess(3).setParameters(p);
        MD.getProcess(3).run;

end