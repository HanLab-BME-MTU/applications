%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Create MovieData
	saveFolder = 'analysis';
	ch1 = Channel('mDiaActin_mini/images_mDia1');
	ch2 = Channel('mDiaActin_mini/images_actin');
	% Constructor needs an array of channels and an output directory (for analysis)
	MD = MovieData([ch1 ch2], saveFolder);
	MD.setPath(saveFolder);
	MD.setFilename('movieData.mat');
	% Check image size and number of frames are consistent.
	% Save the movie if successfull
	MD.sanityCheck; % 
	% Set some additional movie properties
	MD.numAperture_=1.4;
	MD.pixelSize_=215; % in nm after binning
	MD.timeInterval_=5;% in sec
	MD.camBitdepth_=16;
    MD.notes_='Created for test purposes';
	% Save the movieData
	MD.save;
	
    %% Run Threshold
	thresProc = ThresholdProcess(MD);
	MD.addProcess(thresProc);
	% Create a segmentation package
	segPackage = SegmentationPackage(MD);
	MD.addPackage(segPackage);
	% % Associate the threshold process to the package
	MD.packages_{1}.setProcess(1,thresProc);
	% Get the threshold parameters so you can modify them
	params = MD.processes_{1}.funParams_;
	% Set the Thresholding Parameters
	params.MethodIndx = 1;
	params.GaussFilterSigma = 1;
	% Resave the parameters
	parseProcessParams(MD.processes_{end}, params);
	% Run the process
	MD.processes_{1}.run(); %
	MD.save

	%% Refine the segmentation mask
    MD.addProcess(MaskRefinementProcess(MD))
    MD.packages_{1}.setProcess(2,MaskRefinementProcess(MD));   
    MD.processes_{2}.run()
    
    %% Run the protrusion vectors
	MD.addProcess(ProtrusionProcess(MD));
	% Get the Protrusion Parameters so you can modify them
	protParams = MD.processes_{end}.funParams_;
	% Set the Protrusion Parameters
	protParams.ChannelIndex = 1; % 
	protParams.SegProcessIndex = 2; % (refined mask)
	% Resave the parameters (if reconfigured)
	parseProcessParams(MD.processes_{end},protParams);
	% Run the process
	MD.processes_{end}.run();
	MD.save
	
    %% Run the Windowing
	MD.addProcess(WindowingProcess(MD));
% 	% Get the Windowing Parameters so you can modify them
% 	windParams = MD.processes_{end}.funParams_;
% 	% Set the Windowing Parameters
% 	windParams.ChannelIdx = 1; %
% 	windParams.SegProcessIndex = 1;
% 	windParams.PerpSize = 1;
% 	windParams.MinSize = 500; % min size of the object to window
% 	windParams.StartContour = 1;
% 	windParams.StartPointPropag = false; % you might want this on for more normal cells I wanted it to remain constant
% 	windParams.ParaSize = 5;
% 	windParams.MethodName = 'ConstantNumber';
% 	windParams.ReInit = 61;
% 
% 	outName = [windParams.MethodName '_windSize_' num2str(windParams.ParaSize) 'ReInit' num2str(windParams.ReInit)] ;
% 	windParams.OutputDirectory = [MD.outputDirectory_ filesep 'windows_' outName];

% 	parseProcessParams(MD.processes_{end},windParams);
	%
	MD.processes_{end}.run();
	MD.save;

	%% Run protrusion sampling process
	MD.addProcess(ProtrusionSamplingProcess(MD));
	% Get the Windowing Parameters so you can modify them
    protParams = MD.processes_{end}.funParams_;
% 	protParams.OutputDirectory = [MD.outputDirectory_ filesep 'protrusion_samples_' outName];
% 	parseProcessParams(MD.processes_{end},protParams);
	MD.processes_{end}.run();

	MD.save;

	disp('Finish Windowing test run script successfully');