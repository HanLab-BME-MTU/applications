%% Track AP2
% Clathrin Detection & Tracking
% If you are using a mixture model, you need to add gsl in the same
% terminal that you intend to open matlab in, prior to opening matlab.
%
% Good examples on MD/packages/etc:
% 
% Assumes indexLSFMData has been performed on the entire directory, and the
% script then goes through and sequentially performs the same task on all
% of the individual MovieData structures.

clc;
clear;

load('/project/bioinformatics/Danuser_lab/3Dmorphogenesis/raw/kdean/IMCD-wDrug/alpha-eGFP-adaptin/170713/MovieList.mat');

for MDIdx = 8%1:1:length(ML.movieDataFile_)
    %MD=MovieData.load('/project/bioinformatics/Danuser_lab/3Dmorphogenesis/raw/kdean/IMCD-woDrug/alpha-eGFP-adaptin/170725/Cell4/analysis/cropped/analysis/movieData.mat')
    
    %% Movie Data
    load(ML.movieDataFile_{MDIdx})
    %MD.outputDirectory_
    MD.reset();
    
    % Add Package to MD
    Package_ = UTrackPackage3D(MD);
    MD.addPackage(Package_);
    MD.getPackage(1).createDefaultProcess(1); % Load default parameters

    
    %% Specify Detection Settings.
    funParams = MD.getPackage(1).processes_{1}.funParams_;
    funParams.alpha = 0.05;
    funParams.fitMixtures = true;
    funParams.MaxMixtures = 5;
    funParams.filterSigmaXY = [1.5];
    funParams.filterSigmaZ = [1.5];
    funParams.algorithmType = {'pointSourceAutoSigmaMixture'};
    detectionMethod = funParams.algorithmType{1}
    % Options include: pointSourceAutoSigmaLM, pointSourceAutoSigma, pointSourceAutoSigmaFit, pointSourceAutoSigmaMixture
    
    % Commit new prameters
    MD.getPackage(1).getProcess(1).setPara(funParams);
    
    %% Run Detection.
    MD.getPackage(1).processes_{1}.run();
    
    %% Specify Tracking Settings.
    MD.getPackage(1).createDefaultProcess(2);
    funParams = MD.getPackage(1).processes_{2}.funParams_
    
    funParams.gapCloseParam.mergeSplit = 1;
    MD.getPackage(1).getProcess(2).setPara(funParams);
    
    %% Run Tracking
    MD.getPackage(1).processes_{2}.run();
    
    %% Create MIPs
    MD.addProcess(ComputeMIPProcess(MD));
    MD.processes_{3}.run;
    
end