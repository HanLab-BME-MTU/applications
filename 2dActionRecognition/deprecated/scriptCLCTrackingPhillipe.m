
loadMD=0;
runDetection=1;
runTracking=1;
runRegistration=0;
runAmiraTranscript=1;

%processFrames=[1:100];
root='/work/gdanuser/proudot/project/EB3-3D-track/data-analysis/four-phases'

if loadMD
end

MDArray=[MD];% MDM11 MDM12 MDP3];

%% loop
for i=1:length(MDArray)
    MD=MDArray(i);
    %% detection
    outputDirDetect=[MD.outputDirectory_ filesep 'CLCTracking'];mkdir(outputDirDetect);
    if(runDetection)
        [movieInfoPSD]=detectEB3(MD,'type','pointSourceAutoSigmaFit','processFrames',processFrames,'showAll',false,'scales',[1.3,2.5],'channel',1,'printAll',false)
        save([outputDirDetect filesep 'detection.mat'],'movieInfoPSD');
        amiraWriteMovieInfo([outputDirDetect filesep 'amiraVertex' filesep 'detect.am'],movieInfoPSD,'scales',[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_]);
               
        % Write local thresholding results. 
        mkdir([outputDirDetect filesep 'mask']);
%         for tidx=1:length(movieInfoPSD)
%             stackWrite(lab{tidx},[outputDirDetect filesep 'mask' filesep 'detect_T_' num2str(tidx,'%05d') '.tif']);
%         end
        
    else
        load([MD.outputDirectory_ filesep 'CLCTracking' filesep 'detection.mat']);
    end
    
    if(runRegistration)
        regProcess=registerTranslation3D(MD,'show',true,'warp',true,'computeShift',true);
        load(regProcess.outFilePaths_{1}); %bad code
        parfor i=1:MD.nFrames_
            jIdx=find((jumpIdx<i));
            for j=jIdx
                movieInfoPSD(i).xCoord(:,1)=movieInfoPSD(i).xCoord(:,1)+displacements{j}.T(4,1);
                movieInfoPSD(i).yCoord(:,1)=movieInfoPSD(i).yCoord(:,1)+displacements{j}.T(4,2);
                movieInfoPSD(i).zCoord(:,1)=movieInfoPSD(i).zCoord(:,1)+displacements{j}.T(4,3);
            end
        end
        save([MD.outputDirectory_ filesep 'pointSourceKDetect' filesep 'detectionReg.mat'],'movieInfoPSD');
    end
    
    %% Tracking
    outputDir=[MD.outputDirectory_ filesep 'pointSourceKDetect' filesep 'tracks'];mkdir(outputDir);
    if runTracking
        [gapCloseParam,costMatrices,kalmanFunctions,...
         probDim,verbose]=tracker_param();
        watch_KF_iter=0;
      
        saveResults.dir =  outputDir; %directory where to save input and output
        saveResults.filename = 'trackResults.mat'; %name of file where input and output are saved

        [tracksFinal,kalmanInfoLink,errFlag] = ...
            trackCloseGapsKalmanSparse(movieInfoPSD, ... 
                                       costMatrices,gapCloseParam,kalmanFunctions,...
                                       probDim,saveResults,verbose);
    end
    
    %%
    if runAmiraTranscript
        tracks=TracksHandle(tracksFinal);
        save([outputDir filesep 'tracksHandle.mat'],'tracks');
        amiraWriteTracks([MD.outputDirectory_ filesep 'pointSourceKDetect' ...
            filesep 'Amira' filesep 'Reg' filesep 'AmiraTracks' filesep 'KinTracks.am'],tracks,'scales',[MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_]);
    end
    
end