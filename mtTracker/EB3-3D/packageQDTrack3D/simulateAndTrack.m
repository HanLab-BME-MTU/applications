method_nb=1;
compute_U_track=1;

outputDir='.'

% density in spots/1 microm^2 (20 px window)
frequencyRange=[1 3.5];
densityRange=[1 2]; % in 1/um^3
percentageGoodLink=cell(size(frequencyRange,2),1);
percentageBadLink=cell(size(frequencyRange,2),1);
percentageMissingLink=cell(size(frequencyRange,2),1);
cpuTime=cell(size(frequencyRange,2),1);

percentageVecIdx=1;

%% TracksSim parametrization and computation
detectionParam.bitDepth = 16; %Camera bit depth
imSize=[10 10 10]; % 2x2x2 microns  ()    
numF=150;
% lftDist=zeros(1,10);
% lftDist(10)=1;

minSpeed=5;  % 1 microns/sec
maxSpeed=5;  % 1 microns/sec

parfor percentageVecIdx=1:length(frequencyRange)
    frequency=frequencyRange(percentageVecIdx);
    percentageGoodLinkTmp=cell(1,length(densityRange));
    percentageBadLinkTmp=cell(1,length(densityRange));
    percentageMissingLinkTmp=cell(1,length(densityRange));
    cpuTimeTmp=cell(1,length(densityRange));
    for densityIdx=1:length(densityRange)
        density=densityRange(densityIdx);
        disp(['Frequency: ' num2str(frequency)]);
        disp(['density: ' num2str(density)]);
        motionParam=struct('diffCoef',[],'confRad',[],'speed1D',[],'probVelSwitch',[],'fracType',[],'probMS',[]);
        
        % Simulation parameterization and execution
        numP=ceil(density*(imSize(1)*imSize(2)*imSize(3))/(5^3));
        motionParam.diffCoef=[0.01 0.01];
        motionParam.confRad=[0.3 0.8];
        motionParam.speed1D=[minSpeed maxSpeed]/frequency;
        motionParam.probVelSwitch=0;
        motionParam.fracType=[0.1 0.5 0.4];
        motionParam.probMS=[];
        
        sigmaLft=5*frequency;
        meanLft=10*frequency;
        x =1:numF;
        lftDist=exp(-(x-meanLft).^2/(2*sigmaLft^2))/sqrt(2*pi*sigmaLft^2);
        intVec=[500 30];
        
        [simMPM,tracksSim] = simulateRandomPlusLinearMotion3D(imSize,numP,lftDist,numF,intVec,motionParam);
        
        %% data translation and amira file
        % detection GT to movie info data structure
        % MovieInfo data structure allows to assign an implicite ID to each spot
        % We have to store this tracks in the ID form to link kalman filter info to real GT tracks property.
        movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),numF,1);
        for i = 1: numF
            xCoord=simMPM(simMPM(:,4*i-3)~=1.,4*i-3);
            movieInfo(i).xCoord=[ xCoord 0.05*ones(size(xCoord))];
            yCoord=simMPM(simMPM(:,4*i-2)~=1.,4*i-2);
            movieInfo(i).yCoord=[ yCoord 0.05*ones(size(yCoord))];
            zCoord=simMPM(simMPM(:,4*i-1)~=1.,4*i-1);
            movieInfo(i).zCoord=[ zCoord 0.05*ones(size(zCoord))];
            amp=simMPM(simMPM(:,4*i)~=0.,4*i)/(2^detectionParam.bitDepth -1);
            movieInfo(i).amp=[ amp 0.0003*ones(size(amp))];
        end
        
        saveFolder=[outputDir filesep 'simulationData/frequencyVariation/frequency_' num2str(frequency) '/'];
        mkdir(saveFolder);
        % transcripting movieInfo ID on tracksSIm
        cumID = cumsum(simMPM(:,1:4:end)~=1,1);
        for track_idx = 1:size(tracksSim,1)
            tracksSim(track_idx).tracksFeatIndxCG=cumID(track_idx,simMPM(track_idx,1:4:end)~=1);
        end
        
        %save([saveFolder filesep 'tracksSim.mat'],'tracksSim');
        
        tracksSimNF=TracksHandle(tracksSim);
        %save([saveFolder filesep 'tracksNewFormatSim.mat'],'tracksSimNF');
        amiraWriteTracks([saveFolder '/amiraSimu/tracks.am'],tracksSimNF,'fillGaps',false)
        
        method_idx=1;
        
        %% U track
        if compute_U_track
            [gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=tracker_param();
            costMatrices(1).parameters.maxSearchRadius=(maxSpeed+1)/frequency;
            costMatrices(1).parameters.kalmanInitParam.searchRadiusFirstIteration = (maxSpeed+1)/frequency;
            schemeName='U_track'; %directory where to save input and output
            
            
            
            %function name
            costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';
            kalmanFunctions.reserveMem  = 'kalmanResMemLM';
            kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
            kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
            kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
            
            %% recording results
            [percentageGoodLinkTmp{densityIdx}, percentageBadLinkTmp{densityIdx}, ...
                percentageMissingLinkTmp{densityIdx},cpuTimeTmp{densityIdx}, ... 
                tracksFinal,tracksNewFormat]= ...
                track_stat_3D(movieInfo,tracksSim,tracksSimNF, ...
                costMatrices,kalmanFunctions,gapCloseParam,saveFolder, ...
                schemeName,frequency);
            
            amiraWriteTracks([saveFolder '/amiraMeasured/tracksMeasured.am'],tracksNewFormat,'fillGaps',false)
        end
        
    end
    percentageGoodLink{percentageVecIdx}=percentageGoodLinkTmp;
    percentageBadLink{percentageVecIdx}=percentageBadLinkTmp;
    percentageMissingLink{percentageVecIdx}=percentageMissingLinkTmp;
    cpuTime{percentageVecIdx}=cpuTimeTmp;
end

%%
save('trackingPerformancesDensitySimu.mat','percentageGoodLink','percentageBadLink','percentageMissingLink','frequencyRange','densityRange')


%% plotting results
figure();

cmapF=winter(length(densityRange));
h=setupFigure(2,1,2,'AxesWidth',3.5,'AxesWidth',7,'DisplayMode', 'print','AspectRatio',1,'YSpace', [2 1.5 1]);
lineStyles={'--','-',':','+-','.-'};
for dIdx=1:length(densityRange)
    goodLink=cellfun(@(x) x{dIdx}(1),percentageGoodLink);
    plot(h(1),frequencyRange,goodLink,lineStyles{dIdx},'color',cmapF(dIdx,:),'LineWidth',2);
    axis(h(1),[min(frequencyRange),max(frequencyRange),0,100])
    title(h(1),'U-track linking percentage')
    legend(h(1),num2str(densityRange'));
    xlabel(h(1),'Acquisition frequency')
    ylabel(h(1),'correct link percentage')
    hold(h(1),'on');
    
    badLink=cellfun(@(x) x{dIdx}(1),percentageBadLink);
    plot(h(2),frequencyRange,badLink,lineStyles{dIdx},'color',cmapF(dIdx,:),'LineWidth',2);
    axis(h(2),[min(frequencyRange),max(frequencyRange),0,170])
    title(h(2),'U-track wrong linking percentage')
    legend(h(2),num2str(densityRange'));
    xlabel(h(2),'Acquisition frequency')
    ylabel(h(2),'wrong link percentage')
    hold(h(2),'on')
end
hold off

%% plotting results function of density
figure();
cmapD=autumn(length(frequencyRange));
lineStyles={'-','--',':','+-','.-'};
for fIdx=1:length(frequencyRange)
    goodLink=cellfun(@(x) x(1),percentageGoodLink{fIdx});
    subplot(2,1,1);
    plot(densityRange,goodLink,lineStyles{fIdx},'color',cmapD(fIdx,:),'LineWidth',2);
    axis([min(densityRange),max(densityRange),0,100])
    title('Correct link')
    legend(num2str(frequencyRange'));
    xlabel('Density')
    ylabel('correct link percentage')
    hold on
    
    subplot(2,1,2);
    badLink=cellfun(@(x) x(1),percentageBadLink{fIdx});
    plot(densityRange,badLink,lineStyles{fIdx},'color',cmapD(fIdx,:),'LineWidth',2);
    axis([min(densityRange),max(densityRange),0,150])
    title('False link')
    xlabel('Density')
    ylabel('wrong link percentage')
    hold on
end
hold off

%% Computing and plotting Jaccard
jaccardCoeff=cell(length(frequencyRange),1);
for fIdx=1:length(frequencyRange)
jaccardCoeff{fIdx}=cellfun(@(g,b,m) g./(g+b+m),percentageGoodLink{fIdx},percentageBadLink{fIdx},percentageMissingLink{fIdx},'unif',0);
end
%%

h=setupFigure(1,1,1,'AxesWidth',3.5,'AxesWidth',3.5,'DisplayMode', 'print','AspectRatio',1,'YSpace', [2 1.5 1]);
cmapF=winter(length(densityRange));
lineStyles={'--','-',':','+-','.-'};
for dIdx=1:length(densityRange)
    goodLink=cellfun(@(x) x{dIdx}(1),jaccardCoeff);
    plot(h,frequencyRange,goodLink,lineStyles{dIdx},'color',cmapF(dIdx,:),'LineWidth',2);
    axis(h,[min(frequencyRange(1:end)),max(frequencyRange(1:end)),0,1])
    title(h,'Dynamic quantification breakpoint')
    legend(h,num2str(densityRange'));
    xlabel(h,'Acquisition frequency')
    ylabel(h,'JSC')
    hold(h,'on');
end   
%%
h=setupFigure(1,1,1,'AxesWidth',3.5,'AxesWidth',3.5,'DisplayMode', 'print','AspectRatio',1,'YSpace', [2 1.5 1]);
cmapD=autumn(length(frequencyRange));
lineStyles={'-','--',':','+-','.-'};
for dIdx=1:length(frequencyRange)
    goodLink=cellfun(@(x) x(1),jaccardCoeff{dIdx});
    plot(h,densityRange,goodLink,lineStyles{dIdx},'color',cmapD(dIdx,:),'LineWidth',2);
    axis(h,[min(densityRange),max(densityRange),0,1])
    title(h,'Dynamic quantification breakpoint')
    legend(h,num2str(frequencyRange'))
    xlabel(h,'Particle density')
    ylabel(h,'JSC')
    hold(h,'on');
end    

