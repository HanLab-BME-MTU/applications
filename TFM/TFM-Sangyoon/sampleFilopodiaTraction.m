function [] = sampleFilopodiaTraction(MD)
% sampleFilopodiaTraction samples traction along the lines that a user drew
%
% Sangyoon Han
% November 2022
%% ----------- Input ----------- %%

%% --------------- Initialization ---------------%%
% Set up the output directories
outputFilePath = [MD.outputDirectory_ filesep 'FilopodQuantification'];
imgPath = [outputFilePath filesep 'imgs'];
dataPath = [outputFilePath filesep 'data'];
tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];

if ~exist(figPath,'dir')
    mkClrDir(imgPath);
    mkClrDir(dataPath);
    mkClrDir(tifPath);
    mkClrDir(figPath);
end

nChans=numel(MD.channels_);
% e.g. see if there is TFM package
iTFM =  MD.getPackageIndex('TFMPackage');
% e.g. see if the additional channel is used for TFM or not.
%      if there is a channel, get the channel.
if ~isempty(iTFM)
    tfmPack = MD.getPackage(iTFM);
    forceProc = tfmPack.getProcess(4);
    dispProc = tfmPack.getProcess(2);
    
    dispFunParam = dispProc.funParams_; % forceProc.checkChannelOutput;
    iChanTFM = dispFunParam.ChannelIndex;
    % iAdhChan can be mistakenly input. It is reasonable to change it to
    % the other channel if it overlaps with iChanTFM
    if nChans==2
        iAdhChan = setdiff(1:nChans,iChanTFM);
        disp(['Adhesion channel index is automatically changed to ' num2str(iAdhChan) ' since ' ...
            num2str(iChanTFM) ' overlapped with it.']);
    else
            iAdhChan = 2;
    end
    % Read TFM map
    tMap=forceProc.loadChannelOutput('output','tMapUnshifted'); % in Pa per pixel (1pix x 1pix)
else
    disp('No TFM package run. Please run it and run this function again.')
    return
end


psfSigma = getGaussianPSFsigma(MD.numAperture_, 1, MD.pixelSize_*1e-9, MD.getChannel(iAdhChan).emissionWavelength_*1e-9);
if isempty(psfSigma)
    MD.getChannel(iAdhChan).emissionWavelength_ = 568; 
    disp('Adhesion channel emmission wave length is set to be 568 nm. If this value is incorrect, please set up your fluorophore of the adhesion channel')
    psfSigma = getGaussianPSFsigma(MD.numAperture_, 1, MD.pixelSize_*1e-9, MD.getChannel(iAdhChan).emissionWavelength_*1e-9);    
end
%% --------------- Sub-resolution object detection ---------------%%% 
jformat = ['%.' '3' 'd'];
% Changed it for isometric detection for nascent adhesion detection
pixSize = MD.pixelSize_;
convertL = pixSize/1000;
convertArea = (pixSize/1000)^2;
% minEcc = 0.7;

disp(['Adhesion channel was assumed to be in channel ' num2str(iAdhChan) '.'])
disp('Results will be saved under:')
disp(outputFilePath);

for j=2 %1:MD.nFrames_
    I=double(MD.channels_(iAdhChan).loadImage(j));
    curT = tMap(:,:,j);
    try
        iMaskProc = MD.getProcessIndex('MaskRefinementProcess','askUser',false);
        maskProc = MD.getProcess(iMaskProc);
%         mask = maskProc.loadChannelOutput(iPax,j);
        % if there are masks for more than one channels, combine them.
        
        if sum(maskProc.checkChannelOutput)==1
            iMaskChan = find(maskProc.checkChannelOutput);
            mask = maskProc.loadChannelOutput(iMaskChan,j);
        elseif sum(maskProc.checkChannelOutput)>1 && (~isempty(iTFM) && ~ismember(iChanTFM, find(maskProc.checkChannelOutput)))
            %Combine the the multiple masks to one
            maskEach = arrayfun(@(x) maskProc.loadChannelOutput(x,j),find(maskProc.checkChannelOutput),'UniformOutput',false);
            maskAll=reshape(cell2mat(maskEach),size(I,1),size(I,2),[]);
            mask = any(maskAll,3);
        elseif nChans==1 && nChans==iAdhChan
            mask = maskProc.loadChannelOutput(iAdhChan,j); % 1 is CCP channel
        elseif nChans>1 && ~isempty(iTFM) && ismember(iChanTFM, find(maskProc.checkChannelOutput))
            mask = maskProc.loadChannelOutput(iAdhChan,j); % 1 is CCP channel
        else
            mask = maskProc.loadChannelOutput(iAdhChan,j); % 1 is CCP channel
        end
    catch
        mask = true(MD.imSize_);
        noMask=true;
    end

    % Show the image with paxillin
    h0 = figure(1);
    imshow(I,[400 650]); hold on
    boundBand = bwboundaries(mask); 
    nBDs =numel(boundBand);

    for kk=1:nBDs
        boundary = boundBand{kk};
        plot(boundary(:,2), boundary(:,1), 'Color','w', 'LineWidth', 1) % cell boundary
    end

    % Ask user to draw lines.
    
%     disp('Draw a line from a tip of filopodium to its base. One at a time. Doubleclick when done per line.')
    disp('Pick several points along the filopodia from the tip to its base. One at a time. Doubleclick when done per line.')
    
    keepMeasuring = true;
    ii=0;
    
    while keepMeasuring
        ii = ii+1;
        figure(1)
        [xpoints, ypoints] = ginput;
        spcv = cscvn([xpoints, ypoints]');
        hold on
        fnplt(spcv)
        hold off
    %     roi = drawline('SelectedColor','yellow');
        knots = spcv.breaks(1):spcv.breaks(end);
    %     line_xs = xpoints(1):dx:xpoints(end);
        xy = fnval(spcv,knots);
        line_xs = xy(1,:);
        line_ys = xy(2,:);

    %     line_xs = [roi.Position(1,1), roi.Position(2,1)];
    %     line_ys = [roi.Position(1,2), roi.Position(2,2)];
        % extract pixel values along line
        pixvals = improfile(I,line_xs,line_ys);
        tracvals = improfile(curT,line_xs,line_ys);
        h = figure(2); 
        hold off
%         ax1 = axes;
        yyaxis left
        plot(1:length(pixvals),pixvals), axis tight, title('pixel values','FontSize',6)
        ylabel('pixel greyscale values','FontSize',6)
    %     xlabel('relative position from linestart (0) to lineend','FontSize',6)

        yyaxis right
        plot(1:length(pixvals),tracvals)
        ylabel('traction magnitude (Pa)','FontSize',6)
        xlabel('relative position from linestart (0) to lineend','FontSize',6)

        h.Units = 'inch';
        h.Position(3:4) = [2 1];
        savefig(h, [figPath filesep 'profile' num2str(ii,jformat) '.fig'])
    %     hT = findobj(ax1,'Type','text');
        % Tip detection
        tempI = I(round(ypoints(1))-9:round(ypoints(1))+9,round(xpoints(1))-9:round(xpoints(1))+9);
        prmVect = [0, 0, mean(I(:)), 3, min(I(:))];
        [prmVect] = fitGaussian2D(tempI,prmVect,'xyasc');
        tipLocX = xpoints(1)+prmVect(1);
        tipLocY = ypoints(1)+prmVect(2);
        %figure, imshow(I,[]), hold on, plot(tipLocX, tipLocY, 'ro')
        tipAmp = prmVect(3);
        bg = prmVect(5);
        tipForce = improfile(curT,tipLocX,tipLocY);
        % Base force
        % Will compare with cell mask
        idxInMask = arrayfun(@(i,j) maskVectors(i,j,mask),line_xs,line_ys);
        bgCell = min(I(mask));
        baseAmp = mean(pixvals(idxInMask))-bgCell;
        baseForce = mean(tracvals(idxInMask));
        % Shaft force: in between tip and base
        [idxWithinTipCell]=KDTreeBallQuery([tipLocX tipLocY],[line_xs' line_ys'],9*ones(numel(line_xs),1));
        idxWithinTip = ~cellfun(@isempty,idxWithinTipCell)';
        idxShaft = ~idxInMask & ~idxWithinTip;
        shaftAmp = mean(pixvals(idxShaft))-bg;
        shaftForce = mean(tracvals(idxShaft));
        
        % Storing into arrays
        % curve
        cvAll(ii) = spcv;
        intenProfile{ii} = pixvals;
        forceProfile{ii} = tracvals;
        lineXall{ii} = line_xs;
        lineYall{ii} = line_ys;
        
        tipLocXall(ii) = tipLocX;
        tipLocYall(ii) = tipLocY;
        shaftXall{ii} = line_xs(idxShaft);
        shaftYall{ii} = line_ys(idxShaft);
        baseXall{ii} = line_xs(idxInMask);
        baseYall{ii} = line_ys(idxInMask);
        
        ampTipAll(ii) = tipAmp+bg;
        forceTipAll(ii) = tipForce;
        ampBaseAll(ii) = baseAmp+bgCell;
        forceBaseAll(ii) = baseForce;
        ampShaftAll(ii) = shaftAmp+bg;
        forceShaftAll(ii) = shaftForce;
        
        doneDrawing = input('Done drawing ? (0/1):');
        if isempty(doneDrawing)
            doneDrawing = false;
        end
        
        figure(1), hold on
        plot(tipLocX, tipLocY, 'yo')
        plot(line_xs(idxInMask),line_ys(idxInMask),'r','LineWidth',1.5)
        plot(line_xs(idxShaft),line_ys(idxShaft),'w','LineWidth',1.5)
        
        if doneDrawing
            keepMeasuring = false;
        end
    end
    %% analysis of local maxima
    %Local maxes of forces
    lmT = locmax2d(curT,11);
    % Getting the locations
    [row, col] = find(lmT);
    lX = col;
    lY = row;
%     figure, imshow(lmT), hold on
    figure, imshow(curT,[]), colormap jet, hold on
    plot(lX,lY,'co')
    % Now, count how many positions are within reach to these loc maxima
    % tips
    [~, distFromLM]=KDTreeBallQuery([lX lY],[tipLocXall' tipLocYall'], 9*ones(numel(tipLocXall),1));
    TipsWithinLMs = ~cellfun(@isempty,distFromLM)';
    % shafts
    [~, distFromShaft]=cellfun(@(x,y) KDTreeBallQuery([lX lY],[x' y'], ...
        9*ones(numel(x),1)),shaftXall, shaftYall,'unif',false);
    ShaftsWithinLMs = cellfun(@(x) ~isempty(x) & any(cell2mat(x)<9),distFromShaft)';
    plot(shaftXall{1},shaftYall{1},'y','linewidth',2)    
    % bases
    [~, distFromBase]=cellfun(@(x,y) KDTreeBallQuery([lX lY],[x' y'], ...
        10*ones(numel(x),1)),baseXall, baseYall,'unif',false);
    BaseWithinLMs = cellfun(@(x) any(cell2mat(x)<12),distFromBase)';
    figure, plot(lX,lY,'ro'), hold on
    cellfun(@(x,y) plot(x,y,'bo'),baseXall,baseYall)
    h5 = figure;
    locmaxRatio = [sum(TipsWithinLMs)/length(TipsWithinLMs) ...
        sum(ShaftsWithinLMs)/length(ShaftsWithinLMs) ...
        sum(BaseWithinLMs)/length(BaseWithinLMs)];
    
    hb = bar(1:3,[sum(TipsWithinLMs)/length(TipsWithinLMs) ...
        sum(ShaftsWithinLMs)/length(ShaftsWithinLMs) ...
        sum(BaseWithinLMs)/length(BaseWithinLMs)]);
    nameList = {'Tip','Shaft','Base'};
    xticklabels(nameList)
    ylabel({'Ratio of region colocalized'; 'to traction local maxima (1)'},'FontSize',7)
    ylim([0 1])
    savefig(h5, [figPath filesep 'indep_forces.fig'])
    
% %         maskProc = MD.processes_{2};
% %         mask = maskProc.loadChannelOutput(1,j);
%     maskAdhesion = blobSegmentThreshold(I,minSize,0,mask & roiMask(:,:,j));
%     maskAdhesionFine = blobSegmentThreshold(I,0,0,mask & roiMask(:,:,j)); % mask for all adhesions without minSize
end
savefig(h0, [figPath filesep 'Overlay.fig'])

%% All forces
% allForces = {forceTipAll(forceBaseAll>65) forceShaftAll(forceBaseAll>65) forceBaseAll(forceBaseAll>65)};
allForces = {forceTipAll(forceTipAll>420) forceShaftAll forceBaseAll};
allTalins = {ampTipAll ampShaftAll ampBaseAll};
nameList = {'Tip','Shaft','Base'};
h2 = figure;
% yyaxis left
boxPlotCellArray(allForces,nameList,1,1,1,'forceShowP',1,'forceTtest',1)
ylabel('Traction (Pa)')
savefig(h2, [figPath filesep 'allForces.fig'])


h3 = figure;
% yyaxis right
boxPlotCellArray(allTalins,nameList,1,1,1,'forceShowP',1,'forceTtest',1)
ylabel('Talin-GFP amplitude (A.U.)')
ylim([300 800])
savefig(h3, [figPath filesep 'allAmps.fig'])
%% Normalized force
allForcesNormal = {forceTipAll./forceBaseAll forceShaftAll./forceBaseAll forceBaseAll./forceBaseAll};
figure
yyaxis left
boxPlotCellArray(allForcesNormal,nameList,1,1,1)

%% save
save([dataPath filesep 'allData.mat'])

end