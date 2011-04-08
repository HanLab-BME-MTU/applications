function cellMig_part_2_statistics(binPix, minTrackLength, timeWindow, dframes, dt, showTrackerMovie, showEdgeMovie, batchJob, doNuclei)
pixSize_um=0.64;
display(['The pixel size is set to ',num2str(pixSize_um),' um!']);


if nargin < 1
    % define parameters:
    binPix          =input('Specify the band width for density measurements,  binPix=[  100pix]: ');
    minTrackLength  =input('Set the minimal considered track length,  minTrackLength=[12frames]: ');
    timeWindow      =input('Time window for averaging the cell velocity, timeWindow =[ 1frames]: ');
    dframes         =input('Plot only velocities every ...                          =[24frames]: ');
    dt              =input('Specify the time between two frames as fraction of 1h   =[    1/6h]: ');
    showTrackerMovie=input('Do you want to plot a tracker movie?        Yes=1 / No=0=[      No]: ');
    showEdgeMovie   =input('Do you want to show the tracked sheet edge? Yes=1 / No=0=[      No]: ');

    if isempty(binPix)
        binPix=100;
    end
    if isempty(minTrackLength)
        minTrackLength=12;
    end
    if isempty(timeWindow)
        timeWindow=1;
    end
    if isempty(dframes)
        dframes=24;
    end
    if isempty(dt)
        dt=1/6;
    end
    if isempty(showTrackerMovie)
        showTrackerMovie=0;
    end
    if isempty(showEdgeMovie)
        showEdgeMovie=0;
    end
end

if nargin<8 || isempty(batchJob)
    batchJob=0;
end

try
    sglCell;
catch
    sglCell=0;
end

if nargin<9
    doNuclei=1;
end

try
    load('xDetectEdge.mat')
    load('xParameters.mat')
    if doNuclei
        load('xDetectNuclei.mat')
        load('xTracksNuclei.mat')
    end
catch exception
    display('Couldnt find all essential files, please browse to folder')
    return;
end

try
    load('xDetectEdgeHand.mat')
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');    
    display('!!!Overloaded detected edges by hand-drawn edges!!!'); 
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end


if doNuclei
    imageFileListNuclei=[pwd,filesep,'Nuclei'];
    try
        imageFileListNuclei=getFileListFromFolder(imageFileListNuclei);
    catch exception
        [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
            'Select First nuclei Image');
        
        if ~ischar(filename) || ~ischar(pathname)
            return;
        end
        
        imageFileListNuclei = getFileStackNames([pathname filesep filename]);
    end
end

try
    test=toDoList(1);
catch exception
    toDoList=1:length(sheetMask);
end

numFrames=toDoList(end);

if minTrackLength>numFrames
    minTrackLength=min(numFrames,minTrackLength);
    display(['Movie is too short, min. track length was set to: ',num2str(minTrackLength),' frames']);
end

if timeWindow>numFrames-1
    timeWindow=min(numFrames-1,timeWindow);
    display(['Movie is too short, time window was set to: ',num2str(timeWindow),' frames']);
    % Time window for averaging cell velocities. Time window=1 means
    % velocity calc between frame 1 to frame 2;
end

if dframes>=numFrames
    dframes=min(numFrames-1,dframes);
    display(['Movie is too short, plot results for every: ',num2str(dframes),' frames']);
end

try
    dPix;
catch
    dPix=zeros(toDoList(end),1);
end

if sheetMask(1).mat(dPix(1)+1,dPix(1)+1)==1
    signFac=1;
    display('The sheet moves to the right')
else
    signFac=-1;
    display('The sheet moves to the left')
end

%**************************************************************************
% Here starts the display                                                 *
%**************************************************************************

% Show and save the movie:
if showTrackerMovie==1
    % show and save movie:
    movieFileName='mov_trackNuclei.mov';
    
    display('Create tracker movie:...')
    % Set the position of the current figure to screen size:
    scrsz = get(0,'ScreenSize');
    h     = figure();
    set(h,'Position',scrsz);
    
    overlayTracksMovieNew(tracksFinal,[],10,1,movieFileName,[],1,0,0,[],0,1,[],0,0,imageFileListNuclei)
    % move the movie into the current dir:
    if isdir('Nuclei')
        try
            movefile(['Nuclei',filesep,movieFileName],movieFileName);
        catch exception
            display('Couldnt find movie file, you have to find it yourself');
        end
    end
end
display('done!')

% Plot the sheet edge:
if showEdgeMovie==1
    display('Plot the sheet edge:...')
    plotSheetEdge([],sheetBnD,sheetEdge,toDoList,[],batchJob)
    % plotSheetEdge([],[],[],[],[],batchJob)
    display('done!')
end

% calculate the (averaged) traveled distance. This function also performs a
% check if the sheet propagation is not due to edge misdetection:
display('Calculate covered area over time:...')
toDoListAll=toDoList;
[coveredDist,coveredArea,toDoList,badFlag]=calcCoveredArea(sheetMask,dPix,toDoList,1);
if badFlag
    save(['xxBadFlagAtFrame',num2str(toDoList(end)),'.mat'],'toDoList','toDoListAll','-v7.3');
    numFrames=toDoList(end);
end
clear toDoListAll;
display('done!')

if doNuclei
    display('Calculating the distance transform:...')
    %[cellDistFromEdge,distGrad]=createCellDistMat(sheetMask,r,dPix,toDoList);
    [cellDistFromEdge,distGrad]=createCellDistMat(sheetMask,r,dPix,toDoList,sglCell);
    display('done!')
    
    % Calculate the cell area at front and back. This has to be done before
    % filtering out short tracks:
    display('Perform measurement of the cell density:...')
    [tracksMatxCord,tracksMatyCord]=conv2CordMatVelMat(tracksFinal,[],1,1,toDoList);
    [densityMeasurement]=perfDensityMeasurement(cellDistFromEdge,tracksMatxCord,tracksMatyCord,binPix,toDoList);
    display('done!')
    
    
    % sort out good tracks (defined by minimal track length), convert to simple
    % x-y-coordinate matrix, calculate velocities:
    display('Filter tracks for minimal track length:...')
    [tracksMatxCord,tracksMatyCord,velMatxCord,velMatyCord,velMatMag]=conv2CordMatVelMat(tracksFinal,[],minTrackLength,timeWindow,toDoList);
    display(['Calculated velocities are averaged over: ',num2str(timeWindow),' frames!'])
    display('done!')
    
    % calculate the cell to edge distance for all frames.
    display('Calculate individual cell to edge distances:...')
    [cell2edgeDist]=calcCell2EdgeDist(tracksMatxCord,tracksMatyCord,cellDistFromEdge,toDoList);
    display('done!')
    
    % calculate the distance gradient for each nuclei:
    display('Calculate distance gradients at each cell position:...')
    [cellDistGradx,cellDistGrady]=calcCellDistGrad(tracksMatxCord,tracksMatyCord,distGrad,toDoList);
    display('done!')
end

% calculate the edge roughness and the frequency spectrum:
display('Calculate the edge roughness:...')
[roughness,freqSpec]=calcEdgeRoughness(sheetEdge,r,toDoList);
display('done!')

%**************************************************************************
% Here starts the output of the results                                   *
%**************************************************************************

display('Display results:...')

% figure(1)
% for frame=1:dframes:(numFrames-timeWindow) % 1:(numFrames-timeWindow) %[1 numFrames-timeWindow]
%     marker=['ro','bs','m*','c+','gd','yo','ks'];
%     %plot(tracksMatxCord(:,frame),  velMatMag(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
%     plot(tracksMatxCord(:,frame),velMatxCord(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
%     hold on
% end
% hold off

% figure(2)
% frame=numFrames;
% imagesc(cellDistFromEdge(frame).mat);
% hold on
% quiver(tracksMatxCord(:,frame),tracksMatyCord(:,frame),cellDistGradx(:,frame),cellDistGrady(:,frame));
% hold off

if doNuclei
    
    figure()
    label=1;
    for frame=1:dframes:(numFrames-timeWindow)
        marker=['ro','bs','m*','c+','gd','yo','ks'];
        %Remove the NaNs:
        %plot(tracksMatxCord(:,frame),
        %velMatMag(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        plot(cell2edgeDist(:,frame),signFac*velMatxCord(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        hold on
        M{label}=[num2str((frame-1)*dt),'h'];
        label=label+1;
    end
    xlabel('Distance to the wound edge in [pixel]')
    ylabel(['x-component of the velocity [pixel/',num2str(timeWindow*dt),'h]'])
    legend(M);
    saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn');
    clear M
    hold off
    
    % figure(4)
    % % project the cell velocity on the distance gradient (This doesn't help!):
    % for frame=[2 6]
    %     marker=['ro','bs','m*','c+','gd','yo','ks'];
    %     %plot(tracksMatxCord(:,frame),  velMatMag(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
    %     plot(cell2edgeDist(frame).vec,sum(cellDistGrad(frame).vec.*[velMatxCord(:,frame) velMatyCord(:,frame)],2),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
    %     hold on
    % end
    % hold off
    
    figure()
    % plot how the cell density evolves over time in each bin:
    label=1;
    for frame=1:dframes:numFrames
        marker=['ro','bs','m*','c+','gd','yo','ks'];
        %plot(tracksMatxCord(:,frame),  velMatMag(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        [numBins,~]=size(densityMeasurement.density);
        % disregard the last two bins:
        dBin=0;
        plot(densityMeasurement.binPix*(1:numBins-dBin),densityMeasurement.density(1:numBins-dBin,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        hold on
        M{label}=[num2str((frame-1)*dt),'h'];
        label=label+1;
    end
    display(['Disregarded the last',num2str(dBin),' bins for the density measurement'])
    xlabel('Distance to the wound edge in [pixel]')
    ylabel('Cell density in [1/pix]')
    legend(M);
    saveas(gcf,['fig_cell_density_over_dist_to_edge','.tiff'],'tiffn');
    clear M
    hold off
    
end

% plot how the cell density evolves over time in each bin:
figure()
plot(toDoList,coveredDist);
xlabel('Time [frame number]')
ylabel('Traveled distance of the sheet edge [pix]')
saveas(gcf,['fig_distance_over_time','.tiff'],'tiffn');


display('done!')

% save only the important results:
if doNuclei
    save('xResults.mat','signFac','cell2edgeDist','tracksMatxCord','tracksMatyCord','velMatxCord','velMatyCord','velMatMag','densityMeasurement','coveredDist','roughness','freqSpec','minTrackLength','timeWindow','dt','pixSize_um','toDoList','badFlag','-v7.3');
else
    save('xResults.mat','signFac','coveredDist','roughness','freqSpec','minTrackLength','timeWindow','dt','pixSize_um','toDoList','badFlag','-v7.3');
end


% To do:
% 2) Persistence of cell movement (direct link over traveled length)
% 3) Cordination between cells (see PNAS).