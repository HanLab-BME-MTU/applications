function specifyMorphoPoints(ML)


iChan = [1 3 4];%Channels to view for selection


nMov = numel(ML.movies_);

outFile = 'manual morphological feature identification.mat';

for j = 1:nMov
    
    
    currOutFile = [ML.movies_{j}.outputDirectory_ filesep outFile];
    
    if exist(currOutFile,'file')
        bOut = questdlg('Existing Feat. Ident. Overwrite?','Overwrite','Yes','No','No');
    else
        bOut = 'Yes';
    end
    
    currFig = imageViewer(ML.movies_{j},'ChannelIndex',iChan,'Saturate',.03);
    disp('Select the approximate center of cell "back"')
    backPt = impoint(get(currFig,'CurrentAxes'));
    disp('Select clockwise start of lamellipodium')
    lamPt1 = impoint(get(currFig,'CurrentAxes'));
    disp('Select clockwise end of lamellipodium')
    lamPt2 = impoint(get(currFig,'CurrentAxes'));
    
    a = msgbox('Close this when all points are satisfactory');
    uiwait(a);
    featIdent.Back = backPt.getPosition;
    featIdent.Lam1 = lamPt1.getPosition;
    featIdent.Lam2 = lamPt2.getPosition;
    
    hold on
    plot(featIdent.Back(1),featIdent.Back(2),'bx','MarkerSize',15);
    plot(featIdent.Lam1(1),featIdent.Lam1(2),'rx','MarkerSize',15);
    plot(featIdent.Lam2(1),featIdent.Lam2(2),'rx','MarkerSize',15);    
    
    featIdent.Time = datestr(now);
    save(currOutFile,'featIdent');
    
    close(currFig);
    
    disp(['Finished ' num2str(j) ' of ' num2str(nMov)])    
    
    
end

