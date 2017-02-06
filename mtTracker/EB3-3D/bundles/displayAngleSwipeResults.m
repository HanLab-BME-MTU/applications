function displayAngleSwipeResults(fiberCell,randFiberCell)
    plotRange=[0,8];
    numAngle=length(fiberCell);
    [H]=setupFigure(2,ceil((numAngle+1)/2),numAngle+1,'AspectRatio',1,'AxesWidth',4);
    for i=1:numAngle
        %H=handles((2*i-1):2*i);
        handles=H(i);
        kinTracksCell={fiberCell{i},randFiberCell{i}};
        fiberCount=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints);
        kinetochoreCount=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints);
        
        for i=1:length(kinTracksCell)
            kinTracks=kinTracksCell{i};
            for k=1:length(kinTracks)
                fiberCount(i,kinTracks(k).f)=fiberCount(i,kinTracks(k).f)+double(length(kinTracks(k).fiber));
                kinetochoreCount(i,kinTracks(k).f)=kinetochoreCount(i,kinTracks(k).f)+1;
            end
        end
        plot(handles(1),linspace(0,kinTracksCell{1}.numTimePoints,kinTracksCell{1}.numTimePoints), fiberCount./kinetochoreCount);
        xlabel(handles(1),'Frame count');
        ylabel(handles(1),'avg MT per bundle');
        ylim(handles(1),plotRange);
        legend(handles(1),{'kin','rand'});
    end
    
end