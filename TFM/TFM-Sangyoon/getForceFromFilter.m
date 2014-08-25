function [tracksNA] = getForceFromFilter(tracksNA,MD,cropPosition, outputPath,tmax,r1,r2)
% [] = getForceFromFilter(MD,tracksNA,iChan)
% obtain the average force magnidue in the designed fileter per a track and
% save them in the current folder
% 
% input:    
%           tracksNA            tracks
%           MD                     movieData
%           cropPosition        this contains crop-position info
%           outputPath          outputPath designated when running
%                                       colocalizationAdhesionWithTFM.m
%           tmax                max force value used in colocalizationAdhesionWithTFM
%                                   function
%           r1,r2                inner (r1) and outer (r2) radius of filter around the center of track from
%                               which the force value and gathered and averaged
%          
% output:  tracksNA          tracks that contains intensity information in
%                                           forceAround.

% Sangyoon Han August 2014
% exampel: tracksNA =
% getForceFromFilter(tracksNA,MD,cropPosition,'.',10,20)
%% Data Set up
% Get whole frame number
nFrames = MD.nFrames_;
% % name of the intensity of iChan
% nameChan = ['amp' num2str(iChan)];
%% Get the image stack
pathForTheMovieDataFile = MD.getPath;
outputFilePath = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
forcemapPath = [outputFilePath filesep 'fMap'];
iiformat = ['%.' '3' 'd'];

for ii=1:nFrames
    % Get the image in forcemapPath
    tsMap = imread(strcat(forcemapPath,'/force',num2str(ii,iiformat),'max',num2str(tmax),'.tif'));
    forceImg(:,:,ii) = double(tsMap)*3500/(2^15); %converting to Pa
end
%% Get the image in iChan
% get the same xmat and ymat with forceImg
xmin = cropPosition(1);
xmax = cropPosition(1)+cropPosition(3);
ymin = cropPosition(2);
ymax = cropPosition(2)+cropPosition(4);
[x_mat, y_mat]=meshgrid(xmin:xmax,ymin:ymax);

for k=1:numel(tracksNA)
    iStart = tracksNA(k).startingFrame;
    iEnd = tracksNA(k).endingFrame;
    for ii=iStart:iEnd
        % Get the intensity value from tracks
        if ~tracksNA(k).presence(ii) && iStart >= tracksNA(k).startingFrame
            continue
        end
        if ii<tracksNA(k).startingFrame
            p = tracksNA(k).startingFrame;
        else
            p = ii;
        end
        curX = tracksNA(k).xCoord(p);
        curY = tracksNA(k).yCoord(p);
        
        curTrackAroundMask = sqrt((x_mat-curX).^2+(y_mat-curY).^2) > r1 & sqrt((x_mat-curX).^2+(y_mat-curY).^2)  <= r2 ;
        aroundIdx = find(curTrackAroundMask);
        [nrow,ncol] = ind2sub(size(curTrackAroundMask),aroundIdx);
        curAllStress=zeros(length(aroundIdx),1);
        for jj = 1: length(aroundIdx)
            curAllStress(jj) = forceImg(nrow(jj),ncol(jj),ii);
        end
        tracksNA(k).forceAround(ii)  = mean(curAllStress);
        tracksNA(k).forceAroundStd(ii)  = std(curAllStress);
    end
end


