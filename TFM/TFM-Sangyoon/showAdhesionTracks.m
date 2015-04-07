function [F]=showAdhesionTracks(pathForColocalization,idList,varargin)
% this function reads from Colocalization folder and shows the specified
% tracks on selected Channel.
%% input reading
ip =inputParser;
ip.addRequired('pathForColocalization',@ischar)
ip.addOptional('idList','all',@(x)islogical(x)||ischar(x)||isempty(x))
ip.addParamValue('numChan',2,@isscalar); % selcted track ids
ip.parse(pathForColocalization,idList,varargin{:});
pathForColocalization=ip.Results.pathForColocalization;
idList=ip.Results.idList;
numChan=ip.Results.numChan;

%% Read basic data
disp('Loading raw files ...')
tic
tracksNA = load([pathForColocalization filesep 'data' filesep 'tracksNA.mat'],'tracksNA');
tracksNA = tracksNA.tracksNA;
if numChan==1
    imgMap = load([pathForColocalization filesep 'fMap' filesep 'tMap.mat'],'tMap');
    imgMap = imgMap.tMap;
elseif numChan==2
    imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
    imgMap = imgMap.paxImgStack;
end
numFrames = size(imgMap,3);
if ischar(idList) && strcmp(idList,'all')
    startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA)));
    endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA)));
else
    startFrame = max(1, min(arrayfun(@(x) x.startingFrame,tracksNA(idList))));
    endFrame = min(numFrames, max(arrayfun(@(x) x.endingFrame,tracksNA(idList))));
end
toc
%% recording movie
curNumFrames = endFrame-startFrame+1;
F(curNumFrames) = struct('cdata',[],'colormap',[]);
p=0;
h=figure;
for ii=startFrame:endFrame
    p=p+1;
    % actual frame
    curFrame = imgMap(:,:,ii);
    imshow(curFrame,[]), hold on
    if ischar(idList) && strcmp(idList,'all')
        plot(arrayfun(@(x) x.xCoord(ii),tracksNA),arrayfun(@(x) x.yCoord(ii),tracksNA),'ro')
        xmat = cell2mat(arrayfun(@(x) x.xCoord(1:ii),tracksNA,'UniformOutput',false));
        ymat = cell2mat(arrayfun(@(x) x.yCoord(1:ii),tracksNA,'UniformOutput',false));
        plot(xmat',ymat','r')
    else
        plot(arrayfun(@(x) x.xCoord(ii),tracksNA(idList)),arrayfun(@(x) x.yCoord(ii),tracksNA(idList)),'ro')
        xmat = cell2mat(arrayfun(@(x) x.xCoord(1:ii),tracksNA(idList),'UniformOutput',false));
        ymat = cell2mat(arrayfun(@(x) x.yCoord(1:ii),tracksNA(idList),'UniformOutput',false));
        plot(xmat',ymat','r')
    end
    drawnow
    F(p) = getframe;
    hold off
end
close(h)
