function [ idxWindTest ] = GCASubRegionalGetWindowsMovie( movieData,varargin)
% GCASubRegionalGetWindowsMovie
%
%%
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addParameter('TSOverlays',true);
defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep  'GC' filesep 'masks'];
ip.addParameter('subRoiMaskDirIn',defaultInDir); % where the subRoimasks are stored.
ip.addParameter('windMethod','ProtrusionBased'); % what windowing method to use - currently assumes you will use one of the
% windowing methods associated with the MD.
defaultOutDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep 'GC' filesep 'windows'];
ip.addParameter('OutputDirectory_', defaultOutDir);
ip.addParameter('EndFrame',movieData.nFrames_-1);

ip.parse(varargin{:});
%% Initiate
% find the correct window output file.
idxWindProc =  cellfun (@(x) strcmpi(x.name_,'Windowing'),movieData.processes_);
windProcesses =  movieData.processes_(idxWindProc);
toLoad = cellfun(@(x) strcmpi(x.funParams_.MethodName, (ip.Results.windMethod)),windProcesses);
windDir = windProcesses{toLoad}.outFilePaths_;
%windDir = movieData.processes_{toLoad}.outFilePaths_;
listOfWindFiles = searchFiles('.mat',[],windDir,0,'all',1);
if ~isdir(ip.Results.OutputDirectory_)
    mkdir(ip.Results.OutputDirectory_);
end
%%
for iFrame = 1:ip.Results.EndFrame
    % load the windows for the current frame
    load(listOfWindFiles{iFrame});
    
    roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
    
    % load the subRoiMask
    
    [idxWindTestLog] = GCASubRegionalGetWindows(windows,roiMask);
    idxWindTest = find(idxWindTestLog);
    
    % get only winds that likewise existed in the preceding frame
    if iFrame>1
        % rewrite
        idxWindFinal=  intersect(idxWindTest , idxWindFinal);
    else
        idxWindFinal = idxWindTest;
    end
end
save([ip.Results.OutputDirectory_ filesep 'idxWindFinal.mat' ],'idxWindFinal');

% if plot make a folder to test the windows.
% find the window numbers common to all
% chec
if ip.Results.TSOverlays == true
    mkClrDir([ip.Results.OutputDirectory_ filesep 'TSOverlays']);
    % mkClrDir([ip.Results.OutputDirectory_ filesep 'WindowsAll']);
    filenames=movieData.getImageFileNames;
    
    for iFrame = 1:ip.Results.EndFrame
        img = double(imread([movieData.getChannelPaths{1} filesep filenames{1}{iFrame}]));
        roiMask  = logical(imread(([ip.Results.subRoiMaskDirIn filesep 'mask' num2str(iFrame,'%03d') '.tif'])));
        roiYX = bwboundaries(roiMask);
        load(listOfWindFiles{iFrame});
        [ny,nx] = size(img);
        setFigure(nx,ny,'on');
        
        cmap = colormap(jet(length(idxWindFinal)));
        imshow(-img,[]) ;
        hold on
        cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
        for iWind = 1:length(idxWindFinal)
            
            gcaPlotWindows(windows,{'k','FaceAlpha',0},'bandMax',1,'showNum',2);
            gcaPlotWindows(windows(idxWindFinal(iWind)),{cmap(iWind,:),'FaceAlpha',1},'bandMax',1);
            
        end
        saveas(gcf,[ip.Results.OutputDirectory_ filesep 'TSOverlays' filesep num2str(iFrame,'%03d') '.png']);
        close gcf
        %        setFigure(ny,nx,'on');
        %        imshow(-img,[]);
        %        hold on
        %        % plot all the windows as well as a check so can see the number
        %        gcaPlotWindows(windows,'bandMax',1,'showNum',5);
        %        saveas(gcf, [ip.Results.OutputDirectory_ filesep 'WindowsAll' filesep num2str(iFrame,'%03d') '.png']);
        
        
        %        close gcf
    end
end



end


