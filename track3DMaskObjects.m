function movieData = track3DMaskObjects(movieData,paramsIn)
%TRACK3DMASKOBJECTS tracks the postions of objects in the input 3D movie masks
% 
% movieData = track3DMaskObjects(movieData,paramsIn)
% 
% movieData = track3DMaskObjects(movieData)
% 
% This function performs very simple (greedy nearest-neighbor) tracking of
% the objects in the masks for the input movie. This is designed to track
% small numbers of sparsely distributed objects which have already been
% segmented. This is VERY simple tracking, and requires that the number of
% objects in each mask be constant.
% 
% Input:
% 
%   movieData - A MovieData3D object describing the movie to be processed.
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       tracking results to.
%       If not input, the files will be saved to the same directory as the
%       movieData, in a sub-directory called "object_tracking"
%
%       ('ChannelIndex' -> Positive integer scalar) Optional. The integer
%       index of the channel to use masks from for tracking. If not input,
%       the first channel will be used.
%
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - the updated MovieData object with the parameters, paths
%   etc. stored in it, in the field movieData.processes_.
%
%   The tracking results are written to the directory specified by the
%   parameter OuptuDirectory, as a .mat file.
%
%
% Hunter Elliott, 
% 6/2010
%
%% ----- Parameters ----- %%

fName = 'object_tracks'; %String for naming result file

showPlots = false; %Show plots for debugging

%% ---- Input ---- %%

if ~isa(movieData,'MovieData3D')
    error('The first input argument must be a valid MovieData object for a 3D movie!')
end


%Check if the object tracking has been run before
iProc = movieData.getProcessIndex('MaskObjectTrackingProcess',1,0);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(MaskObjectTrackingProcess(movieData));
end

if nargin < 2
    paramsIn = [];
end

%Resolve input and default parameters
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

if numel(p.ChannelIndex)>1
    error('This function only supports single-channel tracking!')
end

%Make sure the movie has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,p.BatchMode);

if isempty(iSegProc)
    error('The movie must be segmented before the branches can be found!')
end

%Make sure the selected channel has valid masks
if ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex)
    error('The seclected channel does not have valid masks! Please check segmentation!');
end


%% ------ Init ------ %%

maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);
maskDir = movieData.processes_{iSegProc}.outMaskPaths_{p.ChannelIndex};
nFrames = movieData.nFrames_;

mkClrDir(p.OutputDirectory)

%Get pixel aspect ratio for distance calculations
pixAspect = movieData.zSpacing_ / movieData.pixelSize_;

%% ----- Detection & Tracking ----- %%



disp('Starting objects tracking...')
if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, tracking objects...');        
end        

for iFrame = 1:nFrames
    
    mask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
            
    %Since we expect sparse objects, label them using maximal connectivity
    %(which is default).
    if iFrame == 1
        [labelMask,nObj] = bwlabeln(mask);
        objTracks = zeros(nFrames,nObj,3);  
        
        if showPlots
    
            objCols = lines(nObj);
            figure
            hold on    
            xlim([0 size(mask,2)]);
            ylim([0 size(mask,1)]);
            zlim([0 size(mask,3)]);                      
            light
            view(-20,60)
        end
        
    else
        [labelMask,nObjCurr] = bwlabeln(mask);
        if nObjCurr ~= nObj
            error(['This function requires that the number of objects in each mask be constant! Different object number found in frame ' num2str(iFrame) ]);            
        end
    end
            
    
    %Get distance transform of mask to find "innermost" point of each
    %object. Remove mask pixels that touch the image boundary first.
    mask([1 end],:,:) = false;
    mask(:,[1 end],:) = false;
    mask(:,:,[1 end]) = false;    
    
    distX = bwdist(~mask);
            
    %Go through each object and find the point furthest from the perimiter,
    %called the "InnerMostPoint" This is used instead of the centroid
    %because I want to track the nucleus of cells with lots of branches etc.
    iNN = zeros(nObj,1);
    unusedObj = true(nObj,1);
    for j = 1:nObj
        
        %This is the lazy way of having the indices come out right when
        %using max - get the distance transform of only this object
        tmp = zeros(size(mask));
        tmp(labelMask == j) = distX(labelMask == j);
        
        [~,iInnermost] = max(tmp(:));
                
        [y,x,z] = ind2sub(size(mask),iInnermost);
        
        %Scale the z coordinate based on the pixel aspect ratio
        z = z * pixAspect;
        
        objTracks(iFrame,j,:) = [x y z];                
        
        if showPlots
           plot3(x,y,z,'.','MarkerSize',10,'Color',objCols(j,:))         
        end
        
        if iFrame > 1                        
            
            iObj = find(unusedObj);
            
            %Find the distance to each object in the previous frame
            [~,iNN(j)] = min(arrayfun(@(x)(sqrt(sum((objTracks(iFrame,j,:) - ...
                                              objTracks(iFrame-1,x,:)) .^2 ))),...
                                              iObj));  
            iNN(j) = iObj(iNN(j));
            unusedObj(iNN(j)) = false;
        
        end
        
        
    end
    
    if iFrame > 1
        objTracks(iFrame,:,:) = objTracks(iFrame,iNN,:);
    end
    
    if ~p.BatchMode && mod(iFrame,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iFrame / nFrames,wtBar)
    end

    
    
    
end



%% ------ Finish ----- %

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


save([p.OutputDirectory filesep fName '.mat'],'objTracks');


movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk


disp('Finished object tracking!')

