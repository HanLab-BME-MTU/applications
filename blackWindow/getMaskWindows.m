function windows = getMaskWindows(maskIn,perpSize,paraSize,varargin)
%GETMASKWINDOWS sub-divides the input 2D mask into "edge-centric" polygonal sampling windows
% 
% Syntax:
%
% windows = getMaskWindows(maskIn,perpSize,paraSize)
% windows = getMaskWindows(...,'OptionName1',optionValue1,'OptionName2,optionValue2,...)
% windows = getMaskWindows(maskIn,perpSize,[],'OptionName1',optionName2)
%
%
% Description: 
%
% The object in the input mask is subdivided into polygonal sampling
% windows. These windows are called "edge-centric" because each column or
% "strip" of windows is associated with a specific portion of the mask
% edge. This association is based on proximity; that is, for each column of
% windows, the portion of the edge they are associated with is the closest
% edge on the mask. Additionally, each row or "band" of windows in the
% returned matrix occupies the same range of distances from the edge of the
% object in the input mask. Similarly, window n + 1 in a given column, or
% "strip" of the window matrix is further away from the edge of the object
% than window n is. More specifically, windows in a given row will all
% cover areas of the cell with the same distance-transform values while
% ascending a given column the distance transform values will increase in a
% maximal-gradient-ascent sense.
% 
% NOTE: The mask must be 2D, have only one object in it, and this object
% may not have holes!
%
%
% Examples: 
% 
% The command
% 
%   windows = getMaskWindows(maskIn,10,15);
% 
% would return a cell-array of windows whose size moving into the mask
% ("perpindicular" to the edge) was 10 pixels, and whose size along the
% mask edge ("parallel" to the edge) was 15 pixels.
% 
% Alternatively, the command
% 
%   windows = getMaskWindows(maskIn,10,[],'NumParallel',20);
% 
% would return a cell-array of windows whose size moving into the mask
% ("perpindicular" to the edge) was 10 pixels, and which contained 20
% strips/columns of windows whose size "parallel" to the mask edge was
% determined by dividing the total length of the mask boundary by 20.
% 
% For more alternative methods of specifying window dimensions, see the
% descriptions of the perpSize and paraSize inputs, and the NumParallel and
% StartPoint options.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Currently, the windowing excludes areas of the object which do not
% strictly fit the distance-transform based criteria described above. This
% means that areas near the image boundary, or which contain ridge-lines of
% the distance transform, will not be windowed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% Required Inputs:
% 
%   maskIn - 2D Binary matrix to be sub-divided into sampling windows. The
%   mask must contain only ONE object to be windowed - that is, it must
%   contain only one topologically distinct area of 'true' values - and
%   this object must not have holes in it.
% 
%   perpSize - Either a scalar specifying the size (in pixels) of each
%   window in the direction perpindicular to the edge of the object, if all
%   windows are to be the same size, or a vector with the size of each
%   window individually (in pixels), if the windows are to have varying
%   size. That is, the strip of windows at the mask edge will occupy
%   distances from the edge of 0 < x <= PerpSize pixels. The second strip
%   will occupy distances of PerpSize+1 < x <= 2xPerpSize and so on. 
%   Default value is 10 pixels.
%   
% Optional Input:
%
%   paraSize - Either a scalar specifying the size (in pixels) of the
%   windows in the direction parallel to the edge of the object if all
%   windows are to be the same size, or a vector specifying the widths of
%   each window strip individually (in pixels), if the windows are to have
%   varying size. Note that the windows will only be exactly this width at
%   the distance from the mask edge specified by the option StartContour.
%   Their width elsewhere within the mask will vary based on it's shape and
%   therefore the structure of the distance transform. It is geometrically
%   impossible for them to be the same width at every distance from the
%   mask edge unless the mask is a perfect rectangle.
%   **NOTE***: If this input is omitted, you must specify either the
%   NumParallel option or an array of start points using the StartPoint
%   option.
%   **NOTE***: If the distance transform isocontour specified by
%   StartContour cannot be evenly divided usin the specified paraSize, then
%   an extra strip of windows will be added to occupy the remaining area.
% 
% 
% Options: 
%
%   These are input as optionName/value pairs:
%   ('OptionName'->possible values)
%
%   ('StartPoint'->Mx2 matrix) If input as a 1x2 vector, this specifies the
%   m,n coordinates of the "origin" of the window polygons in the mask.
%   This will determine the location of the first strip of windows. If the
%   mask object is not completely contained in the image, the start point
%   will always be at the edge of the image and this input will be ignored.
%   Alternatively, if input as a Mx2 matrix, where M>1, this specifies the
%   start points of each of the M strips of windows, and therefore also
%   their widths in the direction parallel to the mask edge. If the
%   location(s) specified by the StartPoint option are not on the mask
%   boundary, the closest point on the mask boundary to each specified
%   point will be used.
%   No default. If not input, the start point will be wherever countourc.m
%   decides to put it. (Usually leftmost point in image).
% 
%   ('NumParallel'->Positive integer scalar) This option is an alternative
%   to specifying the window size in the direction parallel to the object
%   edge with the paraSize input. This option instead specifies the number
%   of window strips. The size of each strip is therefore equal to the
%   length of the object boundary divided by this number. Note that if the
%   paraSize parameter is input, this option will be ignored.
%   No default. If not input, another method of specifying the strip width
%   must be input (paraSize, or StartPoint array);
%
%   ('StartContour'->Positive integer scalar) Unless the mask is perfectly
%   rectangular, the widths of each window in a given strip will vary
%   within the mask (think pizza slices, then think of trying to slice a
%   pizza with a really fucked-up shape). As a result, the window strip
%   widths specified by the user can only be exactly correct at ONE
%   distance from the mask border - one isocontour of the distance
%   transform. This option specifies that isocontour. That is, if
%   StartContour is input as 1, the windows will be the specified width
%   at outer border of the first window band, i.e. the mask edge. If it is
%   input as 2, the windows will be the specified widths at the outer
%   border of the second band of windows, and so on.
%   Default is 2.
%
%   ('ShowPlots'->true/false) If true, a figure will be displayed showing
%   the distance transform, its isocontours, the gradient ascent lines, and
%   the window polygons. This figure is slow to draw, so turning on this
%   option may drastically increase the time required to window.
%
%   ('DoChecks'->true/false) If true, the topology of the mask will be
%   checked prior to starting the windowing to make sure there is only one
%   object and it has no holes in it. Only turn this off if you are going
%   to thoroughly check the masks before passing them to this function. In
%   that case, it can speed processing slightly. 
%   Default is true.
%
%   ('StartPointsOnly'->true/false) If true, the windowing will not be
%   performed and instead the start point locations (the locations of the
%   gradient ascent starts) will be returned returned as an Mx2 matrix.
%   This is used by getMovieWindows.m for propagating these start points
%   and is probably useless for anything else.
%
% Output:
%
%   windows - A 3-layer cell array containing the 2xM matrices giving the
%   x-y coordinates of the polygonal sides for each window. Briefly, this
%   means that for windows{j}{k}{l}, j is the position along the object
%   edge, and k is the position into the object (higher k is further
%   inward), and l is the side number.
%
%       In more detail:
%
%   By "3-layer" I mean a cell-array-of-cell-arrays-of-cell-arrays. I know,
%   it sounds really complicated but bear with me here. The first layer
%   corresponds to position along the edge of the mask object. That is, all
%   the elements of windows{j}{:} correspond to the j-th strip of windows,
%   all of which are associated with a single section of the mask edge. The
%   second layer corresponds to the position away from the edge of the mask
%   object. That is, the window specified by windows{j}{k} is in the j'th
%   strip of windows and the k-th band. All windows with a given k will
%   occupy the same distances from the mask edge. The third layer contains
%   the different sides of the window; that is, windows{j}{k} is a 1xL cell
%   array, where L is the number of sides to the window (either 3 or 4).
%   For L = 4, side 1 is on the isocontour closer to the object edge, side
%   2 follows the gradient ascent to the next contour, side 3 is on the
%   isocontour further from the object edge, and side 4 follows the
%   gradient ascent downhill, back to the contour closer to the object
%   edge. Windows where L=3 are created by the intersection of two gradient
%   ascent lines at a ridge, and therefore have no side on the isocontour
%   further from the object edge. The coordinates of each side are
%   specified by a 2xN matrix of x-y coordinates, where N may be different
%   for each side.
%   The window array MAY contain empty windows - these are areas of the
%   mask that were excluded, and the empty windows are left as
%   place-holders.
%
% Known Issues (areas excluded from windowing):
%
%   -When StartContour is set to 1, and there are small features in the
%   mask approaching the size of ~perpSize in radius, these areas and areas
%   uphill of them on the distance transform will not be windowed.
%   -If a Windows would cross a ridgeline in the distance transform that
%   would cause part of the window to actually be closer to an edge on the
%   other side of the object, this window is excluded.
%   -When the object to be windowed touches the image border, there will be
%   areas near the image boundary which will not be windowed.
%
% Hunter Elliott 
% Re-Written 4/2010
%

%% ---------------- Input ---------------- %%

if nargin < 2 || isempty(maskIn) || isempty(perpSize)
    error('You must include at least 2 inputs: A mask, "maskIn" and a perpindicular size "perpSize"!')
end

if ndims(maskIn) ~= 2 || ~islogical(maskIn)    
    error('The first input "maskIn" must be a 2D, binary matrix!')    
end

[M,N] = size(maskIn);

if any(perpSize < 1)
    error('The perpindicular size specified by the perpSize input must be > 1 pixel!')
end

if any(paraSize < 1)
    error('The parallel size specified by the paraSize input must be > 1 pixel!')
end

%Parse additional options
[startPoint,nPara,startContour,showPlots,doChecks,spOnly] = parseInput(varargin);

if ~isempty(nPara) && (numel(nPara) > 1 || nPara < 2 || round(nPara) ~= nPara)
    error('The number of windows specified by the nPara input must be a positive integer scalar >= 2!')
end

if ~isempty(startPoint) 
    if size(startPoint,2) ~= 2 
        error('The startPoint input must be a Mx2 matrix where M is the number of start points!')
    end
    nStartPts = size(startPoint,1);
else
    nStartPts = 0;
end

%Make sure the user specified the parallel size in at least one of the
%possible ways
if isempty(paraSize) && nStartPts < 2 && isempty(nPara)
    error('You must specify the size of the window strips by either the paraSize input, the nPara input, or by a StartPoint array!')
end

if isempty(startContour)
    startContour = 2;
elseif numel(startContour) > 1 || ~isequal(abs(round(startContour)),startContour)
    error('The startContour option must be a positive integer scalar!')
end

if isempty(doChecks)
    doChecks = true;
end

%We allow mask-checking to be disabled in case masks are pre-checked.
if doChecks
    %Check mask for holes, too many objects.
    CC = bwconncomp(maskIn);        
    if CC.NumObjects > 1 || ~isequal(maskIn,imfill(maskIn,'holes')) ...
            || nnz(maskIn) == 0 || nnz(maskIn) == numel(maskIn)        
        error('Problem with input mask: Mask must contain only one object, and this object must not have holes in it!')
    end
end

if isempty(showPlots)
    showPlots = false;
end

if showPlots
    firstTime = true;
end

if isempty(spOnly)
    spOnly = false;
end

%% --------------- Parameters ----------------- %%
%Hard-coded parameter values. These should NOT be changed unless you really
%know what you are doing.

smNum = .1;%Small number for recognizing zero-value contours - the contourc.m function replaces zero contour values with an arbitrary small number.
collapsedSize = .1; %Small number for recognizing collapsed window strips. If vertex separation along the contours is less than this, thew window is considered collapsed.
gaShift = .51;%How far to shift points near the image edge to avoid immediate termination of gradient ascent.
intErr = .6;%Number for checking if contour intersections are approximately increasing in index. This number allows for numerical error in intersection locations in windows which have partially collapsed.

%% ---------------- Contouring ---------------- &&
%Sub-divides the mask object 'perpindicular' to their edges by finding
%isocontours of the distance transform.

%Get distance transform of mask
distX = bwdist(~maskIn);
maxDist = max(distX(:));


% --- Get IsoValues ---- %

%Get isovalues of distance transform from window size perpSize
if length(perpSize) > 1
    %If a vector of sizes was entered...

    %Check all the sizes
    if min(perpSize) < 1
        error('Window size(s) as specified by perpSize must all be >= 1!!')
    end    
        
    %Convert these sizes into distance-transform isovalues
    distXvals = cumsum([0 perpSize(:)']);
    
    %Remove values that are too high
    distXvals = distXvals(distXvals <= maxDist);
    
else
    %If a single size was entered, convert it into distance-transform
    %isovalues
    distXvals = 0:perpSize:maxDist;
end

% --- Get Contours ---- %

%Find the isocontours of the distance transform at these values
contours = contourc(double(distX),double(distXvals)); %contourc only takes in doubles

% ---- Post-Process Contours ---- %

%Seperate the contours into individual matrices in a cell-array
[contours,contourValues] = separateContours(contours); %We need to retrieve the values also, because a given value may have more than 1 contour
    
%Clean up these contours...
[contours,iKept] = cleanUpContours(contours,4);

%Remove any contour values for contours which were removed in clean-up
contourValues = contourValues(iKept);

%The contourc.m function replaces zero contour values with an arbitrary
%small number (4e-6). Convert these values back to zero
contourValues(contourValues < smNum) = 0;

%We should have only one zero-value contour. If not, this means the mask
%contains more than one object, or the one object has a hole in it.
if nnz(contourValues == 0)>1
    error('Problem with input mask: Mask cannot contain more than one object, and that object may not have holes in it! Check mask!');
else
    iZeroCont = find(contourValues==0);
end

%Determine which contours are closed / open. We check all because if the
%zero level contour is closed, all the contours will be closed. However if
%the zero-level contour is open, we may still have closed higher-level
%contours.
%Additionally, we consider a contour to be closed if its start and endpoint
%are within the same pixel. This is because contourc.m occasionally returns
%start and end points which are very close but not identical (due to
%numerical error?) and because it is impractical to have an open contour
%whose endpoints are less than a pixel apart.

%Get distance between start and end points
startEndDist = cellfun(@(x)(sqrt(sum(diff(x(:,[1 end]),1,2) .^2))),contours);
isClosed = startEndDist < 1;
%Because we allow this leeway in considering a contour closed, we
%completely close them by setting the last and first points equal.
contours(isClosed & startEndDist > 0) = cellfun(@(x)(x(:,[1:end-1 1])),contours(isClosed & startEndDist > 0),'UniformOutput',false);

%If the contour is open, we need to fill in the boundary before checking
%handedness
closedContours = contours;
if any(~isClosed)
    closedContours(~isClosed) = closeContours(contours(~isClosed),distX);        
    isClockwise = cellfun(@(x)(isCurveClockwise(x)),closedContours);
else    
    isClockwise = cellfun(@(x)(isCurveClockwise(x)),contours);
end

%Flip the curves that are counter-clockwise so they run clockwise
contours(~isClockwise) = cellfun(@(x)(x(:,end:-1:1)),contours(~isClockwise),'UniformOutput',false);


%Find the contours which have the start isocontour value
iPerpStart = find(contourValues == distXvals(startContour));
nStart = numel(iPerpStart);


%% --------- Slicing --------- %%
%Divides the cell up perpindicular to it's edge via maximal-gradient-ascent
%of the distance transform.

% ----- Find / Set start point ---- %

%Check if the user specified a start point
if ~isempty(startPoint)        
                            
        %Anonymous function for calc distance to contour
        dToContour = @(x)(sqrt( (x(1,:) - startPoint(1,1) ) .^2 + ... 
                                (x(2,:) - startPoint(1,2) ) .^2 ));

        iParaStart = ones(1,nStart);

        %Find the closest point on all the contours of this value, which will
        %be used as the starting point.
        [~,iParaStart(isClosed(iPerpStart))] = cellfun(@(x)(min(dToContour(x))),...
                                    contours(iPerpStart(isClosed(iPerpStart))));                                                               
                
else            
    iParaStart = ones(1,nStart); %Default is to start at first point of contour(s)   
end


% ---- Determine number and location of slices ---- %

%Loop through contours with this value
distVals = cell(1,nStart);
nWinPara = zeros(1,nStart);
paraVertInd = cell(1,nStart);
paraVert = cell(1,nStart);
slices = cell(1,nStart);

%Calculate length along each start contour beforehand in case nPara rather
%than paraSize was specified
distAlong = cellfun(@(x)([0 cumsum(sqrt(diff(x(1,:)) .^2 + diff(x(2,:)) .^2))]),contours(iPerpStart),'UniformOutput',false);

if nStart > 1 && isempty(paraSize) && ~isempty(nPara)                    
        
        %If the user specified nPara rather than paraSize, and there are
        %multiple start contours, we need to  apportion this number among
        %the start contours according to their length.
        
        %Calculate total start-contour length
        totLen = sum(cellfun(@(x)(x(end)),distAlong));

        %Get fraction of total length corresponding to each contour
        fracLen = cellfun(@(x)(x(end) / totLen),distAlong);

        %Apportion the slices among the start contours according to this
        %fraction, while maintaining the total
        nPara = apportionIntegers(fracLen,nPara);    

        %Remove start contours which are too small for even a single window
        tooShort = nPara == 0;    
        nStart = nStart - nnz(tooShort);
        nPara(tooShort) = [];
        iPerpStart(tooShort) = [];
        
end

if isempty(paraSize) && nStartPts > 1    
    if nStart > 1
    
        %If the user specified a start point array and there are multiple
        %starting contours, we need to assign each start point to the
        %contour it is closest to.                
        iClosestCont = zeros(nStartPts,1);
        %Find which contour each start point is closest to
        for j = 1:nStartPts
            
            [~,iClosestCont(j)] = min(cellfun(@(x)(min(sqrt( ...
                               (x(1,:) - startPoint(j,1)) .^2 + ...
                               (x(2,:) - startPoint(j,2)) .^2))),...
                               contours(iPerpStart)));
                        
        end 
        %Split the start points up among the contours        
        spArray = cell(1,nStart);
        nStartPerCont = zeros(1,nStart);
        for j = 1:nStart
            nStartPerCont(j) = nnz(iClosestCont == j);
            spArray{j} = startPoint(iClosestCont==j,:);            
        end
        
        %Remove start contours which do not have any start points
        %associated with them
        noPoints = nStartPerCont==0;
        nStart = nStart - nnz(noPoints);
        nStartPerCont(noPoints) = [];
        iPerpStart(noPoints) = [];        
        
    else
        spArray = {startPoint};
        nStartPerCont = nStartPts;                        
    end
end

for i = 1:nStart

    %Calculate cumulative distance along the contour
    nPts = length(contours{iPerpStart(i)});
    if iParaStart(i) > 1 %Get the indices using the new start-point, and make sure the contour is completely closed
        contours{iPerpStart(i)} = contours{iPerpStart(i)}(:,[iParaStart(i):nPts 1:iParaStart(i)]);  
        nPts = nPts+1; %Extra point created by closure.
    end

    %Determine window slice locations via one of the possible inputs
    if ~isempty(paraSize)
        %If the user specified a strip size...
        if length(paraSize) == 1
            distVals{i} = 0:paraSize:distAlong{i}(end);                        
        else
            %If the user specified an array of sizes for the strips
            distVals{i} = [0 cumsum(paraSize)];
            %Remove those which go past the end of  the contour.
            distVals{i} = distVals{i}(1:(find(distVals{i}>distAlong{i}(end),1)-1));            
        end
        %If the total contour length is not an even multiple of
            %paraSize, we add an extra slice to close the gap.
        if abs(distVals{i}(end) - distAlong{i}(end)) > collapsedSize
                distVals{i} = [distVals{i} distAlong{i}(end)];
        end            
        %Determine the number of windows there will be parallel to this edge 
        nWinPara(i) = numel(distVals{i})-1; %If the contours are closed, there will be duplication of the vertices on the "seam", but who cares?

    elseif ~isempty(nPara) && nPara(i) > 0
        %If the user specified the number of windows...                
        distVals{i} = linspace(0,distAlong{i}(end),nPara(i)+1);
        nWinPara(i) = nPara(i);
        
    elseif nStartPts > 1
        %If the user specified a start-point array...                
        nWinPara(i) = nStartPerCont(i)-1;
        %If the user specified the slice start points locations directly...
        for j = 1:min(nStartPerCont(i),nPts)%In case the user specified more start points than points on this contour...Anything can happen....
            dToContour = @(x)(sqrt( (x(1,:) - spArray{i}(j,1) ) .^2 + ... 
                                    (x(2,:) - spArray{i}(j,2) ) .^2 ));

            [~,iClosest] = cellfun(@(x)(min(dToContour(x))),...
                                    contours(iPerpStart(i)));
            distVals{i}(j) = distAlong{i}(iClosest);
        end
        %Make sure these are in ascending order                

        %First, if this is a closed contour we may need to change the last
        %index. This is because for start points near the "seam" in the windows
        %min will return the first point rather than the last, even though
        %they are equidistant
        if isClosed(iPerpStart(i)) && distVals{i}(end) == 0
            distVals{i}(end) = distAlong{i}(end);
        end
        %Now we can sort the values safely
        distVals{i} = sort(distVals{i});        
    end                            
        
    %Find the indices which best match these values
    [~,paraVertInd{i}] = arrayfun(@(x)(...
        min(abs(distAlong{i}-distVals{i}(x)))),1:nWinPara(i)+1);

    %Get the coordinates of these indices
    paraVert{i} = contours{iPerpStart(i)}(:,paraVertInd{i});

    %Find start points on image boundary - these are a special case
    sOnBoundA = arrayfun(@(x)(any(ceil(contours{iPerpStart(i)}... %On high boundary
                (:,paraVertInd{i}(x)))==[N,M]')),1:nWinPara(i)+1);
    sOnBoundB = arrayfun(@(x)(any(floor(contours{iPerpStart(i)}...%On low boundary
                (:,paraVertInd{i}(x)))==[1;1])),1:nWinPara(i)+1);    

    % ---- Slice it up! ---- %
    %Compliments of Sylvain! 
    slices{i} = cell(1,nWinPara(i)+1);
    slices{i}(~(sOnBoundA | sOnBoundB)) = gradientDescent(double(max(distX(:)) - distX),... %invert to get gradient ascent
            contours{iPerpStart(i)}(1,paraVertInd{i}(~(sOnBoundA | sOnBoundB))),...
            contours{iPerpStart(i)}(2,paraVertInd{i}(~(sOnBoundA | sOnBoundB))));

    %Start points on the image boundary will terminate immediately during
    %gradient ascent/descent. We just shift these inwards to avoid this.
    if any(sOnBoundA)
        tmpStartsA = contours{iPerpStart(i)}(:,paraVertInd{i}(sOnBoundA))-gaShift;        
        slices{i}(sOnBoundA) = gradientDescent(double(max(distX(:)) - distX),... 
                                         tmpStartsA(1,:),tmpStartsA(2,:));
    end
    if any(sOnBoundB)
        tmpStartsB = contours{iPerpStart(i)}(:,paraVertInd{i}(sOnBoundB))+gaShift;                    
        slices{i}(sOnBoundB) = gradientDescent(double(max(distX(:)) - distX),... 
                                         tmpStartsB(1,:),tmpStartsB(2,:));            
    end


    %Transpose these slices so they match with the contours
    slices{i} = cellfun(@(x)(x'),slices{i},'UniformOutput',false);

    %If the starting contour was not the object boundary, we need to extend
    %the slices to the boundary by doing a gradient descent from the same
    %points
    if contourValues(iPerpStart(i)) > 0

        tmp = cell(1,nWinPara(i)+1);
        tmp(~(sOnBoundA | sOnBoundB)) = gradientDescent(double(distX),...
        contours{iPerpStart(i)}(1,paraVertInd{i}(~(sOnBoundA | sOnBoundB))),...
        contours{iPerpStart(i)}(2,paraVertInd{i}(~(sOnBoundA | sOnBoundB))));    
        if any(sOnBoundA)            
            tmp(sOnBoundA) = gradientDescent(double(distX),... 
                                             tmpStartsA(1,:),tmpStartsA(2,:));
        end
        if any(sOnBoundB)            
            tmp(sOnBoundB) = gradientDescent(double(distX),... 
                                             tmpStartsB(1,:),tmpStartsB(2,:));
        end        
        %Transpose and reverse these...
        tmp = cellfun(@(x)(x(end:-1:1,:)'),tmp,'UniformOutput',false);
        %... and add them to the slices
        slices{i} = arrayfun(@(x)([tmp{x}(:,1:end-1) slices{i}{x}]),1:length(tmp),'UniformOutput',false);                
    end
end

%If the start contour was > 1, and this gave nStart>1, we need to
%combine and order the slices from the different starting contours.
if nStart > 1
    iZeroContInt = cell(1,nStart);
    jZeroContInt = cell(1,nStart);
    for i = 1:nStart
        
        %Determine the indices where these slices intersect the zero-value
        %contour.
        [~,~,iZeroContInt{i},jZeroContInt{i}] = find_intersections(contours{iZeroCont},...
                                                slices{i});        
        
    end   
    
    %Combine and sort the intersection indices, retaining only unique
    %intersectiosn
    [iZeroContInt,sortZeroInt] = unique(cat(1,iZeroContInt{:}));

    
    %The sorting and uniqueness may remove the duplicated slice at the
    %contour start. Replace this.
    if isClosed(iZeroCont) && iZeroContInt(1) ~= iZeroContInt(end)            
        sortZeroInt(end+1) = sortZeroInt(1);
    end    
    
    %Sort the slices using this ordering
    slices = cat(2,slices{:});
    slices = slices(sortZeroInt);
    %Sort the paraVertInd also
    paraVertInd = [paraVertInd{:}];    
    paraVertInd = paraVertInd(sortZeroInt);
    paraVert = [paraVert{:}];
    paraVert = paraVert(:,sortZeroInt);
            
else
    %Un-cell the 1x1 cell arrays
    slices = slices{1};
    paraVertInd = paraVertInd{1};
    paraVert = paraVert{1};
end

if spOnly
    windows = paraVert';
    return
end

%% ------- Windowing ------- %%
%Finds intersections of slices and contours (perpindicular and parallel
%divisions respectively) to determine window polygons


%Since each window depends on two slices, we do the first slice as a special case.

%Find the intersections of the first slice.
% [intXprev,intYprev,iSintPrev,iCintPrev] = find_intersections(slices{1},...                                                
%                                                    contours);                                                                                                                         
% 
% iContIntPrev = find(~isnan(intXprev));
% nBandPrev = numel(iContIntPrev)-1;%One band less than the number of contour intersections
% 

nStrips = numel(slices)-1;%we have one window strip less than the number of slices
windows = cell(1,nStrips);

%Go through each subsequent slice, and find intersections w/ contours.    
for j = 1:(nStrips+1)

    [intXcur,intYcur,iSintCur,iCintCur] = find_intersections(slices{j},...
                                                    contours);

    iContIntCur = find(~isnan(intXcur));

    %If the starting contour was zero, the contour and the slice just
    %barely touch and the intersection algorithm (sometimes) doesn't
    %count this as an intersection. We need to add these back... 
    if ~any(iContIntCur == iZeroCont)        
        if startContour == 1        
            iCintCur(iZeroCont) = paraVertInd(j);            
        else
            %Since there cannot be true local minima within the distance
            %transform, if the slice terminates before reaching the
            %zero-value contour, it either hit the image edge or numerical
            %error caused it to stop just short. We therefore introduce an
            %intersection if the termination was close enough
            
            %Find the point on the zero contour which is closes to the end
            %of the slice
            [minDist,iMinDist] = min(sqrt((slices{j}(1,1) - contours{iZeroCont}(1,:)) .^2 + ...
                                          (slices{j}(2,1) - contours{iZeroCont}(2,:)) .^2));
            %If it's close enough, use this point as the intersection
            if minDist < perpSize                
                iCintCur(iZeroCont) = iMinDist;                                
            end
        end
        iSintCur(iZeroCont) = 1;
        intXcur(iZeroCont) = contours{iZeroCont}(1,iCintCur(iZeroCont));
        intYcur(iZeroCont) = contours{iZeroCont}(2,iCintCur(iZeroCont));
        iContIntCur = find(~isnan(intXcur));
    end
    nBandCur = numel(iContIntCur)-1;                

    %First slice is a special case - a window can only be formed with two
    %slices.
    if j > 1
    
        %The number of bands here is determined by the minimum number of
        %intersections of this slice and the previous slice, since both are
        %needed to create a window.
        windows{j-1} = cell(1,min(nBandCur,nBandPrev));

        for k = 1:nBandPrev

            if k <= nBandCur %This is pointless for now, but will be needed later if I decide to add back windows which meet image boundary.

                %Make sure that both slices intersect the same contours.
                if all(iContIntPrev(k:k+1) == iContIntCur(k:k+1))                                                

                    %Corners run clockwise, starting with lower left
                    a = [intXprev(iContIntPrev(k)) ; intYprev(iContIntPrev(k))];              %Corner A
                    b = [intXprev(iContIntPrev(k+1)) ; intYprev(iContIntPrev(k+1))];          %Corner B
                    c = [intXcur(iContIntCur(k+1)) ; intYcur(iContIntCur(k+1))];              %Corner C
                    d = [intXcur(iContIntCur(k)) ; intYcur(iContIntCur(k))];                  %Corner D

                    %Make sure that the window has not completely collapsed into a
                    %line.
                    if any([sqrt(sum((a-d) .^2)) sqrt(sum((b-c) .^2))] > collapsedSize)                                               

                        ab = slices{j-1}(:,ceil(iSintPrev(iContIntPrev(k))):...
                            floor(iSintPrev(iContIntPrev(k+1))));                                                     %Side A->B


                        if iCintPrev(iContIntPrev(k+1)) <= iCintCur(iContIntCur(k+1)) + intErr                             %Side B->C
                            %The "normal" case, where the previous intersection
                            %has an index lower than the current, or where the
                            %window has collapsed and the intersections are the
                            %same.
                            bc = contours{iContIntPrev(k+1)}(:,ceil(iCintPrev(iContIntPrev(k+1))):...
                                    floor(iCintCur(iContIntCur(k+1)))); 
                        elseif isClosed(iContIntPrev(k+1))

                            %In the case that this window side crosses the
                            %starting point of this contour...
                            bc = contours{iContIntPrev(k+1)}(:,...
                                [ceil(iCintPrev(iContIntPrev(k+1))):end 1:floor(iCintCur(iContIntCur(k+1)))]);                         
                        else
                            %If we get here, something has gone wrong!
                            if doChecks
                                error('What the fuck? this should never happen - report this to Hunter so he can check windowing algorithm!!!')
                            else
                                error('Problem with mask - check mask yourself or enable the mask checking option!')
                            end
                        end

                        cd = slices{j}(:,floor(iSintCur(iContIntCur(k+1))):-1:ceil(iSintCur(iContIntCur(k))));         %Side C->D

                        if iCintCur(iContIntCur(k)) >= iCintPrev(iContIntPrev(k)) - intErr                                             %Side D->A
                            %The "normal" case - this side runs anti-parallel
                            %to the contour so we expect the intersections to
                            %occur in decreasing order.
                            da = contours{iContIntCur(k)}(:,floor(iCintCur(iContIntCur(k))):-1:ceil(iCintPrev(iContIntPrev(k))));
                        elseif isClosed(iContIntCur(k))
                            %This is the case where this window side
                            %crosses the start point of this contour.
                            da = contours{iContIntCur(k)}(:,[floor(iCintCur(iContIntCur(k))):-1:1 end:-1:ceil(iCintPrev(iContIntPrev(k)))]);
                        else %If we get here, something has gone wrong!
                            if doChecks
                                %error('What the fuck? this should never happen - report this to Hunter so he can check windowing algorithm!!!')
                            else
                                error('Problem with mask - check mask yourself or enable the mask checking option!')
                            end
                        end

                        %Combine all the vertices and sides to create the
                        %window
                        windows{j-1}{k}  = {[a ab],[b bc],[c cd],[d da]};                                                                                                                                                                

                        isCollapsed = false;

                        if showPlots

                            if firstTime
                                firstTime = false;
                                fsFigure(.6);
                                hold on                            
                                gX = gradient(distX);%Show gradient of distance transform so ridgelines are obvious
                                imagesc(gX);
                                colormap gray
                                plotDirection(contours,'g');  
                                plotDirection(slices,'g');
                                axis image,axis ij,set(gca,'color','none')
                            end
                            winBorder = [windows{j-1}{k}{:}];
                            fill(winBorder(1,:),winBorder(2,:),'y','FaceAlpha',.5)
                            plot(ab(1,:),ab(2,:),'r.');
                            plot(bc(1,:),bc(2,:),'g.');
                            plot(cd(1,:),cd(2,:),'b.');
                            plot(da(1,:),da(2,:),'m.');
                        end

                    else%If the window has collapsed, we need to terminate this strip of windows.                        
                        windows{j-1} = windows{j-1}(1:k-1);                    
                        isCollapsed = true;
                        break
                    end
                 end


            end

        end       

        %If the two slices end in the same local maxima of the distance
        %transform, and the strip hasn't collapsed, we need to add a final
        %3-corner window to cover the areas inside the innermost contours
        if nBandCur > 0 && nBandCur == nBandPrev && ...
                all(round(slices{j-1}(:,end)) == ...
                    round(slices{j}(:,end))) && ~isCollapsed

            %Corner where slices hit local maxima of distance transform
            e = slices{j-1}(:,end); %Should we instead use the first point where they overlap?

            %Side from previous windows corner b to intersection e
            be = slices{j-1}(:,ceil(iSintPrev(iContIntPrev(nBandCur+1))):end);
            %Side from intersection back to previous windows c corner               
            ec = slices{j}(:,end:-1:ceil(iSintCur(iContIntCur(nBandCur+1))));  
            %Previous windows cb side               
            cb = contours{iContIntPrev(nBandCur+1)}(:,floor(iCintCur(iContIntCur(nBandCur+1))):-1:ceil(iCintPrev(iContIntPrev(nBandCur+1))));    

            windows{j-1}{nBandCur+1} = {[b be],[e ec],[c cb]};

            if showPlots

                    winPoly = [windows{j-1}{nBandCur+1}{:}];

                    fill(winPoly(1,:),winPoly(2,:),'m','FaceAlpha',.5)

                    plot(be(1,:),be(2,:),'r.');
                    plot(ec(1,:),ec(2,:),'g.');
                    plot(cb(1,:),cb(2,:),'b.');                   
            end                          

        end
    end
    intXprev = intXcur;
    intYprev = intYcur;
    iSintPrev = iSintCur;
    iCintPrev = iCintCur;
    nBandPrev = nBandCur;
    iContIntPrev = iContIntCur;
end

    


function [ix,iy,i1,i2] = find_intersections(c,cInt)
%This finds intersections, ensuring that only one intersection is returned
%per contour.
%See below for details of why this is necessary.

    nCon = length(cInt);
        
    ix = nan(nCon,1);
    iy = nan(nCon,1);
    i1 = nan(nCon,1);
    i2 = nan(nCon,1);
    
    
    for j = 1:nCon
        
        %This can be sped up by taking into account that each subsequent
        %intersection will occur further along the slice... HLE
        
            
        [tmpix,tmpiy,tmpi1,tmpi2] = intersectionsHLE( ...
                                              c(1,:),c(2,:), ...
                                              cInt{j}(1,:),cInt{j}(2,:));
        
        if ~isempty(tmpi1)                                    
                        
            if numel(tmpi1) > 1
            
                %Deal with duplicate intersections. In general the slices
                %and contours should only intersect at one point (because
                %they are perpindicular on the surface of the distance
                %transform), but some rare situations cause exceptions to
                %this rule.
                
                if numel(unique(tmpi1)) > 1
                   
                    %This can only be caused (as far as I know) by the
                    %contour running along a rideline and the gradient
                    %ascent just slightly overshooting this ridgeline, and
                    %then intersecting the contour again (gradient ascent
                    %is sub-pixel, while contouring is at the pixel level)
                    %Therefore we need only keep the first intersection(s).
                    [~,i1Unique] = unique(tmpi1);                                                            
                    iKeep = tmpi1 == tmpi1(i1Unique(1));
                    tmpi1 = tmpi1(iKeep);
                    tmpi2 = tmpi2(iKeep);
                    tmpix = tmpix(iKeep);
                    tmpiy = tmpiy(iKeep);
                    
                end
                
                %Check again in case the above check removed the duplicate
                %intersection
                if numel(tmpi1) == 2                                                                                                                                 

                    %This can only be caused (again, as far as I know) by
                    %the contour being colinear with itself because it is
                    %running along a saddle. This will give two
                    %intersections at the exact same x and y and index on
                    %the slice, but with different index on the contour.

                    %We want to take the intersection which corresponds to
                    %the part of the contour that is running in the same
                    %direction as the slices (clockwise). So we check the
                    %directions via the cross-product and use that to
                    %choose the intersection to return.
                    
                    %--Set up the vectors for taking the cross-product--%
                   
                    %Get the vector going in the direction of the slice at
                    %the intersection.
                    %tmpi1(1) and tmpi1(2) are identical due to the
                    %check above, so just use 1
                    if tmpi1(1) == round(tmpi1(1))
                        %If the intersection occurs right at a vertex...
                        if tmpi1(1) == size(c,2)
                            %... and this is the last vertex, use the preceding segment
                            vecSlice = [diff(c(1,end-1:end)) ...
                                        diff(c(2,end-1:end)) ...
                                        0];                                                    
                        else
                            %...otherwise, use the next segment
                            vecSlice = [diff(c(1,tmpi1(1):tmpi1(1)+1)) ...
                                        diff(c(2,tmpi1(1):tmpi1(1)+1))...
                                        0];                                                                                
                        end                    
                    else
                        %... or the vector from the intersection to the next vertex
                        vecSlice = [c(1,ceil(tmpi1(1))) - tmpix(1) ...     
                                    c(2,ceil(tmpi1(1))) - tmpiy(1) ...     
                                    0];                                    
                    end
                    %Get the vector in the direction of the contour at its
                    %first intersection
                    if tmpi2(1) == round(tmpi2(1))                        
                        if tmpi2(1) == size(cInt{j},2)                            
                            vecCont1 = [diff(cInt{j}(1,end-1:end)) ...
                                        diff(cInt{j}(2,end-1:end)) ...
                                        0];                                                    
                        else                            
                            vecCont1= [diff(cInt{j}(1,tmpi2(1):tmpi2(1)+1)) ...
                                       diff(cInt{j}(2,tmpi2(1):tmpi2(1)+1))...
                                        0];                                               
                        end    
                    else
                        %... or the vector from the intersection to the next vertex
                        vecCont1 = [cInt{j}(1,ceil(tmpi2(1))) - tmpix(1) ...     
                                    cInt{j}(2,ceil(tmpi2(1))) - tmpiy(1) ...     
                                    0];                                    
                    end
                    %Get the vector in the direction of the contour at its
                    %second intersection
                    if tmpi2(2) == round(tmpi2(2))                        
                        if tmpi2(2) == size(cInt{j},2)                            
                            vecCont2 = [diff(cInt{j}(1,end-1:end)) ...
                                        diff(cInt{j}(2,end-1:end)) ...
                                        0];                                                    
                        else                            
                            vecCont2= [diff(cInt{j}(1,tmpi2(2):tmpi2(2)+1)) ...
                                       diff(cInt{j}(2,tmpi2(2):tmpi2(2)+1))...
                                        0];                                               
                        end    
                    else
                        %... or the vector from the intersection to the next vertex
                        vecCont2 = [cInt{j}(1,ceil(tmpi2(2))) - tmpix(1) ...     
                                    cInt{j}(2,ceil(tmpi2(2))) - tmpiy(1) ...     
                                    0];                                    
                    end                                        
                    
                    %Take the cross product, and use the intersection where
                    %the cross has a positive z-component - this is the one
                    %that is running clockwise.
                    crossP1 = cross(vecSlice,vecCont1);                    
                    crossP2 = cross(vecSlice,vecCont2);
                    
                    if crossP1(3) > 0
                        iKeep = 1;
                    elseif crossP2(3) > 0
                        iKeep = 2;
                    elseif all(crossP1==0) && all(crossP2 == 0)
                        %If the slice is parallel/antiparallel with the
                        %contour here, we take the parallel one
                        if dot(vecSlice,vecCont1) > 0
                            iKeep = 1;
                        else
                            iKeep = 2;
                        end                                                
                    else                        
                        error('Problem with multiple contour-slice intersections: Ambiguous cross product!')
                    end                                     
                                        
                    i1(j) = tmpi1(iKeep);
                    i2(j) = tmpi2(iKeep);
                    ix(j) = tmpix(iKeep);
                    iy(j) = tmpiy(iKeep);
                        
                elseif numel(tmpi1) == 1
                    
                    ix(j) = tmpix;
                    iy(j) = tmpiy;
                    i1(j) = tmpi1;
                    i2(j) = tmpi2;
                    
                elseif numel(tmpi1) > 2                    
                    error('Problem creating windows: too many slice-contour intersections!')                    
                end
                    
            else                                
                ix(j) = tmpix;
                iy(j) = tmpiy;
                i1(j) = tmpi1;
                i2(j) = tmpi2;
            end
        end
        
    end                                


function [startPoint,nPara,startContour,showPlots,doChecks,spOnly] = parseInput(argArray)

startPoint = [];
nPara = [];
startContour = [];
showPlots = [];
doChecks = [];
spOnly = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg
    
    
    switch argArray{i}
        
        
        case 'StartPoint'
            startPoint = argArray{i+1};
                        
        case 'NumParallel'
            nPara = argArray{i+1};

        case 'StartContour'
            startContour = argArray{i+1};
               
        case 'ShowPlots'
            
            showPlots = argArray{i+1};
            
        case 'DoChecks'
            
            doChecks = argArray{i+1};
            
        case 'StartPointsOnly'
            
            spOnly = argArray{i+1};
            
        otherwise

            error(['"' argArray{i} '" is not a valid option name! Please check input!'])
    end
end                       


