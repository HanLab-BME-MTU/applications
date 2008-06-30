function polyDepolyChooseRoi(runInfo)

if nargin<1
    error('polyDepolyChooseRoi: Not enough input parameters')
end
if ~isstruct(runInfo)
    runInfo=struct;
end

if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('polyDepolyChooseRoi: runInfo should contain fields imDir and anDir');
else
    [anDir] = formatPath(runInfo.anDir);
    [imDir] = formatPath(runInfo.imDir);
end

% norm directory where we will store roiGeom
negEdgeDir=[imDir filesep 'norm' filesep 'negEdgeTifs'];
if isdir(negEdgeDir)
    [listOfImages] = searchFiles('neg_edge',[],negEdgeDir,0);
    img=imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))]);
    [polyDepolyROI, xi, yi] = roipoly(img);
    imwrite(polyDepolyROI,[anDir filesep 'polyDepolyROI.tif']);
    close
else
    error('polyDepolyChooseRoi: path messed up')
end

