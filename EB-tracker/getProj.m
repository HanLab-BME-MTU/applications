function [runInfo]=getProj(projList,targetName,oligoName,movieName,roiName)
% GETPROJ returns paths to projects matching user query

% INPUT:
% projList: nSubprojects x 7 (or more) cell array of the form:
%           [target oligo movie roi roi_imageDir roi_analysisDir origImDir optional ...]
%           'projList' is derived from the setupDirectories function. It
%           assumes that the data is in a directory hierarchy like this:
%           target/oligo/movie/rois/roi_n, where roi(s) are made by copying
%           the original images or querying the user to select a polygonal
%           region of interest. The optional fields can contain data
%           information or pointers to the results once the projects have
%           been analyzed.
%
% targetName: search string for molecular target (or fourth-lowest
%             directory on the tree - see above); must match exactly; case sensitive
% 
% oligoName : search string for siRNA oligo (or the third-lowest
%             directory on the tree); must match exactly; case sensitive
%
% movieName : search string for particular movie (or second-lowest
%             directory on the tree); must match exactly; case sensitive
%
% roiName   : search string for particular roi_n (or lowest directory on the
%             tree); must match exactly; case sensitive
%
% OUTPUT:
% runInfo: structure containing roi_imageDirectory and
%          roi_analysisDirectory for each sub-project fitting the query



if nargin<1 || isempty(projList)
    error('getProj: not enough input arguments');
end

if nargin<2 || isempty(targetName)
    targetNameMatch=ones(size(projList,1),1); % consider all
else
    targetNameMatch=cellfun(@(x) isequal(x,targetName), projList(:,1));
end

if nargin<3 || isempty(oligoName)
    oligoNameMatch=ones(size(projList,1),1);
else
    oligoNameMatch=cellfun(@(x) isequal(x,oligoName), projList(:,2));
end

if nargin<4 || isempty(movieName)
    movieNameMatch=ones(size(projList,1),1);
else
    movieNameMatch=cellfun(@(x) isequal(x,movieName), projList(:,3));
end

if nargin<5 || isempty(roiName)
    roiNameMatch=ones(size(projList,1),1);
else
    roiNameMatch=cellfun(@(x) isequal(x,roiName), projList(:,4));
end

matches=find(targetNameMatch & oligoNameMatch & movieNameMatch & roiNameMatch);

runInfo=[];
for i=1:length(matches)
    runInfo(i,1).imDir=projList{matches(i),5};
    runInfo(i,1).anDir=projList{matches(i),6};
end


