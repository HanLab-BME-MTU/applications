function [ fullpath ] = exchangePaintAnalysis( dir, varargin )
%exchangePaintAnalysisWrapper
%   Takes a base directory and creates a list of .mat files that are
%   associate with a single exchangePaint image, then aligns them and
%   preforms some data analysis
%   
%
% Inputs:
%           dir, base directory for an image should be named "ImageX" where
%                X is a number
%
%Optional:
%
%
%Output: 
%        fullpath: the path to a .mat files containing the results
%
%
%Written by Jeffrey Werbin 
%Harvard Medical School 2013/08/02
%

ip = inputParser;

ip.addRequired('dir',@ischar);

ip.parse(dir,varargin{:});

    
list = findFilesInSubDirs(dir,'*tracking.mat');

%The following takes the list of list extracts some general diagnostic
%information for each file (average number of tracks per frame)

GeneralDiagnostic = cell(size(list));
index =false(size(list));
for i=1:numel(list)
    load(list{i},'tracksFinal','features');
    list(i)   
    trackPerFrame = numel(tracksFinal)/numel(features)
    
    %If this movie is a control movie it will have 'pbs' in its name
    index(i)=isempty(strfind(list{i},'*pbs*'));
    
    GeneralDiagnostic(i)={struct('name',list{i},'tracksPerFrame',tracksPerFrame)};
end

%removes control files from the list
list = list(index);

f = strfind(dir,'/');

name = [dir(f(end-1)+1:f(end)-1),'_',dir(f(end)+1:end)];

fullpath = exchangePaintAlignment(list,name);

load(fullpath,'PointList');

%Make continous density estmates using pairwise L function



%Calculate clustering using mean shift
TotalClust = [];
TotalPnts=[];

%Combine all localizations, while keeping the identity of which channel
%they came from.
n = numel(PointList);
for i = 1:n
    TotalPnts = vertcat(TotalPnts,[PointList{i}.pnts,i*ones(size(PointList{i}.pnts))]);
end


[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.5,'kernel','flat');

for i=1:numel(clusterInfo)

    pnts = TotalPnts(clusterInfo(i).ptIdData,:);

    %finds convex hull and area
    [hull,area]= convhull(pnts(:,1:2));

    % n is the number of pnts that are merged. Here we store the
    % number of points in this merged cluster from each element of the
    % shift array

    composition = hist(pnts(:,3),1:n);

    
    %Add values to clusterInfo
    clusterInfo(i).pnts = pnts;
    clusterInfo(i).area = area*pixelSize^2;
    clusterInfo(i).hull = hull;
    clusterInfo(i).composition = composition;
    
end

%Appends the additional analysis to fullpath file
    save(fullpath, 'clusterInfo','clusterMap','-append');
    
      
    
end


