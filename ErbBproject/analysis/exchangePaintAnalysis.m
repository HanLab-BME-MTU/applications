function [ fullpath ] = exchangePaintAnalysis( dir, varargin )
%exchangePaintAnalysis
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

ip.addOptional('display',0,@isnumeric);
ip.addOptional('pixelSize',62.81,@isnumeric);
ip.addOptional('DoContinuous',true,@islogical);

ip.parse(dir,varargin{:});

display=ip.Results.display;
pixelSize = ip.Results.pixelSize;
doContin = ip.Results.DoContinuous;
    
list = findFilesInSubDirs(dir,'*tracking.mat');

%The following takes the list of list extracts some general diagnostic
%information for each file (average number of tracks per frame)

GeneralDiagnostic = cell(size(list));
index =false(size(list));
for i=1:numel(list)
    load(list{i},'tracksFinal','features');
    list(i)   
    
    tracksPerFrame = numel(tracksFinal)/numel(features)
    
    l=0;
    for j=1:numel(features)
       if ~isempty(features{j})
        l = l+numel(features{j}.x); 
       end
    end
    
    localizationsPerFrame = l/numel(features)
    
    %If this movie is a control movie it will have 'pbs' or 'pre' in its name
    slash = strfind(list{i},filesep);
    name = list{i}(slash(end)+1:end);
    index(i)=isempty(strfind(name,'pbs')) & isempty(strfind(name,'pre'));
    
    GeneralDiagnostic(i)={struct('name',list{i},'tracksPerFrame',tracksPerFrame,'localizationsPerFrame',localizationsPerFrame)};
end

%removes control files from the list
list = list(index);

f = strfind(dir,'/');

name = [dir(f(end-1)+1:f(end)-1),'_',dir(f(end)+1:end)];

fullpath = exchangePaintAlignment(list,name,'ImageDisp',display);

load(fullpath,'PointList');

%Used in several places
n = numel(PointList);

%Make continous density estmates using pairwise L function

if doContin

    ContinuousAnalysis = cell(n);
    Lr = cell([n,1]);
    
    %This loop calculates all the pairwise L(r) functions then modifies it to
    %the H(r) function H(r) = L(r) - r
    %You have to do all as they are not symetric
    
    %these are the radii to use for L(r) cross statistics
    r = 0.2:0.2:5.6;
    
    for i=1:n
        %calculate Lr (self not cross)
        Lr{i} = PointP_Lr_total(PointList{i}.com,r);
        for j =1:n
            %Computes Besag's L and then renormalizes it to the H statistic
            tmp = PointP_Lr_Cross(PointList{i}.com,PointList{j}.com,r);
            tmp2 = tmp./r';
            % This is only a crude way to find the right normalization
            Const = mean(tmp2(15:end));
            ContinuousAnalysis{i,j} = (tmp/Const)-r';
            
        end
    end

end


%Calculate clustering using mean shift
TotalClust = [];
TotalPnts=[];

%Combine all track COMs, while keeping the identity of which channel
%they came from.
n = numel(PointList);
for i = 1:n
    TotalPnts = vertcat(TotalPnts,[PointList{i}.com,i*ones(size(PointList{i}.com))]);
end

TotalPnts=TotalPnts(~isnan(TotalPnts(:,1)),:);

[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.5,'kernel','flat');

for i=1:numel(clusterInfo)

    pnts = TotalPnts(clusterInfo(i).ptIdData,:);

    if numel(pnts(:,1)) > 2
    %finds convex hull and area
    [hull,area]= convhull(pnts(:,1:2));
    else
        hull =[];
        area = NaN;
    end
    
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
    save(fullpath, 'clusterInfo','clusterMap','GeneralDiagnostic','ContinuousAnalysis','Lr','r','-append');
    
      
    
end


