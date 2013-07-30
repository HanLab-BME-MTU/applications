function finished = exchangePaintAlignment(list,name,varargin)
% Exchange Paint Alignment
%
% This code takes analyzed and tracked data from a multicolor paint movie
% It drift corrects ands aligns all the localizations from each movie
% through drift markers.
%
%
% All the .mat files 
% should contain tracksFinal, features and the associated MD
%
%Input: 
%       list: list of full path of .mat files to be included is a cell
%              array
%       
%       name: name to save the resulting analysis under
%
%Optional: 
%       pixelSize, defines the size in nanometers of a pixel
%
%       ImageDisp, if set to true displays the result of the alignment
%
%       dir, the directory to save the output file in. The default is the
%            current directory
%
%Output:
%       finished, returns 1 upon completion
%
%       (PointList), a Cell array of structures that is saved to the disk as name.mat
%                  .pnts, list of pnts drift and shift corrected
%                  .drift, a vector with the x,y shift from time zero at
%                          each frame
%                  .dmark, x,y corrodinates of the drift markers
%                  .name, name of the analyzed file
%                  .shift, x,y shift relative to the first in the series
%
%
% 2013/07/15 Jeffrey Werbin
% Harvard Medical School
%

finished = 1;

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('dir',@iscell);
ip.addRequired('name',@ischar);


ip.addOptional('pixelSize',62.81,@isscalar);
ip.addOptional('ImageDisp',false,@islogical);
ip.addOptional('dir',cd(),@ischar);

ip.parse(dataDirectory,varargin{:});

pixelSize = ip.Results.pixelSize;
ImageDisp = ip.Results.ImageDisp;
dir = ip.Results.dir;

%Get a list of the .mat files in current directory
%We assume that all files
%L = what;
%list = L.mat;



maxGap =10;

MinTrackLen = 2;

num = numel(list);

PointList = cell([num,1]);
movieInfo = struct('pnts',[],'drift',[],'dmark',[],'name',[],'shift',[]);

for j=1:num
    PointList{j} = movieInfo;
    PointList{j}.name = list{j};
    load(list{j});
    
    track = reformTracksFinal(tracksFinal,500);
        
    %Calculates average drift from all the drift markers
    ind = find(vertcat(track.isDrift));
    %MovieLength = MD.nFrames_;
    %MovieLength = 5000; %temporary measure
    MovieLength = numel(features);
    drift = zeros([MovieLength,2]);


    if ~isempty(ind)
        %Keeps track of the number of markers that appear in this frame
        numpnts = drift;

        for i = 1:numel(ind)
            x = track(ind(i)).coord(:,1);
            y = track(ind(i)).coord(:,2);

            
            if sum(isnan(x)) > 0
                x = gapInterpolation(x,maxGap)';
                y = gapInterpolation(y,maxGap)';
            end            

            t = track(ind(i)).timeInfo(1:2); %time range in frames
            drift(t(1)+1:t(2),:) = drift(t(1)+1:t(2),:) + [diff(x),diff(y)];
            numpnts(t(1)+1:t(2),:) = numpnts(t(1)+1:t(2),:) + ones(size(numpnts(t(1)+1:t(2),:)));
        end

        %makes sure that the first point is not neglected
        numpnts(1,:)=1;

        % finds points with no data to avoid divide by zero
        DriftMiss = find(numpnts(:,1) == 0);
        numpnts(DriftMiss,:)=1;

        if ~isempty(DriftMiss)
           ['We have a few (',num2str(numel(DriftMiss)),') missed frames of correction.']
        end

        drift = drift./numpnts;
        drift = cumsum(drift);
    else
        error('file %s is missing drift markers', PointList.name);
    end

    %Here drift correction is applied to all points
      
    for i =1:numel(track)
    frames = track(i).amp(:,end); % the frame in which each point appear
    coord =track(i).coord;
    coord(:,[1,2]) = coord(:,[1,2]) - drift(frames,:);
    track(i).coord = coord;
    
    % Recaculate center of mass calculated as a weigthed mean
    idx=~isnan(coord(:,1));
    tmp = coord(idx,:);
    [nrows,ncols]=size(tmp);
    if( nrows > 1 )
        [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,3:4),'s');
        track(i).com=[wm,ws,sqrt(sum(ws.^2))];
    else
        track(i).com=[tmp(1:2), tmp(3:4), sqrt(sum(tmp(3:4).^2))];
    end
    
    end
        
    PointList{j}.pnts = vertcat(track(vertcat(track.num) > MinTrackLen).coord);
    PointList{j}.drift = drift;
    
    %makes a list of drift markers intial points
    dmark = zeros([numel(ind),2]);
    
    for i = 1:numel(ind)
        dmark(i,:)=track(ind(i)).coord(1,1:2);
    end
    
    %Stores drift mark intial positions after drift correcting
    PointList{j}.dmark = dmark;
end    

% Find the shift between movies using dmark
% assumes x,y shifts only (no rotations or dialations)
%
% uses the first movie as a reference

ref = PointList{1}.dmark;
s_ref = size(ref);
PointList{1}.shift = struct('regParam',[],'Bfit',[],'ErrorStats',[],'trans',[0,0],'clusterInfo',[],'clusterMap',[])

for j = 2:num
    test = PointList{j}.dmark;
    s_test = size(test);
    dm = distMat2(ref,test);
    
    % Assigns drift markers in ref to markers in test
    [OneToTwo,TwoToOne]=lap(dm,-1,0,1);
    
    %Finds the transform that reconciles the two sets of points
    % is dependent on absor code
    tmp = TwoToOne(1:s_test(1));
    
    if s_ref(1) > s_test(1)
       tmp = OneToTwo(1:s_ref(1));
       A = ref(tmp(tmp <= s_test(1)),:);
       B = test;
       [regParam,Bfit,ErrorStats]=absor(A',B');
       shift = mean(A-B)
    else
        tmp = TwoToOne(1:s_test(1));
        A= ref;
        B=test(tmp(tmp <= s_ref(1)),:);
       [regParam,Bfit,ErrorStats]=absor(A',B');
       shift = mean(A-B)
    end
    
    PointList{j}.shift = struct('regParam',regParam,'Bfit',Bfit,'ErrorStats',ErrorStats,'trans',shift);
    
    %Applies Transform to the data
    
    %just x,y translation for now
    tmp = PointList{j}.pnts(:,1:2)+repmat(shift,[numel(PointList{j}.pnts(:,1)),1]);
    PointList{j}.pnts=tmp(~isnan(tmp(:,1)),:);
    
    %Applies mean shift tracking to one type of receptor/protein
    
    [clusterInfo,clusterMap]=MeanShiftClustering(PointList{j}.pnts(:,1:2),0.5,'kernel','flat');
    
    PointList{j}.clusterInfo = clusterInfo;
    PointList{j}.clusterMap = clusterMap;

end

tmp = PointList{1}.pnts(:,1:2);
PointList{1}.pnts=tmp(~isnan(tmp(:,1)),:);
[clusterInfo,clusterMap]=MeanShiftClustering(PointList{1}.pnts(:,1:2),0.5,'kernel','flat');
PointList{1}.clusterInfo = clusterInfo;
PointList{1}.clusterMap = clusterMap;

%save results

save([dir,pathsep,name,'.mat'],'PointList','maxGap','MinTrackLen');


%Create and "Image" of the final merge

cmap = [{'g'},{'r'},{'b'},{'k'},{'m'}];

if ImageDisp

    %first just the aligned tracking markers pre and post
    figure;
    hold;
    for j=1:num
        tmp=PointList{j}.dmark;
        ntmp = numel(tmp(:,1));
        scatter(tmp(:,1),tmp(:,2),[cmap{j},'x']);
        scatter(tmp(:,1)+repmat(PointList{j}.shift.trans(1),[ntmp,1]),tmp(:,2)+repmat(PointList{j}.shift.trans(2),[ntmp,1]),[cmap{j},'s']);
    end
    title('Drift Markers only')

    figure;
    hold;
    for j=1:num
        tmp=PointList{j}.pnts;
        scatter(tmp(:,1),tmp(:,2),[cmap{j},'.']);
    end
    title('5 receptor image') 
    legend('EGFR','ErbB2','ErbB3','IGF1R','Met')

end


end