function allDistances=fsmEstSearchRadius(candsDir,maxRadius)
% fsmEstSearchRadius uses cands structures to estimate the search radius for tracking
%
% fsmEstSearchRadius calculates and stores all the MIN distances between all pairs of particle positions in consecutive
%    frames over the whole range of analyzed images. This does not match exactly the distribution of actual track results, 
%    but should give a good estimate of the expected mean flow speed.
%
% SYNOPSIS      allDistances=fsmEstSearchRadius(candsDir,maxRadius)
%
% INPUT         candsDir     : complete path to the location of cands###.mat files
%               maxRadius    : maximum radius in pixels to be considered (select a value which does
%                              not significantly cut the distribution of distances)
%
% OUTPUT        allDistances : vectors containing all MIN distance between pairs of particles in consecutive frames over the whole
%                              range of analyzed frames
%
% DEPENDENCES   fsmEstSearchRadius uses { createSparseDistanceMatrix }
%               fsmEstSearchRadius is used by { fsmPostProc }
%
% Aaron Ponti, September 10th, 2004

% Current directory
oldDir=cd;

if ~isdir(candsDir)
    uwait(errordlg('Invalid ''cands'' directory.'));
    return
end

% Read directory
candsList=dir(candsDir);

% Remove directories
candsList(find([candsList.isdir]==1))=[];

% Check that there are files left
if length(candsList)==0
    uwait(errordlg('The ''cands'' directory is empty.'));
    return
end

% Extract file names
allCandsList={candsList.name};

% Extract cands file names
pos=strncmp(allCandsList,'cands',5);
allCandsList=allCandsList(find(pos));

% Check that there are files left
if isempty(allCandsList)
    uwait(errordlg('The ''cands'' directory does not contain any ''cands###.mat'' files.'));
    return
end

% Change to cands directory
cd(candsDir);

% Start analysis
c=0;
len=length(allCandsList)
for i=2:len

    % Only in the first round we need to load and extract particles from two cands
    if i==2
        
        % Try do everything
        try
            
            % Extract positions from first cands
            first=load(allCandsList{1});
            candsFirst=first.cands;
            indxFirst=find([candsFirst.status]==1);
            posFirst=reshape([candsFirst(indxFirst).Lmax],2,length([candsFirst(indxFirst).Lmax])/2)';
            % Extract positions from second cands
            second=load([char(allCandsList{2})]);
            candsSecond=second.cands;
            indxSecond=find([candsSecond.status]==1);
            posSecond=reshape([candsSecond(indxSecond).Lmax],2,length([candsSecond(indxSecond).Lmax])/2)';

       catch
            
           % Something went wrong
           uiwait(errordlg('Could not extract speckle coordinates from ''cands''. Quitting.','Error','modal'));
           return
           
        end
        
    else
        
        % In subsequent rounds, the new first cands is the second of the previous round
        posFirst=posSecond;
        
        % Extract positions from the new second cands
        second=load(allCandsList{i});
        candsSecond=second.cands;
        indxSecond=find([candsSecond.status]==1);
        posSecond=reshape([candsSecond(indxSecond).Lmax],2,length([candsSecond(indxSecond).Lmax])/2)';
        
    end
    
    % Calculate all distances (but store only those <= maxRadius)
    D=createSparseDistanceMatrix(posFirst,posSecond,maxRadius);
    
    % In the first round, we initialize the output vector
    if i==2
    
        % Initialize output vector (in case it is not big enough, MATLAB will lengthen it)
        allDistances=zeros(1,size(D,1)*length(allCandsList));
        
    end
        
    % For each speckle, store the lowest distance to any other speckle in the following frame
    %    (unless there is none < maxRadius)
    for j=1:size(D,1)
        cRow=D(j,:);
        indx=find(cRow~=0);
        if ~isempty(indx)
            c=c+1;
            allDistances(c)=fix(min(cRow(indx))); % This sets the very small values replacing 0 in the sparse matrix to zero
        end
    end
    
    % Inform the user
    fprintf('Processed: %d/%d.\n',i-1,len-1);
    
end

% Crop allDistances if too much memory was allocated
if c<length(allDistances)
    allDistances=allDistances(1:c);
end

% Plot results in a histogram and a cumulative histogram
figH=figure;
h1=subplot(2,1,1); 
[n,b]=hist(allDistances,[-0.5:maxRadius+0.5]);
bar(b+0.5,n./sum(n))
mn=mean(allDistances);
sd=std(allDistances);
title(['Distances distribution (X=',num2str(mn),'; s=',num2str(sd),')']);

h2=subplot(2,1,2);
fraction=cumsum(n./sum(n));
h=bar(b+0.5,fraction);
set(h2,'YGrid','on')
set(h2,'YMinorGrid','on')
indx=find(fraction>0.95);
title(['Cumulative distribution: 95% of all distances are <= ',num2str(indx(1)),' pixels']);

% Change some figure attributes
set(figH,'Name','SpeckTackle - Search radius estimation');
set(figH,'NumberTitle','off');

% Change back to previous directory
cd(oldDir);
