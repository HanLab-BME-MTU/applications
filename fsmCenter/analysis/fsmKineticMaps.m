function [polyMap,depolyMap,img2C]=fsmKineticMaps(firstKinScore,imgSize,n,sigma)
% fsmKineticMaps creates maps of polymerization, depolymerization and net assembly rate
%
% SYNOPSIS      [polyMap,depolyMap,img2C]=fsmKineticMaps(firstKinScore,imgSize,n,sigma)
%
% INPUT         firstKinScore  : string containing the name (with complete path) of the 
%                                first kinScore###.mat file
%                                set firstKinScore=[] to have the function open a dialog
%               imgSize        : size of the analyzed images
%               n              : number of frames for time integration (0 for the entire stack)
%               sigma (pixels) : sigma for the low-pass filtering of the maps
%                                (default - sigma=5)
%
% OUTPUT        polyMap        : 2D polymerization map integrated over n frames
%               depolyMap      : 2D depolymerization map integrated over n frames
%               netMap         : 2D net assembly rate map integrated over n frames
%                                The map are always normalized with respect to the highest 
%                                score (positive or negative)
%
% DEPENDENCES   fsmKineticMaps uses { }
%               fsmKineticMaps is used by {  }
%
% Aaron Ponti, September 2th, 2003

if nargin<3 | nargin>4
    error('Three or four input parameters expected');
end

if nargin==3
    sigma=5;
end

if isempty(firstKinScore) | exist(firstKinScore)~=2 % Not a file
    
    % Select kinScore###.mat
    [fName,dirName] = uigetfile(...
        {'*.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select first kinScore matrix');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        polyMap=[];depolyMap=[];img2C=[];
        return 
    end
    
    firstKinScore=[dirName,fName];
end

% Recover all file names from the stack
outFileList=getFileStackNames(firstKinScore);

% Only consider first n frames if n~=0
len=length(outFileList);
if n~=0
    if n>len
        n=len;
    end
    outFileList=outFileList(1:n);
else
    n=len;
end

% Initialize emtpy maps
polyMap=zeros(imgSize);
depolyMap=polyMap;

% Create map
for i=1:n
    
    % Load kinScore
    load(char(outFileList(i)));

    % Read index
    [path,body,indxStr,ext]=getFilenameBody(char(outFileList(i)));
  
    % Copy matrix
    eval(['kinScore=kinScore',indxStr,'; clear kinScore',indxStr,';']);
    
    % Find scores
    indx=find(kinScore(:,2)~=0);
    
    % Copy scores into map
    if ~isempty(indx)
        
        for j=1:length(indx)
            
            % Read score
            score=kinScore(j,4);
            
            switch sign(score)
                case 1, polyMap(kinScore(j,2),kinScore(j,3))=polyMap(kinScore(j,2),kinScore(j,3))+score;
                case -1, depolyMap(kinScore(j,2),kinScore(j,3))=depolyMap(kinScore(j,2),kinScore(j,3))+score;
                otherwise
                    error('A zero score.');
            end
            
        end
        
    end
    
end

% Low-pass filter
if sigma~=0
    polyMap=Gauss2D(polyMap,sigma);
    depolyMap=Gauss2D(depolyMap,sigma);
end

% Create dual-channel image
mx=max([polyMap(:);abs(depolyMap(:))]);

img2C=zeros([imgSize 3]);
img2C(:,:,1)=polyMap/mx;
img2C(:,:,2)=abs(depolyMap)/mx;

