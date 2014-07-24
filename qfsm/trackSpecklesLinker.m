function [MPM,M]=trackSpecklesLinker(M,varargin)
% fsmTrackLinker creates the magic position matrix MPM from M
%
% SYNOPSIS      [MPM,M]=fsmTrackLinker(M)
%
% INPUT         M          : M stack as returned by the tracker functions
%                                  M = [y x y x]   [y x y x]   [y x y x]
%                                         ...    ,    ...    ,    ...
%                                       t1   t2     t2   t3     t3   t4
%                                          1           2           3  
% OUTPUT        MPM        : Magic Position Matrix 
%                                MPM = [ y  x  y  x  y  x ... ]
%                                         t1    t2    t3
%               M          : Rearranged M matrix.
%
% DEPENDENCES   fsmTrackLinker uses { }
%               fsmTrackLinker is used by { fsmTrackMain } 

% Aaron Ponti, 2002
% Sebastien Besson, June 2011
% Copied from fsmTrackLinker

% Check input
ip = inputParser;
ip.addRequired('M',@(x) size(x,3)>1);
ip.addOptional('waitbar',-1,@ishandle);
ip.parse(M,varargin{:});
wtBar= ip.Results.waitbar;

% Initialize waitbar
if ishandle(ip.Results.waitbar)
    waitbar(0,wtBar,'Linking gaps...');
elseif feature('ShowFigureWindows'),
    wtBar=waitbar(0,'Linking gaps...');
end

% Reorganize M
for counter1=1:size(M,3)-1
    
    % Read speckle positions at time point (=img) counter1
    start=M(:,3:4,counter1);
    stop=M(:,1:2,counter1+1);
    
    % Re-arrange stop (and therefore M) to correspond to the sequence of start
    tM=zeros(size(start,1),4);
    
    for counter2=1:size(stop,1)
        t=start(counter2,1)==stop(:,1);
        u=start(counter2,2)==stop(:,2);
        y=find(t & u);
        
        % ANALYSIS
        
        % No matching found -> error!
        if isempty(y)
            fprintf(1,'Time points %d to %d.\n',counter1,counter1+1);
            warning('fsmTrackLinker: Warning! Correspondance not found.');
            tM(counter2,:)=0; % -1;   
        end
        
        % Only one entry found in stop
        if length(y)==1
            tM(counter2,:)=M(y,:,counter1+1);
            stop(y,:)=-3;
        end
        
        % More than one entry found, but either 'no speckle' (0) 
        %   or 'already treated' (-3)
        if length(y)>1 && (start(counter2,1)~=0 || start(counter2,1)~=-3)
            tM(counter2,:)=M(y(1),:,counter1+1);
            stop(y(1),:)=-3;
        end
        
        % More than one repetition of a speckle found (~=0 & ~=-3)
        if  length(y)>1 && start(counter2,1)~=0 && start(counter2,1)~=-3
            error('fsmTrackLinker: Warning! Not all repetitions have been removed.');
        end        
        
    end
    
    % Replace M with re-ordered one
    M(:,:,counter1+1)=tM;
    
    if ishandle(wtBar), waitbar(counter1/(size(M,3)-1),wtBar); end
end

% Remove not needed info
MPM(:,1:2)=M(:,1:2,1);
for counter3=2:size(M,3)
    MPM(:,(counter3-1)*2+(1:2))=M(:,1:2,counter3);
end
MPM(:,counter3*2+(1:2))=M(:,3:4,counter3);

% Close waitbar if not-delegated
if isempty(ip.Results.waitbar), close(wtBar); end