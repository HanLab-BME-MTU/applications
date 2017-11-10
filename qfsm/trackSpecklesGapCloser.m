function [M,gapList]=trackSpecklesGapCloser(M,threshold,invalidCands,varargin)
% fsmTrackgapCloser closes gaps in speckles lifetimes by recovering weak
% local maxima close (in space and time) to a birth and death event
%
% SYNOPSIS      [M,gapClosed]=fsmTrackGapCloser(M,threshold,strg,userPath,firstIndex,frame)
%
% INPUT         M          : M stack as returned by the tracker functions
%                                  M = [y x y x]   [y x y x]   [y x y x]
%                                         ...    ,    ...    ,    ...
%                                       t1   t2     t2   t3     t3   t4
%                                          1           2           3  
%               threshold  : radius of the considered neighborhood
%               invalidCands       : String to correctly format the numeric suffix of saved files
%               userPath   : Work path
%               firstIndex : index of the first corresponding image
%               frame      : optional - by default, the gap closer goes
%                            through the whole M stack and closes all gaps.
%                            Alternatively, it can be used to close exclusively 
%                            the gaps at tn somewhere in the middle of the M matrix.
%                            The user passes, for tn, M(:,:,n:n+1) to the gap
%                            closer and must also set the variable frame=n.
%                            This will be used to load the corresponding
%                            cands###.mat file (see fsm PREPROCESSING module).
%                            If size(M,3) is > 2 the input parameter frame is IGNORED. 
%
% OUTPUT        M          : corrected M stack
%               gapClosed  : number of gaps closed
%

% Aaron Ponti, October 4th, 2002
% Sebastien Besson, June 2011
% Adapted from fsmTrackGapCloser

% Check input
ip = inputParser;
ip.addRequired('M',@(x) size(x,3)>1);
ip.addRequired('threshold',@isscalar);
ip.addRequired('invalidCands',@iscell);
ip.addParamValue('vectors',[],@iscell);
ip.addParamValue('waitbar',[],@ishandle);
ip.parse(M,threshold,invalidCands,varargin{:});
vectors = ip.Results.vectors;

%Initialize counters
matching=0;
notMatching=0;

% Initialize gapList
gapList = cell(1,size(M,3)+1);
gapList{1} = [0 0];
gapList{end} = [0 0];

if ~isempty(ip.Results.waitbar)
    h=ip.Results.waitbar;
    waitbar(0,h,'Closing gaps...');
else
    h=waitbar(0,'Closing gaps...');
end

for i=1:size(M,3)-1       
    % Load cands structure
    cands = invalidCands{i+1};
    % No insignificant local maxima, no need to try and close gaps
    if isempty(cands), gapList{i+1}=[0 0]; continue;  end

    if ~isempty(vectors)
        % Extract source and target positions to be (back)propagated
        source=M(:,1:2,i);
        target=M(:,3:4,i+1);
        
        % Keep track of the row in source and target which contain the speckles
        fS=find(source(:,1)~=0);
        fT=find(target(:,1)~=0);
        
        % Forward propagate
        pSource=propagateSpecklePositions(source(fS,:),vectors{i},'forward');
        % Backward propagate
        pTarget=propagateSpecklePositions(target(fT,:),vectors{i+1},'backward');
        
        % Create the copy of M (only three time-points)
        Mcopy=M(:,:,i:i+1);
        Mcopy(fS,1:2,1)=pSource;
        Mcopy(fT,3:4,2)=pTarget;
            
        % Pass the copy to the gap closing subfunction
        [pMcopy,glist,match,notMatch] = closeGap(Mcopy,1,cands,threshold);
        
        % Replace propagated positions with original positions and write back
        % to M (with gap closed)
        M(:,:,i:i+1)=pMcopy;
        M(:,1:2,i)=source;
        M(:,3:4,i+1)=target;
    else
        [M,glist,match,notMatch] = closeGap(M,i,cands,threshold);
    end
    % Update counters and gapList cell array
    gapList{i+1} = glist;
    matching = matching+match;
    notMatching = notMatching+notMatch;
    waitbar(i/(size(M,3)-1),h);
end

% Display results
tot=matching+notMatching;
if tot~=0
    fprintf(1,'Closed gaps: %d/%d (%d%%)\n',matching,tot,fix(matching/tot*100));
else
    fprintf(1,'No gaps found!!!\n');
end

% Close waitbar if not-delegated
if isempty(ip.Results.waitbar), close(h); end

function [M,gapList,matching,notMatching] = closeGap(M,iFrame,cands,threshold)

% Initialize counters
matching=0;
notMatching=0;
gapList=[];

% Read current time-point (first = 2)
P=M(:,1:2,iFrame);   % previous time-point          ( P )
C=M(:,3:4,iFrame);   % current time-point           ( C )
T=M(:,1:2,iFrame+1); % to be modified to close gaps ( T )
N=M(:,3:4,iFrame+1); % next time-point              ( N )
    
% Positions of non-speckles (in C)
posC=C(:,1)==0;
posC=find(posC);

% Check for existance of speckles at these positions in P
% posP=P(posC,1)~=0;
% posP=find(posP);
posP=posC(P(posC,1)~=0);

% Read the speckle coordinates in P
pM=P(posP,:);

% Positions of non-speckles (in T)
posT=T(:,1)==0;
posT=find(posT);

% Check for existance of speckles at these positions in N
posN=N(posT,1)~=0;
posN=posT(posN);

% Read the speckle coordinates in N
nM=N(posN,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check if some of the speckles in N at the positions given
%     by pos match those given by posP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(pM) && ~isempty(nM)    % There are candidate speckles
    
    % Create distance matrix
    D=createDistanceMatrix(pM,nM);
    
    % Row
    E=zeros(size(D));
    for i=1:size(D,1)
        t=D(i,:);
        t=t==min(t);
        E(i,:)=t;
    end
    
    % Column
    F=zeros(size(D));
    for i=1:size(D,2)
        t=D(:,i);
        t=t==min(t);
        F(:,i)=t;
    end
    
    % Thresholding
    H=D<=threshold;
    
    % Resulting selected distance matrix
    G=(E & F & H);
    
    % Now find matching pairs: y = speckle in pM, x = speckle in nM
    [y x]=find(G);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If there exist some speckles from N matching some speckles from P, build a
    %    matching matrix (PM)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(y) && ~isempty(x)   % Speckles have to be paired
        
        % Sort y and x (preserving coupling!) to avoid future problems
        PM=[];
        PM(1:length(y),1)=y;
        PM(1:length(x),2)=x;
        PM=sortrows(PM,1:2);
        y=PM(:,1);
        x=PM(:,2);
        PM=[];
        
        % Now build matching pairs
        len=length(unique(sort(y))); % Must be sorted for unique, otherwise you will lose some of the speckles!
        
        counter=1;
        
        for i=1:len
            nY=y;                      % Create a copy of y for boolean testing
            posY=y(counter);
            nY=nY==posY;
            indx=find(nY);
            PM(i,1)=posY;
            PM(i,2:1+length(indx))=x(indx)';
            counter=counter+length(x(indx)');
        end
        clear counter;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Check whether a speckle in nM is bound to more than one speckle in pM
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        RepetitionList=[];
        numberOfRepetitions=0;
        
        source=unique(sort(PM(:,2)));
        target=PM(:,2);
        if length(source)<length(target)
            counter=0;
            for cn=1:length(source)
                %
                ss=source(cn);
                ss=ss==target(:);
                ss=find(ss);
                if length(ss)>1
                    % Keep the first repeating speckle, mark all others to be deleted
                    ss(1,:)=[];
                    counter=counter+1;
                    RepetitionList(counter:counter+length(ss)-1,1)=ss;
                    counter=counter+length(ss)-1;
                    numberOfRepetitions=numberOfRepetitions+1;
                end
            end
            % Remove marked repeating speckles
            if RepetitionList~=0
                PM(RepetitionList,:)=[];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Start closing gaps
        %    sD = [y x] : coordinates of the speckle to be written into M
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clear D E F G H
        
        %
        % Go through PM
        %
        
        for c2=1:size(PM,1)
            
            % Initialize sD
            sD=[];
            
            % Go through speckles
            t=PM(c2,2:end)>0;
            t=find(t);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % MORE THAN TWO SPECKLES FROM nM MATCH THE SPECKLE FROM pM
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if length(t)>2
                
                fprintf(1,'\nWarning -> gapCloser: one speckle death is matched by %d births.\n',length(t));
                fprintf(1,'One will be selected without further analysis.\n');
                % This problem is reduced to a 1 to 1 case
                PM(c2,3:end)=0; % Remove index of death events 2 and above
                t=1;                     % Set the number of death events found to 1
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % ONLY ONE SPECKLE FROM nM MATCHES THE SPECKLE FROM pM
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if length(t)==1
                
                % Only one matching speckle
                st=pM(PM(c2,1),:);
                
                % Look for a loc max in the cands structure close to st
                tmp=[cands.Lmax];
                tmp=reshape(tmp,2,length(tmp)/2)';
                D=createDistanceMatrix(tmp,st);
                
                % Find minimum distance (which has also to be < a given threshold)
                E=(D<threshold & D==min(D));
                E=find(E);
                
                if ~isempty(E)
                    
                    % Check that only one speckle is at minimum distance
                    if length(E)>1
                        
                        for i=1:length(E)
                            sD(i,1:2)=cands(E(i)).Lmax;
                        end
                        
                        % Use also the speckle from N to build a distance criterium
                        sp=nM(PM(c2,2),:);
                        
                        % Select one of the speckles from the cands structure
                        for i=1:length(E)
                            D2(i,1)=sqrt((sD(i,1)-st(1))^2+(sD(i,2)-st(2))^2)*... % Distance to st
                                sqrt((sD(i,1)-sp(1))^2+(sD(i,2)-sp(2))^2);         % Distance to sp
                        end
                        
                        % Select min distance
                        D2=D2==min(D2);
                        
                        % Check that these speckles are not already coupled (in T)
                        w=[];
                        for i=1:length(E)
                            v2=[];
                            t2=sD(i,1)==T(:,1);
                            u2=sD(i,2)==T(:,2);
                            v2=find(t2 & u2);
                            if isempty(v2)
                                w(i,1)=1;
                            end
                            if length(v2)==1
                                w(i,1)=0;
                            end
                            if length(v2)>1
                                error('This SHOULDN''T OCCUR!');
                            end
                        end
                        
                        % Insert this info into the criterion
                        D2=w.*D2;
                        
                        % Get speckle
                        D2=find(D2);
                        
                        % Pick out speckle
                        if isempty(D2), sD=[]; end
                        if length(D2)==1
                            sD=sD(D2,:);
                            % "Remove" assigned speckle from the cands structure
                            cands(E(D2)).Lmax=[Inf Inf];
                        end
                        if length(D2)>1
                            D2=D2(1);
                            sD=sD(D2,:);
                            
                            % "Remove" assigned speckle from the cands structure
                            cands(E(D2)).Lmax=[Inf Inf];
                        end
                        
                    else
                        
                        % Valid match - Gap successfully closed
                        sD=cands(E).Lmax;
                        
                        % "Remove" assigned speckle from the cands structure
                        cands(E).Lmax=[Inf Inf];
                        
                    end
                    
                else
                    
                    %disp('No candidate speckle found in the cands structure');
                    sD=[];
                    
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % TWO SPECKLES FROM nM MATCH THE SPECKLE FROM pM
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if length(t)==2
                
                %
                % SPECKLE 1
                %
                
                % First matching speckle
                st=pM(PM(c2,1),:);
                
                % Look for a loc max in the cands structure close to st
                tmp=[cands.Lmax];
                tmp=reshape(tmp,2,length(tmp)/2)';
                D=createDistanceMatrix(tmp,st);
                
                % Find minimum distance (which has also to be < a given threshold)
                E=(D<threshold & D==min(D));
                E=find(E);
                
                if ~isempty(E)
                    
                    % Check that only one speckle is at minimum distance
                    if length(E)>1
                        
                        for i=1:length(E)
                            sD(i,1:2)=cands(E(i)).Lmax;
                        end
                        
                        % Use also the speckle from N to build a distance criterium
                        sp=nM(PM(c2,2),:);
                        
                        % Select one of the speckles from the cands structure
                        for i=1:length(E)
                            D2(i,1)=sqrt((sD(i,1)-st(1))^2+(sD(i,2)-st(2))^2)*... % Distance to st
                                sqrt((sD(i,1)-sp(1))^2+(sD(i,2)-sp(2))^2);         % Distance to sp
                        end
                        
                        % Select min distance
                        D2=D2==min(D2);
                        
                        % Check that these speckles are not already coupled (in T)
                        w=[];
                        for i=1:length(E)
                            v2=[];
                            t2=sD(i,1)==T(:,1);
                            u2=sD(i,2)==T(:,2);
                            v2=find(t2 & u2);
                            if isempty(v2)
                                w(i,1)=1;
                            end
                            if length(v2)==1
                                w(i,1)=0;
                            end
                            if length(v2)>1
                                error('This SHOULDN''T OCCUR!');
                            end
                        end
                        
                        % Insert this info into the criterion
                        D2=w.*D2;
                        
                        % Get speckle
                        D2=find(D2);
                        
                        % Pick out speckle
                        if length(D2)==0
                            sD1=[]; E1=[]; sD=[];
                        end
                        if length(D2)==1
                            sD1=sD(D2,:); E1=E(D2); sD=[];
                        end
                        if length(D2)>1
                            D2=D2(1); sD1=sD(D2,:); E1=E(1); sD=[];
                        end
                        
                    else
                        
                        % First candidate speckle (sD1)
                        sD1=cands(E).Lmax;
                        E1=E;
                        
                    end
                    
                else
                    
                    % NO CANDIDATE SPECKLE FOUND
                    sD1=[];
                    E1=E;
                    
                end
                
                %
                % SPECKLE 2
                %
                
                clear D D2 E F G H
                sD=[];
                
                % First speckle is always the same
                st=pM(PM(c2,1),:);
                
                % Look for a loc max in the cands structure data close to st
                tmp=[cands.Lmax];
                tmp=reshape(tmp,2,length(tmp)/2)';
                D=createDistanceMatrix(tmp,st);
                
                % Find minimum distance (which has also to be < a given threshold)
                E=(D<threshold & D==min(D));
                E=find(E);
                
                if ~isempty(E)
                    
                    % Check that only one speckle is at minimum distance
                    if length(E)>1
                        
                        for i=1:length(E)
                            sD(i,1:2)=cands(E(i)).Lmax;
                        end
                        
                        % Use also the speckle from N to build a distance criterium
                        sp=nM(PM(c2,3),:);
                        
                        % Select one of the speckles from the cands structure
                        for i=1:length(E)
                            D2(i,1)=sqrt((sD(i,1)-st(1))^2+(sD(i,2)-st(2))^2)*... % Distance to st
                                sqrt((sD(i,1)-sp(1))^2+(sD(i,2)-sp(2))^2);         % Distance to sp
                        end
                        
                        % Select min distance
                        D2=D2==min(D2);
                        
                        % Check that these speckles are not already coupled (in T)
                        w=[];
                        for i=1:length(E)
                            v2=[];
                            t2=sD(i,1)==T(:,1);
                            u2=sD(i,2)==T(:,2);
                            v2=find(t2 & u2);
                            if isempty(v2)
                                w(i,1)=1;
                            end
                            if length(v2)==1
                                w(i,1)=0;
                            end
                            if length(v2)>1
                                error('This SHOULDN''T OCCUR!');
                            end
                        end
                        
                        % Insert this info into the criterion
                        D2=w.*D2;
                        
                        % Get speckle
                        D2=find(D2);
                        
                        % Pick out speckle
                        if isempty(D2)
                            sD2=[]; E2=[]; sD=[];
                        end
                        if length(D2)==1
                            sD2=sD(D2,:); E2=E(D2); sD=[];
                        end
                        if length(D2)>1
                            D2=D2(1); sD2=sD(D2,:); E2=E(1); sD=[];
                        end
                        
                    else
                        
                        % Second candidate speckle (sD2)
                        sD2=cands(E).Lmax;
                        E2=E;
                        
                    end
                    
                else
                    
                    % NO CANDIDATE SPECKLE FOUND
                    sD2=[];
                    E2=E;
                    
                end
                
                %
                % SELECT BETWEEN SPECKLE 1 AND 2
                %
                
                % Check that speckle 1 is not already coupled (in T)
                if ~isempty(sD1)
                    t1=sD1(1,1)==T(:,1);
                    u1=sD1(1,2)==T(:,2);
                    v1=(t1 & u1);
                    
                    if ~isempty(find(v1,1)), sD1=[]; end
                end
                
                
                % Check that speckle 2 is not already coupled (in T)
                if ~isempty(sD2)
                    t2=sD2(1,1)==T(:,1);
                    u2=sD2(1,2)==T(:,2);
                    v2=(t2 & u2);
                    
                    if ~isempty(find(v2,1)), sD2=[]; end
                end
                
                % Select
                if isempty(sD1) && isempty(sD2)
                    sD=[];
                    E=[];
                end
                if isempty(sD1) && ~isempty(sD2)
                    sD=sD2;
                    E=E2;
                end
                if ~isempty(sD1) && isempty(sD2)
                    sD=sD1;
                    E=E1;
                end
                if ~isempty(sD1) && ~isempty(sD2)
                    D12(1,1)=sqrt((sD1(1,1)-st(1))^2+(sD1(1,2)-st(2))^2);
                    D12(2,1)=sqrt((sD2(1,1)-st(1))^2+(sD2(1,2)-st(2))^2);
                    D12=D12==min(D12);
                    D12=find(D12);
                    if length(D12)==1
                        if D12==1
                            sD=sD1;
                            E=E1;
                        end
                        if D12==2
                            sD=sD2;
                            E=E2;
                        end
                    end
                    if length(D12)==2
                        % Simply take one of them
                        sD=sD1;
                        E=E1;
                    end
                end
                
                % Remove from the cands structure
                if ~isempty(E)
                    cands(E).Lmax=[Inf Inf];
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % SELECTED SPECKLE HAS NOW TO BE WRITTEN INTO M
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isempty(sD)
                
                % Check that this speckle is not already coupled (in T)
                %    (this is redundant in the case that more than one speckle were found in the cands structure)
                t=sD(1,1)==T(:,1);
                u=sD(1,2)==T(:,2);
                v=(t & u);
                
                if isempty(find(v,1))
                    
                    % Another speckle to be added to the list of "closed gaps"
                    matching=matching+1;
                    
                    % Valid match - write st to M (to close gap)
                    M(posP(PM(c2,1)),3:4,iFrame)=sD;
                    M(posN(PM(c2,2)),1:2,iFrame+1)=sD;

                    % Write this speckle as 'gap'
                    gapList(matching,1:2)=sD;      
                else
                    % Update counter
                    notMatching=notMatching+1;    
                end     
            else
                % Update counter
                notMatching=notMatching+1;      
            end         
        end
    end    
end