function [M,gapClosed]=fsmTrackGapCloser(M,threshold,strg,userPath,firstIndex,frame)
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
%               strg       : String to correctly format the numeric suffix of saved files
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
% DEPENDENCES   fsmTrackgapCloser uses { }
%               fsmTrackgapCloser is used by { fsmTrackMain } 
%
% Aaron Ponti, October 4th, 2002

% Set flag
ENHANCED=1;
if nargin==5
    ENHANCED=0;
end

% Initialize counters
matching=0;
notMatching=0;

% Number of iterations
n=size(M,3);

if ENHANCED==0
    % Initializing waitbar
    h=waitbar(0,'Closing gaps...');
end

%
% START
%

% Constant
Mdepth=size(M,3);

% if Mdepth==2 & isempty(frame)
%     error('A value for the parameter frame is expected');
% end

for c1=1:Mdepth-1
    
    % Current index
    currentIndex=c1+firstIndex-1;

    % Reset number of gaps per time point and gapList
    gapN=0; gapList=[];
    
    clear D E F G H t u v
    PM=[];
   
    % Which file should be loaded?
    if ENHANCED==0
        frame=currentIndex+1;
    end
        
    % Load cands structure
    indxStr=sprintf(strg,frame);
    %eval(strcat('load cands',filesep,'cands',indxStr,'.mat;'));
    eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']); % Cands workspace

    % Only consider statistically insignificant local maxima
    toBeExtr=find([cands.status]==0);
    if ~isempty(toBeExtr)
        cands=cands(toBeExtr);
    else
        cands=[];
    end
    
    if isempty(cands) % No insignificant local maxima, no need to try and close gaps
        % Save gap list
        gapList=[0 0];
        %eval(strcat('save gapList',filesep,'gapList',indxStr,'.mat gapList;'));
        eval(['save ',userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat gapList;']); 
        gapList=[];
        if n>2
            % Update waitbar
            waitbar(c1/n,h);
        end
        continue;   % Jump to next iteration
    end
    
    % Read current time-point (first = 2)
    P=M(:,1:2,c1);   % previous time-point          ( P )
    C=M(:,3:4,c1);   % current time-point           ( C )
    T=M(:,1:2,c1+1); % to be modified to close gaps ( T )
    N=M(:,3:4,c1+1); % next time-point              ( N )
    
    % Positions of non-speckles (in C)
    posC=C(:,1)==0;
    posC=find(posC);
    
    % Check for existance of speckles at these positions in P
    posP=P(posC,1)~=0;
    posP=find(posP);
    posP=posC(posP);
    
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
    
    if ~isempty(pM) & ~isempty(nM)    % There are candidate speckles
        
        % Create distance matrix
        D=createDistanceMatrix(pM,nM);
        
        % Row
        for i=1:size(D,1)
            t=D(i,:);
            t=t==min(t);
            E(i,:)=t;
        end
        
        % Column
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
        
        if ~isempty(y) & ~isempty(x)   % Speckles have to be paired
            
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
                            if length(D2)==0
                                sD=[];
                            end
                            if length(D2)==1
                                sD=sD(D2,:);
                                % "Remove" assigned speckle from the cands structure
                                cands(E(D2)).Lmax=[Inf Inf];
                            end
                            if length(D2)>1
                                %disp(c1);
                                %disp('Even after further selection, there are still more than one matching speckle!');
                                %disp('First one selected.');
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
                            if length(D2)==0
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
                        v1=find(v1);
                        
                        if ~isempty(v1)
                            sD1=[];
                        end
                    end
                    
                    
                    % Check that speckle 2 is not already coupled (in T)
                    if ~isempty(sD2)
                        t2=sD2(1,1)==T(:,1);
                        u2=sD2(1,2)==T(:,2);
                        v2=(t2 & u2);
                        v2=find(v2);
                        
                        if ~isempty(v2)
                            sD2=[];
                        end
                    end
                    
                    % Select
                    if isempty(sD1) & isempty(sD2)
                        sD=[];
                        E=[];
                    end
                    if isempty(sD1) & ~isempty(sD2)
                        sD=sD2;
                        E=E2;
                    end
                    if ~isempty(sD1) & isempty(sD2)
                        sD=sD1;
                        E=E1;
                    end
                    if ~isempty(sD1) & ~isempty(sD2)
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
                    v=find(v);
                    
                    if isempty(v)
                        
                        % Another speckle to be added to the list of "closed gaps"
                        gapN=gapN+1;
                        
                        % Valid match - write st to M (to close gap)
                        M(posP(PM(c2,1)),3:4,c1)=sD;
                        M(posN(PM(c2,2)),1:2,c1+1)=sD;
                        
                        % Update counter
                        matching=matching+1;
                        
                        % Write this speckle as 'gap'
                        gapList(gapN,1:2)=sD;
                        
                    else
                        
                        %disp('Candidate speckle already assigned');
                        
                        % Update counter
                        notMatching=notMatching+1;
                        
                    end
                    
                else
                    
                    % Update counter
                    notMatching=notMatching+1;
                    
                end
                
            end
            
        else
            
            %disp('y and/or x empty');
            
        end
        
    else
        
        %disp('pM and/or nM empty');
        
    end
    
    % Save gap list
    if isempty(gapList)
        gapList=[0 0];
    end
    %eval(strcat('save gapList',filesep,'gapList',indxStr,'.mat gapList;'));
    eval(['save ',userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat gapList;']); 
    gapList=[];
    
    if ENHANCED==0
        % Update waitbar
        waitbar(c1/n,h);
    end
    
end

if ENHANCED==0
    % Close waitbar
    close(h);
end

% The following is only needed when the whole M is processed
if ENHANCED==0

    % Save also gapList for the first and last timepoints
    gapList=[0 0];
    indxStr=sprintf(strg,firstIndex);
    gapListFirstFileName=[userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat'];
    if ~exist(gapListFirstFileName); % Do not overwrite if it exists from a previous experiment
        eval(['save ',gapListFirstFileName,' gapList;']); 
    end
    
    indxStr=sprintf(strg,currentIndex+2);
    gapListLastFileName=[userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat'];
    if ~exist(gapListLastFileName); % Do not overwrite if it exists from a previous experiment
        eval(['save ',gapListLastFileName,' gapList;']); 
    end
    
end

% Display
tot=matching+notMatching;
if tot~=0
    fprintf(1,'Closed gaps: %d/%d (%d%%)\n',matching,tot,fix(matching/tot*100));
else
    fprintf(1,'No gaps found!!!\n');
end

% Return info
if matching==0
    % No gap closed
    gapClosed=0;
else
    gapClosed=1;
end