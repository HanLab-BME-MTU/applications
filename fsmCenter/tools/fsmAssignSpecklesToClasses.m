function speckleClasses=fsmAssignSpecklesToClasses(speckleArray)
% fsmAssignSpecklesToClasses assigns each speckle in speckleArray to one of 13 classes
% 
% SYNOPSIS   speckleClasses=fsmAssignSpecklesToClasses(speckleArray)
%
% The speckleArray structure is analyzed and a new, smaller structure, speckleClasses is returned 
% with following fields:
%
% speckleClasses
%                 .bTime   : time of actual birth (not of 'b' speckles, tb+1)
%                 .dTime   : time of actual death (not of 'd' speckles, td-1)
%                 .first   : index of 'b' speckle in speckleArray
%                 .first   : index of 'd' speckle in speckleArray
%                 .pos     : [y x     Speckle coordinate at actual birth
%                             y x]    Speckle coordinate at actual death
%                 .class   : one of ten possible classes (see below)
%                 .network : indicates the network (lp or la) to which the speckle belongs
%                            this field is NOT FILLED in this function, but in assignSpecklesToNetworks
%
% .class   
%  1 : significant BIRTH due to POLYMERIZATION,   significant DEATH due to DEPOLYMERIZATION
%  2 : significant BIRTH due to POLYMERIZATION,   significant DEATH due to POLYMERIZATION
%  3 : significant BIRTH due to DEPOLYMERIZATION, significant DEATH due to DEPOLYMERIZATION
%  4 : significant BIRTH due to DEPOLYMERIZATION, significant DEATH due to POLYMERIZATION
%  5 : significant BIRTH due to POLYMERIZATION,   non significant DEATH
%  6 : significant BIRTH due to DEPOLYMERIZATION, non significant DEATH
%  7 : non significant BIRTH,                     significant DEATH due to DEPOLYMERIZATION
%  8 : non significant BIRTH,                     significant DEATH due to POLYMERIZATION
%  9 : non significant BIRTH,                     non significant DEATH
% 10 : ghost speckle, a subset of class 9.
% 11 : speckle already present in frame 1 which dies within the movie
% 12 : speckle born within the movie and still present at movie end
% 13 : speckle which lives for the whole movie
%
% Aaron Ponti, 04/09/2004

% Check input
if nargin~=1
    error('One input parameter expected');
end

% Total number of speckles
total=length(speckleArray);

% Pre-allocate memory
tot=0.5*(length(find([speckleArray.status]=='b'))+length(find([speckleArray.status]=='d'))+length(find([speckleArray.status]=='f'))+length(find([speckleArray.status]=='l')));
speckleClasses=repmat(struct('bTime',0,'dTime',0,'first',0,'last',0,'pos',0,'class',0,'network',0),1,tot);
% Remark: the .network field is not used in this function, it will be used in assignSpecklesToNetworks

% Initialization
count=0;
lastEv='n';
currentB=0;

% Calculate some longer step size (not to spend too much time updating the waitbar)
step=0.5*10^fix(log10(total)-1); if step<1, step=1; end

% Initialize the waitbar
wH=waitbar(0,'Please wait...');

% Start
for i=1:total
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % INTERCEPT SPECKLE BEGINNING
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Look for a birth   
    if speckleArray(i).status=='b'
        
        % If last speckle was not complete, it will be discarded by this one
        currentB=i; % Mark this as the latest birth met
        switch speckleArray(i).activity
            case  1, BP=1;
            case  0, BP=0;
            case -1, BP=-1;
            otherwise
                error('Invalid activity');
        end            
        lastEv='b';
        
    end
    
    % Look for a speckles without birth   
    if speckleArray(i).status=='f'
        
        currentB=i; % Mark this as the latest "birth" met
        BP=-2; % 'f' speckles have no activity
        lastEv='f';
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % INTERCEPT SPECKLE END
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Look for a death
    if speckleArray(i).status=='d'
        
        if lastEv=='b' % Last event was a birth, this is a complete speckle
            switch speckleArray(i).activity
                case  1, DP=1;
                case  0, DP=0;
                case -1, DP=-1;
                otherwise
                    error('Invalid activity');
            end
            
            % Define class of speckle
            switch BP
                case 1
                    if DP==-1, class=1; end % BP, DD
                    if DP== 0, class=5; end % BP, insignificant death
                    if DP== 1, class=2; end % BP, DP
                case 0
                    if DP==-1, class=7; end % insignificant birth, DD
                    if DP== 0, class=9; end % insignificant birth, insignificant death
                    if DP== 1, class=8; end % insignificant birth, DP
                case -1
                    if DP==-1, class=3; end % BD, DD  
                    if DP== 0, class=6; end % BD, insignificant death
                    if DP== 1, class=4; end % BD, DP           
                otherwise
                    error('Wrong class of birth.');
                end
                
                % Discriminate Ghost speckles
                if class==9 & (i-currentB)==2
                    class=10;
                end
                                
        elseif lastEv=='f' % Last event was a 'f' this is a speckle with no birth
                
            class=11;
              
        else
                
            error('A ''d'' speckle can be preceded only by either a ''b'' or an ''f'' speckle.');
                
        end
            
        % Save current information into the speckle structure
        count=count+1;
        if count>tot
            disp('Not enough space allocated. Reallocating...'); % This should never happen
        end
        speckleClasses(count).bTime = speckleArray(currentB).timepoint+1; % Store actual birth time (not of 'b' speckle)
        speckleClasses(count).dTime = speckleArray(i).timepoint-1; % Store actual death time (not of 'd' speckle)      
        speckleClasses(count).first = currentB;
        speckleClasses(count).last  = i;
        speckleClasses(count).pos   = [speckleArray(currentB).spPos; speckleArray(i).spPos]; %reshape([speckleArray(currentB:i).spPos],2,1+(i-currentB))'; 
        speckleClasses(count).class = class;

        % Mark last event as 'd'
        lastEv='d';    
        
    end
        
    if speckleArray(i).status=='l'
        
        if lastEv=='b' % This is a speckle with no birth
            class=12;
        elseif lastEv=='f' % This is a speckle living for the whole movie
            class=13;
        else
            error('An ''l'' speckle can be preceded only by either a ''b'' or an ''f'' speckle.');
        end
        
        % Save current information into the speckle structure
        count=count+1;
        if count>tot
            disp('Not enough space allocated. Reallocating...'); % This should never happen
        end
        speckleClasses(count).bTime = speckleArray(currentB).timepoint+1; % Store actual birth time (not of 'b' speckle)
        speckleClasses(count).dTime = speckleArray(i).timepoint-1; % Store actual death time (not of 'd' speckle)      
        speckleClasses(count).first = currentB;
        speckleClasses(count).last  = i;
        speckleClasses(count).pos   = [speckleArray(currentB).spPos; speckleArray(i).spPos]; %reshape([speckleArray(currentB:i).spPos],2,1+(i-currentB))';
        speckleClasses(count).class = class;
        
        % Mark last event as 'l'
        lastEv='l';    
        
    end
    
    % Update waitbar if needed
    if step>1
        if mod(i,step)==1
            waitbar(i/total,wH);
        end
    else
        waitbar(i/total,wH);
    end        
    
end

% Close waitbar
close(wH);
