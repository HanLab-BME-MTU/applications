function stats=fsmInfoSpeckleArrayStatistics(speckleArray,writeToFile)
% fsmInfoSpeckleArrayStatistics performs statistical analysis on the speckleArray structure
%
% SYNOPSIS   stats=fsmInfoSpeckleArrayStatistics(speckleArray,writeToFile)
%
% INPUT      speckleArray       :  structure returned by the builder module
%            writeToFile [ 0|1 ]:   if 0, writes statistics to console
%                                   if 1, writes to console and to disk (file 'stats.log'
%                                   in the current directory)
%
% OUTPUT     stats              : structure containing following fields
%
%                                 .bSpeckle        : number of speckles without death
%                                 .dSpeckle        : number of speckles without birth
%                                 .gCounter        : number of 'g' speckles (gaps)
%                                 .sCounter        : number of 's' speckles
%                                 .numberOfGhost   : number of ghost speckles 
%                                                    (lifetime = 1)
%                                 .numberOfSpeckle : number of speckles with lifetime > 1
%                                 .polyScore       : strongest polymerization score
%                                 .depolyScore     : strongest depolymerization score
%                                 .pScores         : vector of polymerization scores
%                                 .dScores         : vector of depolymerization scores
%                                 .meanLifeTime    : mean speckle lifetime (ghost speckles
%                                                    excluded)
%                                 .events          : number of 'b' and 'd' speckles
%                                 .weakEvents      : number of 'b' and 'd' speckles
%                                                    corresponding to insignificant local 
%                                                    maxima
%                                 .speckleDI       : vector containing all deltaI for
%                                                    long-living speckles (lifetime > 1)
%                                 .ghostDI         : vector containing all deltaI for
%                                                    ghost speckles
% 
% DEPENDENCES   fsmInfoSpeckleArrayStatistics uses {}
%               fsmInfoSpeckleArrayStatistics is used by { fsmMain }
%
% Aaron Ponti, May 26th, 2003


% Initialize result structure
stats=struct('complete',0,...
    'bSpeckle',0,...
    'dSpeckle',0,...
    'gCounter',0,...
    'sCounter',0,...
    'numberOfGhost',0,...
    'numberOfSpeckle',0,...
    'polyScore',0,...
    'depolyScore',0,...
    'pScores',0,...
    'dScores',0,...,
    'meanLifeTime',0,...
    'events',0,...
    'weakEvents',0,...
    'speckleDI',0,...
    'ghostDI',0);    

% Initialization: (1)
last='n';

% Initialization: (3)
polCounter=0;
depolCounter=0;

% Initialization: (4)
count=0;
lastEv='n'; posB=0; 
lifetime=[];

for i=1:length([speckleArray.timepoint])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % (1) Number of complete speckle objects, of speckles without birth and speckles without death
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if speckleArray.status(i)=='l'
        % Speckle without death
        stats.bSpeckle=stats.bSpeckle+1;
        last='l';
    end
    if speckleArray.status(i)=='f'
        % Speckle without birth
        stats.dSpeckle=stats.dSpeckle+1;
        last='f';
    end
    if speckleArray.status(i)=='d'
        if last=='b'
            % A complete speckle
            stats.complete=stats.complete+1;
        end
        last='d';
    end
    if speckleArray.status(i)=='b'
        last='b';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % (2) Number of 'g' (gaps) and 's','f','l' speckles
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if speckleArray.status(i)=='g'
        stats.gCounter=stats.gCounter+1;
    end
    if speckleArray.status(i)=='s' | speckleArray.status(i)=='f' | speckleArray.status(i)=='l'
        stats.sCounter=stats.sCounter+1;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % (3) Strongest poly and depoly scores
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % POLY
    if speckleArray.activity(i)==1
        % Store a positive score
        polCounter=polCounter+1;
        stats.pScores(polCounter)=speckleArray.score(i);
        % Check whether this score is the strongest positive score
        if speckleArray.score(i)>stats.polyScore
            stats.polyScore=speckleArray.score(i);
        end
    end
    % DEPOLY
    if speckleArray.activity(i)==-1
        % Store a negative score
        depolCounter=depolCounter+1;
        stats.dScores(depolCounter)=speckleArray.score(i);
        % Check whether this score is the strongest negative score
        if speckleArray.score(i)<stats.depolyScore
            stats.depolyScore=speckleArray.score(i);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % (4) Mean speckle lifetime
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if speckleArray.status(i)=='b' & speckleArray.activity(i)~=0
        if lastEv=='b' | lastEv=='f' | lastEv=='l'
            % Forget this
            lastEv='n';
            posB=-1;
        else
            lastEv='b';
            posB=i;
        end
    end
    if speckleArray.status(i)=='d'
        if speckleArray.activity(i)~=0
            if lastEv=='b' & posB~=-1
                if ((i-1)-posB)>1 % Count them only if the trajectory is more than one frame (i.e. not a ghost speckle)
                    count=count+1;
                    lifetime(count)=(i-1)-posB;
                end
            elseif lastEv=='d' | lastEv=='f' | lastEv=='l'
                % A speckle without birth - forget it
                posB=-1;
            end
            lastEv='d';
        end
        lastEv='d';
    end
    if speckleArray.status(i)=='f' |speckleArray.status(i)=='l'
        lastEv='n';
        posB=-1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % (5) Events and weak events
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % For this analysis we only consider SIGNIFICANT events
    if (speckleArray.status(i)=='b' | speckleArray.status(i)=='d') & speckleArray.activity(i)~=0
        stats.events=stats.events+1;
        if speckleArray.lmEvent(i)~=0
            stats.weakEvents=stats.weakEvents+1;
        end
    end
    
end

% Calculate mean lifetime
stats.meanLifeTime=mean(lifetime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ghost speckles vs. longer-living speckles - distributions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last event check
foundB=0; foundD=0;
for i=length([speckleArray.timepoint]):-1:1
    if speckleArray.status(i)=='b' & foundB==0
        lastB=i;
        foundB=1;
        continue;
    end
    if speckleArray.status(i)=='d' & foundD==0
        lastD=i;
        foundD=1;
        continue;
    end
    if foundB==1 & foundD==1
        break;
    end
end
lastE=min(lastB,lastD);

i=1; gC=0; nC=0;
while i<=lastB % To stay within speckleArray with the check
	if i>lastE
        break;
    end
    clear tb td;
    % Found a birth
    if speckleArray.status(i)=='b'
        % Store tb
        tb=i;
        %
        i=i+1;
        n=0;
        deltaIs=0;
        
        % Go along the lifetime
        while speckleArray.status(i)=='s' | speckleArray.status(i)=='g'
            n=n+1;
            i=i+1;
        end
        if speckleArray.status(i)=='d'
            % End of the speckles reached
            
            % Store td
            td=i;
            
            if n==1
                % Ghost speckle
                gC=gC+1;
                stats.ghostDI(gC)=speckleArray.deltaI(tb+1);
                stats.numberOfGhost=stats.numberOfGhost+1;
            end
            if n>1
                nC=nC+1;            
                
                % Take all deltaIs
                jPos=0;
                for j=tb+1:td-1
                    jPos=jPos+1;
                    deltaIs(jPos)=speckleArray.deltaI(j);
                end
                stats.speckleDI(nC:nC+length(deltaIs)-1)=deltaIs;
                stats.numberOfSpeckle=stats.numberOfSpeckle+1;
                % Update pointer
                nC=nC+length(deltaIs)-1;
                
            end
            
            % Jump to the time point after the end of the speckle
            i=td+1;
            
        end
        
        if i<=length([speckleArray.timepoint])-1
            if speckleArray.status(i)=='b' 
                % Speckle with no 'd'
                i=i-1;
            end
        end
        
    else
        i=i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write to console
writeResults(stats,1);

if writeToFile==1

    % Open file 'log.txt'
	fid=fopen('log.txt','a');
	
	% Could the file be opened successfully?
	if fid==-1
		error('Couldn''t open the file.');
	end

    % Write to file
    writeResults(stats,fid);
    
    % Close file
	if fclose(fid)==-1;
		error('File  could not be closed!');
    end
    
end    

function writeResults(stats,fid)

fprintf(fid,'\n\nSummary from speckleArray (%s)\n-------------------------------------------------\n\n',datestr(now));
fprintf(fid,'Total number of speckles [s|g|f|l]  : %d\n',stats.sCounter+stats.gCounter);
fprintf(fid,'Number of events                    : %d\n',stats.events);
fprintf(fid,'Number of (complete) speckles       : %d\n',stats.complete);
fprintf(fid,'Number of speckles without birth    : %d\n',stats.dSpeckle);
fprintf(fid,'Number of speckles without death    : %d\n',stats.bSpeckle);
fprintf(fid,'Number of gaps closed               : %d\n',stats.gCounter);
fprintf(fid,'Number of ghost speckles            : %d\n',stats.numberOfGhost);
fprintf(fid,'COMPLETE speckles with lifetime more\nthan one frame (trajectories)       : %d\n',stats.complete-stats.numberOfGhost);
fprintf(fid,'Mean Life Time                      : %f\n',stats.meanLifeTime);
fprintf(fid,'Ghost speckles of all speckles      : %d%%\n',round(stats.numberOfGhost/stats.complete*100));
fprintf(fid,'Number of ''b'' or ''d'' events\ncorresponding to insignificant\nlocal maxima                        : %d/%d (%.2f%%)\n',stats.weakEvents,stats.events,100*stats.weakEvents/stats.events);
fprintf(fid,'\n');
fprintf(fid,'Average/median POLY score           : %f / %f\n',mean(stats.pScores),median(stats.pScores));
fprintf(fid,'Strongest POLY score                : %f\n',stats.polyScore);
fprintf(fid,'Average/median DEPOLY score         : %f / %f\n',mean(stats.dScores),median(stats.dScores));
fprintf(fid,'Strongest DEPOLY score              : %f\n',stats.depolyScore);
fprintf(fid,'\n');
fprintf(fid,'deltaI for long-living speckles (complete: %5d)\n',stats.numberOfSpeckle);
fprintf(fid,'-------------------------------------------------\n');
fprintf(fid,'Average/median deltaI               : %f / %f\n',mean(stats.speckleDI),median(stats.speckleDI));
fprintf(fid,'Standard deviation on deltaI        : %f\n',std(stats.speckleDI));
fprintf(fid,'Strongest deltaI                    : %f\n',max(stats.speckleDI));
fprintf(fid,'\n');
fprintf(fid,'deltaI for ghost speckles (complete: %5d)\n',stats.numberOfGhost);
fprintf(fid,'-------------------------------------------\n');
fprintf(fid,'Average/median deltaI               : %f / %f\n',mean(stats.ghostDI),median(stats.ghostDI));
fprintf(fid,'Standard deviation on deltaI        : %f\n',std(stats.ghostDI));
fprintf(fid,'Strongest deltaI                    : %f\n',max(stats.ghostDI));
fprintf(fid,'\n');
fprintf(fid,'Ratio f(deltaI(speckles))/f(deltaI(ghosts))\n');
fprintf(fid,'-------------------------------------------\n');
fprintf(fid,'f := mean                           : %f\n',mean(stats.speckleDI)/mean(stats.ghostDI));
fprintf(fid,'f := median                         : %f\n',median(stats.speckleDI)/median(stats.ghostDI));
fprintf(fid,'f := standard deviation             : %f\n',std(stats.speckleDI)/std(stats.ghostDI));
