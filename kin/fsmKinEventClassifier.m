function speckleArray=fsmKinEventClassifier(speckleArray,fsmParam)
%
% SYNOPSIS   speckleArray=fsmKinEventClassifier(speckleArray,fsmParam)
%
% INPUT      speckleArray : structure containing all speckle information from a movie
%            fsmParam     : general parameter structure
%
% OUTPUT     speckleArray : speckleArray with added kinetic information

% DEBUGGING
displayInfo=0;

if nargin~=2
    error('Two parameters expected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ NEEDED PARAMETERS FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timepoints=3; % Forced to be 3 (3-point fitting)
time=fsmParam.specific.imageNumber;
k=fsmParam.main.noiseParam(5);
bleachRed=fsmParam.kin.bleachRed;

% Current directory
oldDir=cd;


if length(speckleArray)==1;
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% START PROCESSING SPECKLES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=length(speckleArray);            % End
c1=0;                              % Start

numberOfEvents=0;
numberOfSignificantEvents=0;
numberOfInsignificantEvents=0;
numberOfBleachEvents=0;

% Initializing progress bar
h = waitbar(0,'Classifying speckles');

while c1<=(l-1)
	
	% Reset variables
	intensity=[]; background=[]; sigmaMax=[]; sigmaMin=[]; sigmaDiff=[]; event=[];
	
	% Update counter
	c1=c1+1;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Check speckle
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if speckleArray(c1).status=='b'
		
		numberOfEvents=numberOfEvents+1;
		
		% COLLECT DATA TO HANDLE A BIRTH EVENT
		
		birthTime=speckleArray(c1+1).timepoint;
		
		event='birth';
		
		if c1<=(l-timepoints+1)
			
			intensity(1)=speckleArray(c1).intensity;
			background(1)=speckleArray(c1).background;
			sigmaMax(1)=speckleArray(c1).sigmaSp;
			sigmaMin(1)=speckleArray(c1).sigmaBg;
			sigmaDiff(1)=sqrt(sigmaMax(1)^2+sigmaMin(1)^2);
			
			for c2=1:timepoints-1
				
				if (speckleArray(c1+c2).status=='s' | speckleArray(c1+c2).status=='g' | speckleArray(c1+c2).status=='l') & speckleArray(c1+c2).timepoint==speckleArray(c1).timepoint+c2
					intensity(1+c2)=speckleArray(c1+c2).intensity;
					background(1+c2)=speckleArray(c1+c2).background;
					sigmaMax(1+c2)=speckleArray(c1+c2).sigmaSp;
					sigmaMin(1+c2)=speckleArray(c1+c2).sigmaBg;
					sigmaDiff(1+c2)=sqrt(sigmaMax(1+c2)^2+sigmaMin(1+c2)^2);
					
				else
					
					birthTime=[];
					intensity=[];
					event=[];
					break;
					
				end
				
			end
			
		end
		
	end 
	
	if speckleArray(c1).status=='d' 
		
		numberOfEvents=numberOfEvents+1;

		% COLLECT DATA TO HANDLE A DEATH EVENT
		deathTime=speckleArray(c1).timepoint-1;
		
		event='death';
		
		if c1>=timepoints
			
			intensity(timepoints)=speckleArray(c1).intensity;
			background(timepoints)=speckleArray(c1).background;
			sigmaMax(timepoints)=speckleArray(c1).sigmaSp;
			sigmaMin(timepoints)=speckleArray(c1).sigmaBg;
			sigmaDiff(timepoints)=sqrt(sigmaMax(timepoints)^2+sigmaMin(timepoints)^2);
			
			for c2=1:timepoints-1
				
				if (speckleArray(c1-c2).status=='s' | speckleArray(c1-c2).status=='g' | speckleArray(c1-c2).status=='f') & speckleArray(c1-c2).timepoint==speckleArray(c1).timepoint-c2
					intensity(timepoints-c2)=speckleArray(c1-c2).intensity;
					background(timepoints-c2)=speckleArray(c1-c2).background;
					sigmaMax(timepoints-c2)=speckleArray(c1-c2).sigmaSp;
					sigmaMin(timepoints-c2)=speckleArray(c1-c2).sigmaBg;
					sigmaDiff(timepoints-c2)=sqrt(sigmaMax(timepoints-c2)^2+sigmaMin(timepoints-c2)^2);
					
				else
					
					deathTime=[];
					intensity=[];
					event=[];
					
					break;
					
				end
				
			end
			
		end
		
	end 
	
	if speckleArray(c1).status=='s' | speckleArray(c1).status=='g' | speckleArray(c1).status=='f' | speckleArray(c1).status=='l'
		
		% We are inside a speckle - nothing to classify
		intensity=[];
		event=[];
		
	end
	
	if ~isempty(intensity) % An event has been found
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%
		% NOW PERFORM ANALYSIS ON CURRENT (SELECTED) SPECKLE
		%
		% LEAST-SQUARE WITH KNOWN COVARIANCE: LSCOV
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%
		% (1) Check if it must be a 2- or 3-point fitting
		%
		
		%&*&%                                                                                           %&*&%
		%&*&%     PROBLEM: THE SYSTEM MUST BE OVERDETERMINED FOR LSCOV - THEREFORE NO 2x2 A ALLOWED     %&*&%  
		%&*&%                                                                                           %&*&%
		
		deltaI=intensity-background;   % The slope of deltaI will be the score
		
		%
		% (2) Fit intensity, background and deltaI
		%
		%    Classification of events:
		%    The value slope(deltaI)/sigma(slope(deltaI)) decides if the event is statistically significant
		%    The value slope(speckle)/sigma(slope(speckle)) and slope(background)/sigma(slope(background))
		%       decide whether speckle slope and background slopes are significant, respectively.
		%       If yes: critS=slope(speckle) and/or critB=slope(background)
		%       If no:  critS=0, and/or critB=0;
				
		% Time
		t=[1:timepoints]';
		% A
		A=[t,ones(timepoints,1)];
		
        % DELTAI
		K=deltaI(1:timepoints)';
        X=diag(sigmaDiff(1:timepoints).^2);
        [D,dD]=lscov(A,K,X);                          % Slope D(1) and error dD(1) on deltaI
        if dD(1)==0 | ~isreal(dD(1))
            relD=0; % Make it nonsignificant -> the event will not be classified
            if displayInfo==1
                clc;
                fprintf(1,'LSCOV returned complex numbers. The event will not be classified.\n');
                pause;
            end
        else
            relD=D(1)/dD(1);
        end
        		
		% If relD is not significant, directly jump to the next event
		if abs(relD)<k
			numberOfInsignificantEvents=numberOfInsignificantEvents+1;
            continue;
        else	
        
            % INTENSITY
            I=intensity(1:timepoints)';
            V=diag(sigmaMax(1:timepoints).^2);
            [S,dS]=lscov(A,I,V);                          % Slope S(1) and error dS(1) on intensity
            if dS(1)==0 | ~isreal(dS(1))
                relS=0;		                              % This check is needed because sometimes LSCOV returns complex numbers in case of zero slopes
                critS=0;                                  % Set the slope for comparison to zero
            else
                relS=S(1)/dS(1);
                if abs(relS)>=k                        % The event is statistically significant
                    critS=S(1);                          % The slope can be used for comparison
                else
                    critS=0;
                end
            end
            
            % BACKGROUND
            J=background(1:timepoints)';
            W=diag(sigmaMin(1:timepoints).^2);
            [B,dB]=lscov(A,J,W);                          % Slope B(1) and error dB(1) on background
            if dB(1)==0 | ~isreal(dB(1))
                relB=0;		                              % This check is needed because sometimes LSCOV returns complex numbers in case of zero slopes
                critB=0;                                  % Set the slope for comparison to zero			
            else
                relB=B(1)/dB(1);
                if abs(relB)>=k                        % The event is statistically significant
                    critB=B(1);                          % The slope can be used for comparison
                else
                    critB=0;                               % Set the slope for comparison to zero
                end
            end
            
            
            % Analyze this
            numberOfSignificantEvents=numberOfSignificantEvents+1;						
			
            %
			% (3) Assign the event to a class 
			%

			if displayInfo==1         
				plot([1:timepoints],intensity,'r*');
				title(event);
				hold on;
				plot([1:timepoints],S(1)*[1:timepoints]+S(2),'k-');
				plot([1:timepoints],background,'b*');
				plot([1:timepoints],B(1)*[1:timepoints]+B(2),'k-');
				hold off;
				
				clc;
				switch event
				case 'birth', fprintf(1,'Birth time: %d\n',birthTime);
				case 'death', fprintf(1,'Death time: %d\n',deathTime);
				otherwise
				end
				if critS~=0
					signS='significant';
				else 
					signS='non significant';
				end
				fprintf(1,'Speckle    [  red ]: slope %6f relS %8f [%s]\n',S(1),relS,signS);
				if critB~=0
					signB='significant';
				else 
					signB='non significant';
				end
				fprintf(1,'Background [ blue ]: slope %6f relB %8f [%s]\n',B(1),relB,signB);
				switch abs(relD)>k
				case 0
					valid='non significant';
				case 1
					valid='significant';
				otherwise
					valid='check!!!';
				end
				fprintf(1,'deltaI             : slope %6f relD %8f [%s]\n',D(1),abs(relD),valid);
			end
			
			
			
			% Reset
			polScore=0; depolScore=0; balance=0;
			
			switch event
				
			case 'birth'
				
				switch sign(critS) %relS
					
				case -1
					
					%%%
					
					switch sign(critB) %relB
						
					case -1
						
						if abs(critS)<abs(critB)                 % BIRTH DUE TO BACKGROUND DISSOCIATION
							
							polScore=0;
							depolScore=1;
							balance=B(1); % Score is the negative background slope
							
						end
						
						if abs(critS)==abs(critB)                % PATHOLOGICAL CASE -> IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
						if abs(critS)>abs(critB)                 % PATHOLOGICAL CASE ->IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
					case 0                                     % PATHOLOGICAL CASE -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					case 1                                     % PATHOLOGICAL CASE -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					otherwise
					end
					
					%%%
					
				case 0                                        % The SPECKLE SLOPE is 0
					
					%%%
					
					switch sign(critB)
						
					case -1                                   % Slopes are diverging (D(1)>0), background slope negative
						                                          % SPECKLE BIRTH DUE TO BACKGROUND DEPOLYMERIZATION
						polScore=0;
						depolScore=1;
						balance=B(1);  % The score is the negative background slope
						
					case 0                                    % Both slopes are 0
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					case 1			                          % PATHOLOGICAL CASE -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					otherwise
					end
					
					%%%
					
				case 1                                        % The SPECKLE SLOPE is positive
					
					%%%
					
					switch sign(critB)
						
					case -1                                   % Slopes are diverging (D(1)>0), backgorund slope negative
						
						if critS>abs(critB)                     % BIRTH DUE TO SPECKLE POLYMERIZATION
							
							polScore=1;
							depolScore=0;
							balance=S(1); % The score is the positive speckle slope
							
						end
						
						if critS==abs(critB)                    % Arbitrarily assign score to polymerization
							
							polScore=1;
							depolScore=0;
							balance=S(1); % The score is the positive speckle slope
							
						end
						
						if critS<abs(critB)                     % BIRTH DUE TO BACKGROUND DEPOLYMERIZATION
							
							polScore=0;
							depolScore=1;
							balance=B(1); % The score is the negative background slope
							
						end
						
					case 0                                    % BIRTH DUE TO SPECKLE POLYMERIZATION
						
						polScore=1;
						depolScore=0;
						balance=S(1); % The score is the positive speckle slope
						
					case 1			                          % Both slopes are POSITIVE
						
						if critS>critB                          % BIRTH DUE TO SPECKLE POLYMERIZATION
							
							polScore=1;
							depolScore=0;
							balance=S(1); % The score is the positive speckle slope
							
						end
						
						if critS==critB                         % PATHOLOGICAL CASE -> IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
						if critS<critB                          % PATHOLOGICAL CASE -> IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
					otherwise	
					end
					
					%%%
					
				otherwise
				end
				
				%%%
				
			case 'death'
				
				switch sign(critS)
					
				case -1
					
					%%%
					
					switch sign(critB)
						
					case -1                                    
						
						if abs(critS)<abs(critB)                 % PATHOLOGICAL CASE -> IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
						if abs(critS)==abs(critB)                % PATHOLOGICAL CASE ->IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
						if abs(critS)>abs(critB)                 % DEATH DUE TO SPECKLE DEPOLYMERIZATION
							
							polScore=0;
							depolScore=1;
							balance=S(1); % The score is the negative speckle slope
							
						end
						
					case 0                                     % DEATH DUE TO SPECKLE DEPOLYMERIZATION
						
						polScore=0;
						depolScore=1;
						balance=S(1); % The score is the negative speckle slope
						
					case 1
						
						if abs(critS)>critB                 % DEATH DUE TO SPECKLE DEPOLYMERIZATION
							
							polScore=0;
							depolScore=1;
							balance=S(1); % The score is the negative speckle slope
							
						end
						
						if abs(critS)==critB                % Arbitrarily assign score to speckle depolymerization
							
							polScore=0;
							depolScore=1; 
							balance=S(1); % The score is the negative speckle slope
							
						end
						
						if abs(critS)<critB                 % DEATH DUE TO BACKGROUND POLYMERIZATION
							
							polScore=1;
							depolScore=0;
							balance=B(1); % Score is the positive background slope
							
						end
						
					otherwise
					end
					
					%%%
					
				case 0
					
					switch sign(critB)
						
					case -1                                    % PATHOLOGICAL CASE -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					case 0                                     % Pathological case -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					case 1                                     % DEATH DUE TO BACKGROUND POLYMERIZATION
						
						polScore=1;
						depolScore=0;
						balance=B(1); % Score is the positive background slope
						
					otherwise
					end
					
					%%%
					
				case 1
					
					switch sign(critB)
						
					case -1                                    % PATHOLOGICAL CASE -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					case 0                                     % PATHOLOGICAL CASE -> IGNORE IT
						
						polScore=0;
						depolScore=0;
						balance=0;
						
					case 1
						
						if critS>critB                 % PATHOLOGICAL CASE -> IGNORE IT
							
							polScore=0;
							depolScore=0;
							balance=0;
							
						end
						
						if critS==critB                % PATHOLOGICAL CASE -> IGNORE IT
							
							polScore=0;
							depolScore=0; 
							balance=0;
							
						end
						
						if critS<critB                 % DEATH DUE TO BACKGROUND POLYMERIZATION
							
							polScore=1;
							depolScore=0;
							balance=B(1); % The score is the positive background slope
							
						end
						
					otherwise 
					end
					
					%%%
					
				otherwise
				end
				
				%%%
				
			otherwise
			end
			
			% LAST CHECK: BLEACH
			if abs(balance)<bleachRed
			
				% WRITE RESULTS FOR A BLEACH-SCORE
				speckleArray(c1).score=0;
				speckleArray(c1).activity=0;
				numberOfBleachEvents=numberOfBleachEvents+1;
				
			else
				
				% WRITE RESULTS
				speckleArray(c1).score=balance;
				if polScore==1 & depolScore==0
					speckleArray(c1).activity=1;
				elseif polScore==0 & depolScore==1
					speckleArray(c1).activity=-1;
				end
			end
			
			% Check
			if polScore==1 & depolScore==1
				error('A score can be due EITHER to polymerization OR to depolymerization');
			end			
			
			if displayInfo==1
				if polScore==1 
					fprintf('SCORED AS POLYMERIZATION (activity = %d; score = %f)\n',polScore,balance);
				elseif depolScore==1
					fprintf('SCORED AS DEPOLYMERIZATION (activity = %d; score = %f)\n',-depolScore,balance);
				else
				fprintf(1,'NO SCORE/ACTIVITY ASSIGNED TO THIS EVENT (activity = %d; score = %f)\n',polScore,balance);	
				end					
				pause;
			end
			
		end
		
	end
	
	% Update waitbar
	waitbar(c1/(l-1),h);

end

% Close waitbar
close(h);

% Write info to console and 'settings.log'
fid=fopen('settings.log','a');

% Could the file be opened successfully?
if fid==-1
	error('Couldn''t open the file.');
end

% To console
fprintf(1,'\n\nSpeckle classification - Summary\n--------------------------------\n\n');
fprintf(1,'[1] Total number of events                         : %d\n',numberOfEvents);
fprintf(1,'[2] Number of VALID events                         : %d\n',numberOfSignificantEvents+numberOfInsignificantEvents);
fprintf(1,'[3] Number of STATISTICALLY SIGNIFICANT events     : %d\n',numberOfSignificantEvents);
fprintf(1,'[4] Number of NON STATISTICALLY SIGNIFICANT events : %d\n',numberOfInsignificantEvents);
fprintf(1,'[5] Number of BLEACH events (a subset of [3])      : %d\n',numberOfBleachEvents);
fprintf(1,'[6] Number of RESULTING SCORES                     : %d\n\n\n',numberOfSignificantEvents-numberOfBleachEvents);

% To file
fprintf(fid,'#BEGIN\n\n');
fprintf(fid,'\n\nSpeckle classification (%s)\n',datestr(now));
fprintf(fid,'Summary\n-----------------------------------------------\n\n');
fprintf(fid,'[1] Total number of events                         : %d\n',numberOfEvents);
fprintf(fid,'[2] Number of VALID events                         : %d\n',numberOfSignificantEvents+numberOfInsignificantEvents);
fprintf(fid,'[3] Number of STATISTICALLY SIGNIFICANT events     : %d\n',numberOfSignificantEvents);
fprintf(fid,'[4] Number of NON STATISTICALLY SIGNIFICANT events : %d\n',numberOfInsignificantEvents);
fprintf(fid,'[5] Number of BLEACH events (a subset of [3])      : %d\n',numberOfBleachEvents);
fprintf(fid,'[6] Number of RESULTING SCORES                     : %d\n\n\n',numberOfSignificantEvents-numberOfBleachEvents);
fprintf(fid,'#END\n\n');

% Close file
if fclose(fid)==-1;
	error('File  could not be closed!');
end



