function [experiment] = loadConditionDataDualChannel()
% loadConditionData loads the relevant information for all the data
% available for a specific experiment condition; this requires a specific
% dircetory structure and nomenclature (see below) 
%
% SYNOPSIS [experiment] = loadConditionData()
%
% INPUT    
%
% OUTPUT   experiment: structure with the fields
%                       .source = pathname of the data/movie
%                       .date = date when movie was taken
%                       .framerate = framerate of the movie (2s or 0.4s)
%
% 
%
% Dinah Loerke, January 24th, 2008


% select directory where all data for this condition are located

directory_or = cd;
[directory_cond] = uigetdir(directory_or,'Please select the folder for this condition'); 
cd(directory_cond);

% get files in specified directory
files_cond = dir(directory_cond);
% files(i).name contains the string name of the file
% files(i).isdir contains the true/false information about whether the file


% identify the usable directories that aren't '.' or '..'
goodlist_cond = zeros(length(files_cond),1);
for g=1:length(files_cond)
    if ( length(files_cond(g).name)>2 ) & ( files_cond(g).isdir == true )
        goodlist_cond(g) = 1;
    end
end
usedir_cond = find(goodlist_cond);

% initialize counter
ct = 1;

if length(usedir_cond)>0
    
    % loop over approriate directories
    for i=1:length(usedir_cond)
        
        % current folder name
        dirname_date = files_cond(usedir_cond(i)).name;
        % search for date numbers, either at right or left margin of the
        % folder name
        idxs=length(dirname_date);
        % if last (or first) char is a number
        if (uint8(dirname_date(idxs))>47 & uint8(dirname_date(idxs))<58)
            idc = idxs;
            while(uint8(dirname_date(idc))>47 & uint8(dirname_date(idc))<58)
                idc = idc-1;
            end
            currDate = dirname_date(idc+1:idxs);
            
        elseif (uint8(dirname_date(1))>47 & uint8(dirname_date(1))<58)
            idc = 1;
            while(uint8(dirname_date(idc))>47 & uint8(dirname_date(idc))<58)
                idc = idc+1;
            end
            currDate = dirname_date(1:idc-1);
        
        else
            % default for date = 010101
            currDate = '010101';
        end
        
        
        % move to current date directory
        cd(dirname_date);
        directory_date = cd;

        % look for the individual cell data in this folder
        files_cell = dir(cd);

        
        % identify the directories that aren't '.' or '..'
        goodlist_cell = zeros(length(files_cell),1);
        for g=1:length(files_cell)
            if ( length(files_cell(g).name)>2 ) & ( files_cell(g).isdir == true )
                goodlist_cell(g) = 1;
            end
        end
        usedir_cell = find(goodlist_cell);

        % loop over all cell folders in this directory

        if length(usedir_cell)>0
            
            for k=1:length(usedir_cell)

                % current cell directory
                dirname_cell = files_cell(usedir_cell(k)).name;

                % extract framerate from the name of the cell folder
                % NOTE: if there's no specific identification for fast, then
                % the default is slow

                inFast1 = findstr(dirname_cell, 'fast');
                inFast2 = findstr(dirname_cell, '400ms');

                if ( (length(inFast1)>0) | (length(inFast2)>0) )
                    currFramerate = 0.4;
                else 
                    currFramerate = 2;
                end
                
                
                % look for the individual channels in this cell folder
                cd(dirname_cell);
                files_channel = dir(cd);
                
                if k==1
                    dirname_channel1curr = uigetdir(cd,'select first (master) channel (e.g. CCP channel)');
                    dirname_channel2curr = uigetdir(cd,'select second (slave) channel (e.g. other protein)');
                    % separate final folder name
                    k1 = findstr(dirname_channel1curr, filesep);
                    k2 = findstr(dirname_channel2curr, filesep);
                    p1 = max(k1); l1 = length(dirname_channel1curr);
                    p2 = max(k2); l2 = length(dirname_channel2curr);
                    dirname_channel1root = dirname_channel1curr(p1+1:l1);
                    dirname_channel2root = dirname_channel2curr(p2+1:l2);
                else
                    
                    if exist(dirname_channel1root)==7
                        dirname_channel1curr = [cd,filesep,dirname_channel1root];                        
                    else
                        dirname_channel1curr = uigetdir(cd,'select first (master) channel (e.g. CCP channel)');                        
                    end
                    
                    if exist(dirname_channel2root)==7
                        dirname_channel2curr = [cd,filesep,dirname_channel2root];
                    else
                        dirname_channel2curr = uigetdir(cd,'select second (slave) channel (e.g. other protein)');
                    end
                    
                end
                

                % enter data
                experiment(ct).source = dirname_channel1curr;
                experiment(ct).channel1 = dirname_channel1curr;
                experiment(ct).channel2 = dirname_channel2curr;
                experiment(ct).date = currDate;
                experiment(ct).framerate = currFramerate;

                ct = ct+1;

                cd(directory_date);
            end % of for k-loop
            
        end % of if length(usedir_cell)>0
        
        cd(directory_cond);
        
    end % of for i-loop

else
    error('no usable data in directory');
end % of if

end % of function
            