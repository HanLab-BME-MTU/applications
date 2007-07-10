function NA = naFromLensID(lensID,ask)
%NAFROMLENSID is a lookup table to find the NA for a given lensID
%
% SYNOPSIS  NA = naFromLensID(lensID,ask)
%
% INPUT     lensID : ID of the lens according to deltaVision
%           ask    : Whether or not to ask for lensID if not found in list
%                    Default: 1
%
% OUTPUT    NA for the specified lens
%
% REMARKS   After the program asks for a lensID, it will write it into the
%             .m-file. Alternatively, you can write it into the list at the
%             end of the code yourself.
%
% c: 8/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TEST INPUT
if nargin == 0
    error('not enough input arguments')
end
if nargin < 2 || isempty(ask)
    ask = 1;
end
NA = [];

% LOOKUP LENSID
lensIDList = lookupLensIDs;

lensIdx = find(lensIDList(:,1) == lensID);

% ASSIGN OR ASK
if ~isempty(lensIdx)
    NA = lensIDList(lensIdx(1),2);
elseif ask
    ans = inputdlg(...
        {sprintf(['Please input NA for lens %i'],lensID)},...
        'Manual NA');
    if ~isempty(ans)
        NA = str2double(ans{1});
    end
    if isempty(NA) || ~isnumeric(NA)
        NA = [];
    else
        % add the lensID to this file
        % read naFromLensID.m
        file = which('naFromLensID.m');
        mFile = textread(file,'%s',...
            'delimiter','\n','whitespace','');
        % now add a line just befor the last line
        newLine = sprintf('     %i,  %1.2f;...',lensID,NA);
        mFile{end+1} = mFile{end};
        mFile{end-1} = newLine;
        
        % finally, write the mFile
        fid = fopen(file,'w');
        for i=1:length(mFile)
            fprintf(fid,'%s\n',mFile{i});
        end
        fclose(fid);
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lensID = lookupLensIDs
% add lensID's and NA's by adding a line like:
% '9999,  1.25;...' into the table below. Do not add any lines after the
% closing of the table '];'

lensID = ...
    [10002,  1.4;...
     10006,  1.4;...
     10612,  1.4;...
     12003,  1.4;...
     0,  1.40;...
     10005,  1.35;...
     10105,  1.40;...
     10007,  1.40;...
     ];
