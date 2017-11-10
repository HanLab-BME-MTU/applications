function status=openDVExperiment(fullPathToFile,expID)
% OPENDVEXPERIMENT loads DV library to access DV experiment files for reading
%
% Input
%    fullPathToFile: string, complete path to file to be opened
%             expID: int, experiment ID for later assignment/access
% Output
%            status: int, file could be opened successfully or not
%
% US, 01/02/2012
%

if( expID < 0 )
    error('In openDVExperiment: expID must be larger than zero!');
end

%if( fullPathToFile(1) ~= '/' )
%    error('In openDVExperiment: please provide complete path to DV file!');
%end

% open DV Image I/O library
% DVImgLibOpen(0);
% check supplied file for errors
err=DVImgCheckFile(fullPathToFile);
if( err ~= 0 )
    error('In openDVExperiment: identified corrupt file!');
end
% open supplied file
status=DVImgOpen(expID,fullPathToFile,'ro');