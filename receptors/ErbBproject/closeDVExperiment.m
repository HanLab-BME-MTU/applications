function status=closeDVExperiment(expID)
% CLOSEDVEXPERIMENT closes DV experiments and the DV I/O library
%
% Input
%     expID: int, experiment ID of file to be closed
% Output
%     status: int, file could be closed successfully or not
%
% US, 01/02/2012
%

if( expID < 0 )
    error('In closeDVExperiment: expID must be larger than zero!');
end

% close supplied file
status=DVImgClose(expID);
% close DV Image I/O library
DVImgLibClose();