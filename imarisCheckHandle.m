function imarisApplication = imarisCheckHandle(handleOrIndex)
%IMARISCHECKHANDLE returns the imarisApplication handle for a program that has to be called with an active Imaris
%
% SYNOPSIS imarisApplication = imarisCheckHandle(handleOrIndex)
%
% INPUT    handleOrIndex: Either a handle to an imarisApplication, or the
%                         index to the current Imaris as communicated to
%                         the server. 
%
% OUTPUT   imarisApplication : handle to the current Imaris
%
% REMARKS  if the input is neither handle nor index, the program throws an
%          error. Use it in a try...catch loop if you want to start a new
%          session in that case.
%
% c: jonas 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===============
% TEST INPUT
%===============

if nargin == 0 || isempty(handleOrIndex)
    error('Please specify an input argument for IMARISCHECKHANDLE!')
end

%===============

%===================
% FIND APPLICATION
%===================

% check for handle
type = whos('handleOrIndex');

if strcmpi(type.class, 'COM.Imaris_Application')
    % we have an application
    imarisApplication = handleOrIndex;
else
    % start the server
    imaServer = actxserver('Imaris.Server');
    % check the id
    imarisApplication = imaServer.GetObject(handleOrIndex);
    
    % make sure the application is good
    type = whos('imarisApplication');
    
    if strcmpi(type.class, 'COM.Imaris_Application')
        % all good
    else
        % we have a problem
        error('IMARISCHECKHANDLE: No valid index or handle!')
    end
end

% if we came this far, we have a handle. Make it visible
imarisApplication.mVisible = 1;