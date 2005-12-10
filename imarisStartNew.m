function imarisApplication = imarisStartNew
%IMARISSTARTNEW starts a new Imaris and stores the handle in the workspace
% naming it imarisApplication, or imarisApplication#, where # is a number
%
% SYNOPSIS imarisApplication = imarisStartNew
%
% OUTPUT   imarisApplication : Handle to the new imaris session
%
%c: jonas 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no input to test, therefore start imaris

imarisApplication = actxserver('Imaris.Application');
imaAppName = 'imarisApplication';
% make sure we do not accidentially overwrite an imaris application
if evalin('base','exist(''imarisApplication'',''var'')')
    num = 2;
    while evalin('base',['exist(''imarisApplication',num2str(num),''',''var'')'])
        num = num+1;
    end
    imaAppName = [imaAppName, num2str(num)];
end
imarisApplication.mVisible = 1;
assignin('base',imaAppName,imarisApplication);
disp(sprintf(['The handle to the current imaris is ''%s''\n',...
    'Delete it to close this session of Imaris'],imaAppName));

