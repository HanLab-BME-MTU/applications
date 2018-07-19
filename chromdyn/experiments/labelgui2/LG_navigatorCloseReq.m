function LG_navigatorCloseReq
%LG_navigatorCloseReq is the closerequest function for the navigator window

% if loadedFromOutside, we still close. This is sort of equivalent to
% aborting the job.

% loop through handles in naviHandles.movieWindowH and successively try to
% close.


% get handles
naviHandles = LG_getNaviHandles;

% loop
nWindows = length(naviHandles.movieWindowH);
status = ones(nWindows,1);
for i=1:nWindows
    if ishandle(naviHandles.movieWindowH(i))
        status(i) = close(naviHandles.movieWindowH(i));
    else
        naviHandles.movieWindowH(i) = 0;
        status(i) = 1;
    end
end

if isempty(status) || all(status == 1)
    % empty status: there are no windows at all
    delete(naviHandles.LG_navigator);
else
    % don't close navigatorWindow
end