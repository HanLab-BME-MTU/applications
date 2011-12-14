function labelgui2_makeMovie
%LABELGUI2_MAKEMOVIE creates movies from the labelgui2-window
%
% SYNOPSIS: labelgui2_makeMovie
%
% INPUT 
%
% OUTPUT 
%
% REMARKS
%
% SEE ALSO labelgui2
%
% EXAMPLE
%

% created with MATLAB ver.: 7.13.0.564 (R2011b) on Mac OS X  Version: 10.7.2 Build: 11C74 
%
% created by: Jonas Dorn
% DATE: 05-Dec-2011
%
% Last revision $Rev: 1973 $ $Date: 2011-05-31 18:50:22 -0400 (Tue, 31 May 2011) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, movieWindowHandles] = LG_getNaviHandles;

if isempty(movieWindowHandles)
    return
end

nTimepoints = movieWindowHandles.dataProperties.nTimepoints;

[fname,pname] = uigetfile('*.avi','select movie name');

fps = 3; 

% code lifted from showMovieGui - may contain some errors
if verLessThan('Matlab','7.11')
    
    aviObj = avifile(fullfile(pname,fname),'compression','none','fps',fps);
    % set 30 minutes per frame
    aviObj.fps = fps;
    aviObj.quality = 75; % no effect if no compression
    
    %for the moment, don't use wmv3.
    aviObj.compression = 'none';
    legacy = true;
else
    aviObj = VideoWriter(fullfile(pname,fname),'Uncompressed AVI')';
    aviObj.FrameRate = fps;
   % aviObj.Quality = 90;
    aviObj.open;
    legacy = false;
end


% loop to add frames
cData = [];
for t = 1:nTimepoints
    LG_gotoFrame(t);
    figure(movieWindowHandles.movieWindow)
    %fullscreen(fh)
    frame = getframe(movieWindowHandles.movieWindow);
    % make sure that all the frames have the same size
    if isempty(cData)
        cData = NaN(size(frame.cdata));
    else
        tmp = cData;
        [x,y,z] = size(frame.cdata);
        [xt,yt,zt] = size(tmp);
        tmp(1:min(xt,x),1:min(yt,y),1:min(zt,z)) = ...
            frame.cdata(1:min(xt,x),1:min(yt,y),1:min(zt,z));
        frame.cdata = uint8(tmp);
    end
    if legacy
        aviObj = addframe(aviObj,frame);
    else
        aviObj.writeVideo(frame);
    end
end
% finish up
if legacy
    aviObj=close(aviObj); %#ok<NASGU>
else
    aviObj.close;
end

fprintf('movie %s%s%s has been successfully created\n',pname,filesep,fname);
