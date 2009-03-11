function [results] = CorrelateData2Pos(positions, data, rdist, type)
% CorrelateData2Pos correlates data in variable format to the positions of 
% interest contained in the positions file in mpm format; rdist specified
% the distance around the positions
%
% SYNOPSIS [results] = CorrelateData2Pos(positions, data, rdist, type)
%
% INPUT     positions:  mpm file containing positions of tracked features 
%                       in x,y-columns
%           data:       variable input data; can be intensity image,
%                       parameter files, tracking data
%           rdist:      radius of interest around the positions
%           type:       type of analysis
%                       'density'
%                       'kscores'
%                       'cscores'
%                       'intensity'
%                       'motion'
%                       'correspondence'
%
% OUTPUT:   results:        
%
% last modified: Dinah Loerke, July 24, 2008



% check inputs
[px,py] = size(positions);
if mod(py,2)>0
    error('positions matrix has odd number of rows - needs to be MPM-type file');
end


switch type

    % for density, motion, and correspondence, the input needs to be a
    % tracking-type matrix in mpm format
    case {'density','motion','correspondence'}
        [dx,dy] = size(data);
        
        if mod(dy,2)>0
            disp('Analysis requires MPM-type data input');
            error('data input has odd number of rows');
        end

    % for kscores and cscores, the input needs to be a structure of files,
    % with the length corresponding to the number of frames
    case {'kscores','cscores'} 
        if ~isstruct(data) & ~iscell(data)
            disp('Analysis requires structure or cell array as data input');
            error('data input is not a structure or cell array');
        end
    % for intensity, the input needs to be an image
    case 'intensity'
        
        [ix,iy,iz] = size(data);
        if ( min(ix,iy)==1 ) 
            disp('Analysis requires an image (or image stack) as data input');
            error('data input does not have the correct format');
        end
        
    otherwise
      disp('Unknown method.')
end



switch type

    % for density, motion, and correspondence, the input needs to be a
    % tracking-type matrix in mpm format
    case 'density'
        [results] = CorrelateData2Pos_density(positions, data, rdist);
    case 'motion'
        [results] = CorrelateData2Pos_motion(positions, data, rdist);
    case 'correspondence'
        [results] = CorrelateData2Pos_correspondence(positions, data, rdist);
    case 'kscores'
        [results] = CorrelateData2Pos_kscores(positions, data, rdist,3);
    case 'cscores'
        [results] = CorrelateData2Pos_cscores(positions, data, rdist);
    case 'intensity'
        [results] = CorrelateData2Pos_intensity(positions, data, rdist);          
    otherwise
      disp('method not yet implemented')
end

% display module

end % of function
