function [channelPlugins layerPlugins] = getPlugins()

%
% CHANNELS 
%

channelPlugins = struct('desc', {},  'filterSpec', {}, 'load', {});

% Raw images
channelPlugins(1).desc = 'Raw Images';
channelPlugins(1).filterSpec = {'*.tif';'*jpg';'*.png'};
channelPlugins(1).load = @imread;

% [qFSM] Speed Map
channelPlugins(2).desc = '[qFSM] Speed Map';
channelPlugins(2).filterSpec = {'*.mat'};
channelPlugins(2).load = @loadSpeedMap;

% [qFSM] Poly Map
channelPlugins(3).desc = '[qFSM] Poly Map';
channelPlugins(3).filterSpec = {'*.mat'};
channelPlugins(3).load = @loadPolyMap;

% [qFSM] Depoly Map
channelPlugins(4).desc = '[qFSM] Depoly Map';
channelPlugins(4).filterSpec = {'*.mat'};
channelPlugins(4).load = @loadDepolyMap;

%
% LAYERS
%

layerPlugins = struct('desc', {}, 'filterSpec', {}, 'display', {});

% [qFSM] Speckles
layerPlugins(1).desc = '[qFSM] Speckles';
layerPlugins(1).filterSpec = {'*.mat'}; 
layerPlugins(1).display = @displaySpeckles;

% [Filter] Orientation
layerPlugins(2).desc = '[Filter] Orientation';
layerPlugins(2).filterSpec = {'*mat'};
layerPlugins(2).display = @displayOrientation;

end