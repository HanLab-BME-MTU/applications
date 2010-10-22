function [channelPlugins layerPlugins] = getPlugins()

%
% CHANNELS 
%

channelPlugins = struct('desc', {},  'filterSpec', {}, 'loadFunc', {});

% Raw image
channelPlugins(1).desc = 'Raw Image';
channelPlugins(1).filterSpec = {'*.tif';'*jpg';'*.png'};
channelPlugins(1).loadFunc = @imread;

% Scalar Map
channelPlugins(2).desc = 'Scalar Map';
channelPlugins(2).filterSpec = {'*.mat'};
channelPlugins(2).loadFunc = @loadScalarMap;

% [qFSM] Poly Map
channelPlugins(3).desc = '[qFSM] Poly Map';
channelPlugins(3).filterSpec = {'*.mat'};
channelPlugins(3).loadFunc = @loadScalarMap;

% [qFSM] Depoly Map
channelPlugins(4).desc = '[qFSM] Depoly Map';
channelPlugins(4).filterSpec = {'*.mat'};
channelPlugins(4).loadFunc = @loadDepolyMap;

%
% LAYERS
%

layerPlugins = struct('desc', {}, 'filterSpec', {}, 'display', {});

% [qFSM] Speckles
layerPlugins(1).desc = '[qFSM] Speckles';
layerPlugins(1).filterSpec = {'*.mat'}; 
layerPlugins(1).displayFunc = @displaySpeckles;

% [PANDA] Window
layerPlugins(2).desc = '[PANDA] Windows';
layerPlugins(2).filterSpec = {'*.mat'}; 
layerPlugins(2).displayFunc = @plotWindowsFSM;

end