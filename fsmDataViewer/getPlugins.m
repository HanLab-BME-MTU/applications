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

% Vector Field
layerPlugins(1).desc = 'Vector Field';
layerPlugins(1).filterSpec = {'*.mat'};
layerPlugins(1).displayFunc = @displayVectorField;

% Points
layerPlugins(2).desc = 'Points';
layerPlugins(2).filterSpec = {'*.mat'};
layerPlugins(2).displayFunc = @displayPoints;

% [qFSM] Speckles
layerPlugins(3).desc = '[qFSM] Speckles';
layerPlugins(3).filterSpec = {'*.mat'}; 
layerPlugins(3).displayFunc = @displaySpeckles;

% [IF] Init Curves
layerPlugins(4).desc = '[IF] Init Curves';
layerPlugins(4).filterSpec = {'*.mat'}; 
layerPlugins(4).displayFunc = @displayInitCurves;

% [IF] Curve Models
layerPlugins(5).desc = '[IF] Curve Models';
layerPlugins(5).filterSpec = {'*.mat'}; 
layerPlugins(5).displayFunc = @displayCurveModels;

% [PANDA] Window
layerPlugins(6).desc = '[PANDA] Windows';
layerPlugins(6).filterSpec = {'*.mat'}; 
layerPlugins(6).displayFunc = @plotWindowsFSM;

end