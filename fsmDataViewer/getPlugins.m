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

% 'displayFunc' function displays layer information at a given frame into
% an axes. The prototype of 'displayFunc' should follows:
%
% @displayFunc(hAxes, tag, layer, layerColor) 
%
% There are 2 types of layers that can be used in this viewer:
% - 1 bundle file containing layer information for the whole movie
% - multiple file containing layer information per frame, 1 file / frame.
%
% 'dispatchFunc' function is responsible for dispatching data per frame,
% whether it comes from 1 bundle file or multiple files. Builins functions
% are provided to handle multiple senarios:
%
% @dispatchFilesToFrames:    1 file per frame. Each file should contains
%                            only 1 variable.
%
% @dispatchStructToFrames:   1 file containing a structure array where the
%                            number of elements is equal to the number of
%                            frames. The bundle file should contains only 1
%                            struct variable.
%
% @dispatchCellToFrames:     1 file containing a cell array where the
%                            number of elements is equal to the number of
%                            frames. The bundle file should contains only 1
%                            cell variable.
%
% @dispatchMatrix1ToFrames:  1 file containaing a matrix where the 1st
%                            dimension corresponds to the data of each
%                            frame. The bundle file should contains only 1
%                            matrix.
%
% @dispatchMatrix2ToFrames:  1 file containing a matrix where the 2nd
%                            dimension corresponds to the data of each
%                            frame. The bundle file should contains only 1
%                            matrix.
%
% @dispatchMatrix3ToFrames:  1 file containing a matrix where the 3rd
%                            dimension corresponds to the data of each
%                            frame. The bundle file should contains only 1
%                            matrix.

layerPlugins = struct(...
    'desc', {},...
    'filterSpec', {},...
    'display', {},...
    'inBundleFile', {});

% Vector Field
layerPlugins(1).desc = 'Vector Field';
layerPlugins(1).filterSpec = {'*.mat'};
layerPlugins(1).displayFunc = @displayVectorField;
layerPlugins(1).dispatchFunc = @dispatchCellToFrames;

% Iso Features
layerPlugins(2).desc = 'Iso Features';
layerPlugins(2).filterSpec = {'*.mat'};
layerPlugins(2).displayFunc = @displayIsoFeatures;
layerPlugins(2).dispatchFunc = @dispatchCellToFrames;

% Aniso Features
layerPlugins(3).desc = 'Aniso Features';
layerPlugins(3).filterSpec = {'*.mat'};
layerPlugins(3).displayFunc = @displayAnisoFeatures;
layerPlugins(3).dispatchFunc = @dispatchAnisoFeatures;

% [qFSM] Speckles
layerPlugins(4).desc = '[qFSM] Speckles';
layerPlugins(4).filterSpec = {'*.mat'}; 
layerPlugins(4).displayFunc = @displaySpeckles;
layerPlugins(4).dispatchFunc = @dispatchFilesToFrames;

% [panda] Window
layerPlugins(5).desc = '[panda] Windows';
layerPlugins(5).filterSpec = {'*.mat'}; 
layerPlugins(5).displayFunc = @plotWindowsFSM;
layerPlugins(5).dispatchFunc = @dispatchMatrix3ToFrames;

% [Khuloud tracker] Features
layerPlugins(6).desc = '[Khuloud Tracker] Features Info';
layerPlugins(6).filterSpec = {'*.mat'};
layerPlugins(6).displayFunc = @displayIsoFeatures;
layerPlugins(6).dispatchFunc = @dispatchFeaturesInfo;

% [Khuloud tracker] Tracks
layerPlugins(7).desc = '[Khuloud Tracker] Tracks Final';
layerPlugins(7).filterSpec = {'*.mat'};
layerPlugins(7).displayFunc = @displayTracksFinal;
layerPlugins(7).dispatchFunc = @dispatchTracksFinal;

% [Graph] Classified Tracks
layerPlugins(8).desc = '[Graph] Classified Tracks';
layerPlugins(8).filterSpec = {'*.mat'};
layerPlugins(8).displayFunc = @displayClassifiedTracks;
layerPlugins(8).dispatchFunc = @dispatchClassifiedTracks;

end