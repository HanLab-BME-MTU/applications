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

% Points
layerPlugins(2).desc = 'Points';
layerPlugins(2).filterSpec = {'*.mat'};
layerPlugins(2).displayFunc = @displayPoints;
layerPlugins(2).dispatchFunc = @dispatchCellToFrames;

% [qFSM] Speckles
layerPlugins(3).desc = '[qFSM] Speckles';
layerPlugins(3).filterSpec = {'*.mat'}; 
layerPlugins(3).displayFunc = @displaySpeckles;
layerPlugins(3).dispatchFunc = @dispatchFilesToFrames;

% [panda] Window
layerPlugins(4).desc = '[panda] Windows';
layerPlugins(4).filterSpec = {'*.mat'}; 
layerPlugins(4).displayFunc = @plotWindowsFSM;
layerPlugins(4).dispatchFunc = @dispatchMatrix3ToFrames;

% [Khuloud tracker] Tracks
layerPlugins(5).desc = '[Khuloud Tracker] Feature Info';
layerPlugins(5).filterSpec = {'*.mat'};
layerPlugins(5).displayFunc = @displayPoints;
layerPlugins(5).dispatchFunc = @dispatchKhuloudFeatures;

% [Khuloud tracker] Tracks
layerPlugins(6).desc = '[Khuloud Tracker] Tracks';
layerPlugins(6).filterSpec = {'*.mat'};
layerPlugins(6).displayFunc = @displayKhuloudTracks;
layerPlugins(6).dispatchFunc = @dispatchKhuloudTracks;

end