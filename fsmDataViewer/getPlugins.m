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
layerPlugins(2).dispatchFunc = @dispatchIsoFeatures;

% Aniso Features
layerPlugins(3).desc = 'Aniso Features';
layerPlugins(3).filterSpec = {'*.mat'};
layerPlugins(3).displayFunc = @displayAnisoFeatures;
layerPlugins(3).dispatchFunc = @dispatchAnisoFeatures;

% Segment2D
layerPlugins(4).desc = 'Segment2D';
layerPlugins(4).filterSpec = {'*.mat'};
layerPlugins(4).displayFunc = @displayAnisoFeatures;
layerPlugins(4).dispatchFunc = @dispatchSegment2D;

% [qFSM] Speckles
layerPlugins(5).desc = '[qFSM] Speckles';
layerPlugins(5).filterSpec = {'*.mat'}; 
layerPlugins(5).displayFunc = @displaySpeckles;
layerPlugins(5).dispatchFunc = @dispatchFilesToFrames;

% [qFSM] Flow Analysis
layerPlugins(6).desc = '[qFSM] Flow Analysis';
layerPlugins(6).filterSpec = {'*.mat'}; 
layerPlugins(6).displayFunc = @displayVectorField;
layerPlugins(6).dispatchFunc = @dispatchFlowAnalysis;

% [panda] Window
layerPlugins(7).desc = '[panda] Windows';
layerPlugins(7).filterSpec = {'*.mat'}; 
layerPlugins(7).displayFunc = @plotWindowsFSM;
layerPlugins(7).dispatchFunc = @dispatchMatrix3ToFrames;

% [Khuloud tracker] Features
layerPlugins(8).desc = '[Khuloud Tracker] Features Info';
layerPlugins(8).filterSpec = {'*.mat'};
layerPlugins(8).displayFunc = @displayIsoFeatures;
layerPlugins(8).dispatchFunc = @dispatchFeaturesInfo;

% [Khuloud tracker] Tracks
layerPlugins(9).desc = '[Khuloud Tracker] Tracks Final';
layerPlugins(9).filterSpec = {'*.mat'};
layerPlugins(9).displayFunc = @displayTracksFinal;
layerPlugins(9).dispatchFunc = @dispatchTracksFinal;

% [Graph] Pair Tracks
layerPlugins(10).desc = '[Graph] Pair Tracks';
layerPlugins(10).filterSpec = {'*.mat'};
layerPlugins(10).displayFunc = @displayAnisoFeatures;
layerPlugins(10).dispatchFunc = @dispatchCellToFrames;

% [Graph] Colored Tracks
layerPlugins(11).desc = '[Graph] Colored Tracks';
layerPlugins(11).filterSpec = {'*.mat'};
layerPlugins(11).displayFunc = @displayClassifiedTracks;
layerPlugins(11).dispatchFunc = @dispatchClassifiedTracks;

% [Graph] Connected Component Models
layerPlugins(12).desc = '[Graph] Colored Segments';
layerPlugins(12).filterSpec = {'*.mat'};
layerPlugins(12).displayFunc = @displayClassifiedSegments;
layerPlugins(12).dispatchFunc = @dispatchCellToFrames;

% Theta map
layerPlugins(13).desc = 'Theta Map';
layerPlugins(13).filterSpec = {'*.mat'};
layerPlugins(13).displayFunc = @displayVectorField;
layerPlugins(13).dispatchFunc = @dispatchThetaMap;

end