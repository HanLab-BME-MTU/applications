function [imageData, mask, metaData] = handleData(handles)
%handleData extracts relevant parameters from handles for afterward
%analysis functions, like findNulei and singleNucleusSpotDetection

%   handles is a parameter from InvivoCytometer_2.0_source_code/code_package/CellSegmentationQualityAnnotator.m
%   03/2016 Ning

% Multi-channel 3D stack
imageData = handles.data.imageData; 
% 3D stack for Dapi channel
mask = handles.data.imLabelCellSeg;
% Multi-channel 3D stack information
metaData = handles.data.metadata;


end

