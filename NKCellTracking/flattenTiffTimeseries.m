function [dataOut,metadata] = flattenTiffTimeseries(filename)
% FLATTENTIFFIMAGE takes a .tiff timeseries of 3D image stacks and
% averages the z-dimension to produce a timeseries of 2D images.
% 
%   Usage: []=flattenTiffTimeseries(filename)
%   Input:
%           filename:       filename of the multipage tiff to be flattened
%           newFilename:    filename of the resulting tiff image to be
%                           written into the same directory as the original
%   Output:
%           flatData:       3D matrix representing the timeseries of 2D images
%
%'/storage/disk1/sehaarma/NKcell Data Flattening & Tracking/high vs low cytotoxicity/low cytotoxicity-HT29.tif'

%read OME-TIFF file using bio-formats
data = bfopen(filename);

metadata = data{1,2};
metadataKeys = metadata.keySet().iterator();
for i=1:metadata.size()
  key = metadataKeys.nextElement();
  value = metadata.get(key);
  if isnumeric(value)
      value = num2str(value);
  end
  fprintf('%s = %s\n', key, value)
end
desiredmeta = metadata;

disp(data{1,1}{1,2})
%numFrames = input('How many time frames does the series have?');
numFrames = 273;

dataOut = zeros(512,512,1,4,numFrames,'uint16');
%dataOut = zeros(512,512,1,4,273,'uint8');

g = 1;
for i = 1:numFrames
    %load each z-frame for each color from current time frame
    z1c1 = data{1,1}{g,1};
    z1c2 = data{1,1}{g+1,1};
    z1c3 = data{1,1}{g+2,1};
    z1c4 = data{1,1}{g+3,1};
    z2c1 = data{1,1}{g+4,1};
    z2c2 = data{1,1}{g+5,1};
    z2c3 = data{1,1}{g+6,1};
    z2c4 = data{1,1}{g+7,1};
    z3c1 = data{1,1}{g+8,1};
    z3c2 = data{1,1}{g+9,1};
    z3c3 = data{1,1}{g+10,1};
    z3c4 = data{1,1}{g+11,1};
    %concatenate z-frames into 3D matrix
    c1 = cat(3,z1c1,z2c1,z3c1);
    c2 = cat(3,z1c2,z2c2,z3c2);
    c3 = cat(3,z1c3,z2c3,z3c3);
    c4 = cat(3,z1c4,z2c4,z3c4);
    %average across the 3rd matrix dimension to flatten
    c1 = mean(c1,3);
    c2 = mean(c2,3);
    c3 = mean(c3,3);
    c4 = mean(c4,3);
    %insert each color into data output matrix
    dataOut(:,:,1,1,i) = c1;
    dataOut(:,:,1,2,i) = c2;
    dataOut(:,:,1,3,i) = c3;
    dataOut(:,:,1,4,i) = c4;
    g = g + 12;
end

metadata = createMinimalOMEXMLMetadata(dataOut);
%camera bit depth
significantBits = ome.xml.model.primitives.PositiveInteger(java.lang.Integer(16));
metadata.setPixelsSignificantBits(significantBits,0);
%time of experiment
sizeT = ome.xml.model.primitives.PositiveInteger(java.lang.Integer(numFrames));
metadata.setPixelsSizeT(sizeT,0);
%timeIncrement = ome.units.quantity.Time(java.lang.Double(317.116), ome.units.UNITS.S);
metadata.setPixelsTimeIncrement(ome.units.quantity.Time(java.lang.Double(317.1158447265625), ome.units.UNITS.S), 0);
%physical pixel size
voxelSize = ome.units.quantity.Length(java.lang.Double(2.4748064701340358), ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeX(voxelSize, 0);
metadata.setPixelsPhysicalSizeY(voxelSize, 0);
voxelSizeZ = ome.units.quantity.Length(java.lang.Double(18.025), ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeZ(voxelSizeZ, 0);
%interleaved
metadata.setPixelsInterleaved(java.lang.Boolean("false"), 0);
%image name
metadata.setImageName(java.lang.String("high cytotoxicity-LoVo_2Dtimeseries.tif"), 0);

desiredmeta = data{1,4};
desiredmeta.setPixelsSizeZ(ome.xml.model.primitives.PositiveInteger(java.lang.Integer(1)), 0);
desiredmeta.setImageName(java.lang.String("high cytotoxicity-LoVo_2Dtimeseries.tif"), 0);

bfsave(dataOut,strcat(filename(1:end-4),'_2Dtimeseries.tif'),'metadata',metadata);
end
