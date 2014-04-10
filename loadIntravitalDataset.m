function [ imageSeries ] = loadIntravitalDataset( dataFilePath )   
    
    data_bfopen = bfopen( dataFilePath );    
    
    if ndims( data_bfopen ) ~= 2 || size(data_bfopen,2) ~= 4
        error( '\n\nERROR: bfopen couldnt read the data properly. The data should be read as m x 4 cell array where m is the number of series in the dataset\n\n' );
    end
    
    imageSeries = struct;
    for sid = 1:size(data_bfopen,1) 

        metadata_bfopen = data_bfopen{sid,4};
        
        % get the metadata information needed
        metadata.numChannels = metadata_bfopen.getPixelsSizeC(0).getValue;
        metadata.numTimePoints = metadata_bfopen.getPixelsSizeT(0).getValue;
        metadata.volSize = [ metadata_bfopen.getPixelsSizeX(0).getValue, metadata_bfopen.getPixelsSizeY(0).getValue metadata_bfopen.getPixelsSizeZ(0).getValue ];
        metadata.voxelSpacing = [ metadata_bfopen.getPixelsPhysicalSizeX(0).getValue metadata_bfopen.getPixelsPhysicalSizeY(0).getValue metadata_bfopen.getPixelsPhysicalSizeZ(0).getValue ];    

        metadata.channelNames = cell(1,metadata.numChannels);
        metadata.channelExcitationWavelength = [];
        for c = 1:metadata.numChannels 
            metadata.channelNames{c} = (metadata_bfopen.getChannelID(sid-1,c-1).toCharArray)';
            if ismethod(metadata_bfopen, 'getChannelExcitationWavelength') && ~isempty( metadata_bfopen.getChannelExcitationWavelength(0,c-1) )
                metadata.channelExcitationWavelength = [metadata.channelExcitationWavelength, metadata_bfopen.getChannelExcitationWavelength(0,c-1).getValue()];
            end
        end
        
        if numel( metadata.channelExcitationWavelength ) < metadata.numChannels
            metadata.channelExcitationWavelength = [];
        end        
        % extract and store data in a conveniently accessible form -- cell
        % array of image stacks
        imageData = cell( metadata.numTimePoints, metadata.numChannels );  
        for t = 1:metadata.numTimePoints         
            for c = 1:metadata.numChannels 

                indx = (t-1) * metadata.volSize(3) * metadata.numChannels + [ c : metadata.numChannels : metadata.volSize(3) * metadata.numChannels ];
                imageData{t,c} = cat(3, data_bfopen{sid,1}{indx, 1});

            end
        end   
        
        imageSeries(sid).metadata = metadata;
        imageSeries(sid).imageData = imageData;
        
    end 
    
end