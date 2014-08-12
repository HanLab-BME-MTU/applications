classdef BioImageReader

    methods (Access = public)
        
        function this = BioImageReader(dataFilePath)
            
            this.dataFilePath = dataFilePath;
            this.r = bfGetReader(dataFilePath);    
            this.metadataList = this.extractMetadata();
            
        end
        
        function numSeries = getSeriesCount(this)
            numSeries = this.r.getSeriesCount();
        end
            
        function metadata = getMetadata(this, seriesId)
            
           if ~exist('seriesId', 'var')  
              metadata = this.metadataList;
           else
              metadata = this.metadataList(seriesId); 
           end
            
        end

        function imageData = getImageData(this, seriesId, varargin)

            p = inputParser;
            p.CaseSensitive = false;
            
            p.addRequired('seriesId', @(x) (isnumeric(x) && isvector(x) && x >= 1 && x <= this.getSeriesCount()) );
            metadata = this.metadataList(seriesId);
            
            p.addParamValue('timepointIds', [], @(x) (isnumeric(x) && isvector(x) && all(ismember(x,1:metadata.numTimePoints))) );
            p.addParamValue('channelIds', [], @(x) (isnumeric(x) && isvector(x) && all(ismember(x,1:metadata.numChannels))) );
            p.addParamValue('zPlanes', [], @(x) (numel(metadata.imageSize) > 1 && isnumeric(x) && isvector(x) && all(ismember(x,1:metadata.imageSize(3)))) );
            p.parse(seriesId, varargin{:});
            
            PARAMETERS = p.Results;
            
            if isempty(PARAMETERS.timepointIds)
                timepointIds = 1:metadata.numTimePoints;
            else
                timepointIds = PARAMETERS.timepointIds;
            end
               
            if isempty(PARAMETERS.channelIds)
                channelIds = 1:metadata.numChannels;
            else
                channelIds = PARAMETERS.channelIds;
            end
            
            if isempty(PARAMETERS.zPlanes)
                if numel(metadata.imageSize) > 1
                    zPlanes = 1:metadata.imageSize(3);
                else
                    zPlanes = 1;
                end
            else
                zPlanes = PARAMETERS.zPlanes;
            end
            
            imageData = cell(numel(timepointIds), numel(channelIds)); 

            for i = 1:numel(timepointIds)           
               tid = timepointIds(i); 
               for j = 1:numel(channelIds)
                    chid = channelIds(j);
                    imageData{i,j} = [];
                    for k = 1:numel(zPlanes)
                        pid = zPlanes(k);
                        planeId = this.r.getIndex(pid-1, chid-1, tid-1) + 1;
                        imPlane = bfGetPlane(this.r, planeId);
                        imageData{i,j} = cat(ndims(imPlane)+1, imageData{i,j}, imPlane);
                   end
               end
            end
            
        end
        
    end
    
    methods(Access = private)
    
        function [metadataList] = extractMetadata(this)
            
            metadataList = [];
            
            for sid = 1:this.getSeriesCount()
                
                % basic metadata
                this.r.setSeries(sid-1);

                metadata.version = '2.0.0';
                metadata.dataFilePath = this.dataFilePath;
                
                metadata.format = char( this.r.getFormat() );                

                metadata.numSeries = this.getSeriesCount();
                metadata.seriesId = sid;
                
                metadata.numTimePoints = this.r.getSizeT();
                
                metadata.numChannels = this.r.getSizeC();
                metadata.channelNames = cell(1, metadata.numChannels);
                for chid = 1:metadata.numChannels
                    cName = char( this.r.getMetadataStore().getChannelName(sid-1, chid-1) );
                    if isempty(cName)
                        cName = sprintf('ch%d', chid);
                    end
                    metadata.channelNames{chid} = cName;
                end

                if this.r.getSizeZ() > 1
                    metadata.imageSize = [this.r.getSizeX(), this.r.getSizeY(), this.r.getSizeZ()]; 
                    try
                        metadata.pixelSize = [this.r.getMetadataStore().getPixelsPhysicalSizeX(0).getValue(), ...
                                              this.r.getMetadataStore().getPixelsPhysicalSizeY(0).getValue(), ...
                                              this.r.getMetadataStore().getPixelsPhysicalSizeZ(0).getValue()];
                    catch ME
                        metadata.pixelSize = ones(1, 3);
                    end
                                              
                else
                    metadata.imageSize = [this.r.getSizeX(), this.r.getSizeY()];    
                    try
                        metadata.pixelSize = [this.r.getMetadataStore().getPixelsPhysicalSizeX(0).getValue(), ...
                                              this.r.getMetadataStore().getPixelsPhysicalSizeY(0).getValue()];
                    catch ME
                        metadata.pixelSize = ones(1, 2);
                    end
                end


                metadata.pixelType = char(loci.formats.FormatTools.getPixelTypeString( this.r.getPixelType() ));
                metadata.bitsPerPixel = this.r.getBitsPerPixel();   
                
                % metadata.omeMetadata = this.r; % this is not serializable
                
                metadataList = [metadataList; metadata];
                
            end
            
        end
        
    end
    
    properties (SetAccess = private)
       
        r
        metadataList
        dataFilePath
        
    end

end