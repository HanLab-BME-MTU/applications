classdef AdhBoundaryDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        LineStyle = '-';
        LineWidth = 0.5;
        Color = 'm';  
    end
    methods
        function obj=AdhBoundaryDisplay(varargin)            
            obj = obj@MovieDataDisplay(varargin{:});
        end
        function h=initDraw(obj, data, tag, varargin)
            
            if isempty(data), h=[]; return; end            
            h = gobjects(size(data,1),1);
            for i=1:numel(data) 
                adhBoundary = data(i).adhBoundary;
                try
                    h(i) = line(adhBoundary(:,2), adhBoundary(:,1),...
                                'Color',obj.Color,'LineStyle',obj.LineStyle,...
                                varargin{:});
                catch
                    disp('error');
                end
            end
            set(h,'Tag',tag);
        end

        function updateDraw(obj, h, data, varargin)
            tag = get(h(1),'Tag');
            nTracks = size(data,1);
            
            % Delete tracks
            delete(h(nTracks+1:end));
            h(nTracks+1:end) = [];
            if nTracks==0, return; end
            
            existingTracks = false(nTracks,1);
            existingTracks(1:min(numel(h),nTracks)) = true;
            
            for j = 1:nTracks
                
                adhBoundary = data(j).adhBoundary;
                
                if existingTracks  
                    set(h(j),'XData',adhBoundary(:,2), 'YData', adhBoundary(:,1),...
                             'Color',obj.Color,'LineStyle',obj.LineStyle,varargin{:});
                else
                    h(j) = line(adhBoundary(:,2), adhBoundary(:,1),...
                                'Color',obj.Color,'LineStyle',obj.LineStyle,...
                                varargin{:});               
                end
            end
            set(h,'Tag',tag);          
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@(x)(ischar(x) || (numel(x)==3 && isnumeric(x)));
            params(2).name='LineStyle';
            params(2).validator=@ischar;
            params(3).name='LineWidth';
            params(3).validator=@isscalar;
        end

        function f=getDataValidator() 
            f=@isstruct;
        end
    end    
end