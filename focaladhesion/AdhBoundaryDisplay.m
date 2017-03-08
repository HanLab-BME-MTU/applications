classdef AdhBoundaryDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        LineStyle = '-';
        LineWidth = 0.5;
        Color = 'k';  
    end
    methods
        function obj=AdhBoundaryDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj, data, tag, varargin)
            
            
            if isempty(data.number), h=[]; return; end            
            h = gobjects(size(data.number,1),1);
            
            for i=1:numel(data.number), 
                adhBoundary = data.adhBoundary{i};
                h(i) = plot(adhBoundary(:,2), adhBoundary(:,1), varargin{:});
            end
            set(h,'Tag',tag);
        end


        function updateDraw(obj, h, data)
            tag=get(h(1),'Tag');
            obj.initDraw(data,tag,h);
            return;

            
        % function updateDraw(obj,h,data)
        %     tag=get(h(1),'Tag');
        %     nTracks = size(data.x,1);
            
        %     % Delete tracks
        %     delete(h(nTracks+1:end));
        %     h(nTracks+1:end)=[];
        %     if nTracks==0, return; end
            
        %     %
        %     existingTracks=false(nTracks,1);
        %     existingTracks(1:min(numel(h),nTracks))=true;
        %     uniqueTypes=unique(data.trackType);
        %     for i=1:numel(uniqueTypes), 
        %         index=data.trackType==uniqueTypes(i);
        %         index1=index & existingTracks;
        %         for j=find(index1)'
                    
        %             set(h(j),'XData',data.x(j,max(1,end-obj.dragtailLength):end),...
        %                 'YData',data.y(j,max(1,end-obj.dragtailLength):end)',...
        %                 'Color',obj.Color(uniqueTypes(i)),'LineStyle',obj.LineStyle(uniqueTypes(i)));
        %         end
                
        %         index2= index & ~existingTracks;
        %         h(index2) = plot(data.x(index2,max(1,end-obj.dragtailLength):end)',...
        %             data.y(index2,max(1,end-obj.dragtailLength):end)',...
        %             'Color',obj.Color(uniqueTypes(i)),'LineStyle',obj.LineStyle(uniqueTypes(i)),'LineWidth',obj.LineWidth);
        %     end
            
        %     % Set tag
        %     set(h,'Tag',tag);          
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