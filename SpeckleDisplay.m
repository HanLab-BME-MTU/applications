classdef SpeckleDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Markers={'o','^','s','d'};
        LineStyle = '-'
        Color='r';        
    end
    methods
        function obj=SpeckleDisplay(varargin)
            obj@MovieDataDisplay(varargin{:})
        end
        function h=initDraw(obj,data,tag,varargin)
            pos = vertcat(data.Lmax);
            order = [data.speckleType];
            h=-ones(numel(unique(order)),1);
            for i=unique(order)
                h(i)=plot(pos(order==i,2),pos(order==i,1),obj.Markers{i},...
                   'MarkerEdgeColor',obj.Color);%,'MarkerFaceColor',colors(i,:));
            end
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            pos = vertcat(data.Lmax);
            order = [data.speckleType];
            orderMax=max(unique(order));
            for i=1:min(numel(h),orderMax)
                set(h(i),'XData',pos(order==i,2),'YData',pos(order==i,1),...
                    'Marker',obj.Markers{i});%,'MarkerFaceColor',colors(i,:))
            end
            
            for i=numel(h)+1:orderMax
                h(i)=plot(pos(order==i,2),pos(order==i,1),obj.Markers{i},...
                    'MarkerEdgeColor',obj.Color);%,'MarkerFaceColor',colors(i,:));
                set(h,'Tag',tag);
            end           
            
            delete(h(orderMax+1:end));
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Markers';
            params(1).validator=@iscell;
            params(2).name='Color';
            params(2).validator=@ischar;
            params(3).name='LineStyle';
            params(3).validator=@ischar;
        end

        function f=getDataValidator()
            f=@isstruct;
        end
    end    
end