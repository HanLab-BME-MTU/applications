classdef SpeckleDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'o';
        LineStyle = '-'
        Color='r';        
        maxOrder=4;
    end
    methods
        function obj=SpeckleDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            pos = vertcat(data.Lmax);
            order = [data.speckleType];
            colors = hsv(obj.maxOrder);
            h=-ones(numel(unique(order)),1);
            for i=unique(order)
                h(i)=plot(pos(order==i,2),pos(order==i,1),'o',...
                   'MarkerEdgeColor',colors(i,:));%,'MarkerFaceColor',colors(i,:));
            end
            set(h,'Tag',tag,'Marker',obj.Marker);
        end
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            pos = vertcat(data.Lmax);
            order = [data.speckleType];
            orderMax=max(unique(order));
            colors = hsv(obj.maxOrder);
            for i=1:min(numel(h),orderMax)
                set(h(i),'XData',pos(order==i,2),'YData',pos(order==i,1),...
                    'MarkerEdgeColor',colors(i,:));%,'MarkerFaceColor',colors(i,:))
            end
            
            for i=numel(h)+1:orderMax
                h(i)=plot(pos(order==i,2),pos(order==i,1),'o',...
                    'MarkerEdgeColor',colors(i,:));%,'MarkerFaceColor',colors(i,:));
                set(h,'Tag',tag,'Marker',obj.Marker);
            end           
            
            delete(h(orderMax+1:end));
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@ischar);
            ip.addParamValue('Marker',obj.Marker,@ischar);
            ip.addParamValue('LineStyle',obj.LineStyle,@ischar);  
        end 
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.Marker=ip.Results.Marker;
            obj.LineStyle=ip.Results.LineStyle;
        end
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isstruct;
        end
    end    
end