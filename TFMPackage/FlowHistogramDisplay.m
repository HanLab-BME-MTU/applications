classdef FlowHistogramDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        Linewidth = 2;
        nBins=20;
    end
    methods
        function obj=FlowHistogramDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            % Read data number
            nData=numel(data);
            dataIdx=[1 round(nData/2) nData];
            colors =hsv(numel(dataIdx));
            
            % define small and large fonts     
            tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
            sfont = {'FontName', 'Helvetica', 'FontSize', 18};
            lfont = {'FontName', 'Helvetica', 'FontSize', 22};
            
            % Generate plot
            hold on;
            h=-1*ones(numel(dataIdx),1);
            for i=1:numel(dataIdx),
                [n,x]=hist(data{dataIdx(i)},obj.nBins);
                h(i)=plot(x,n,'Color',colors(i,:),'Linewidth',obj.Linewidth);
            end
            legend(arrayfun(@(x) ['Frame ' num2str(x)],dataIdx,'UniformOutput',false),...
                'Location', 'NorthEast', tfont{:});
            xlabel('Flow (pixels/frame)',lfont{:});
            ylabel('Number',lfont{:});
            set(gca, 'LineWidth', 1.5, sfont{:});
            
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            tag = get(h(1),'Tag');
            cla(get(get(h(1),'Parent')))
            obj.initDraw(data,tag);
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Marker',obj.Marker,@ischar); 
            ip.addParamValue('Linewidth',obj.Linewidth,@isscalar);  
            ip.addParamValue('nBins',obj.nBins,@isscalar);  
        end 
        function setProperties(obj,ip)
            obj.Marker=ip.Results.Marker;
            obj.Linewidth=ip.Results.Linewidth;
            obj.nBins=ip.Results.nBins;
        end
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@iscell;
        end
    end    
end