classdef FASpotDisplay < MovieDataDisplay
    % Concrete display class for displaying FA detection points
    % Andrew R. Jamieson - Mar. 2017

    properties
        % FA_types = {'BA','NA','FC','FA'};
        ColorDict = containers.Map({'BA','NA','FC','FA'}, {'g','r', 'o', 'b'})
        Marker = 'o';
        MarkerSize = 6; 
        LineStyle = '-';
        LineWidth = 1;
        XLabel='';
        YLabel='';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
        ButtonDownFcn=[];
    end
    methods
                        
        function obj=LineDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)

            if isempty(data.x), h=[]; return; end            
            h = gobjects(size(data.x,1),1);
                        
            for FAtype = obj.ColorDict.keys  
                rows = data.state == 'FC';
                vars = {'xCoord','yCoord'};
                h(index) = initDraw@LineDisplay(obj, data{rows,vars}, tag, 'Color', obj.ColorDict(FAtype), varargin{:});
            end

            set(h,'Tag',tag);
        end
        
        function updateDraw(obj,h,data)
            % Update handle xData and yData
            set(h,'XData',data(:,1),'YData',data(:,2));
            obj.setLineProperties(h);
            obj.setAxesProperties();
        end
        
        function setLineProperties(obj, h)
            set(h, 'MarkerSize', obj.MarkerSize,...
                'Color', obj.Color, 'Marker',obj.Marker,...
                'Linestyle', obj.LineStyle, 'LineWidth', obj.LineWidth,...
                'ButtonDownFcn', obj.ButtonDownFcn);
        end
        
        function setAxesProperties(obj)
            % Set labels and fonts
            if ~isempty(obj.XLabel),xlabel(obj.XLabel,obj.lfont{:}); end
            if ~isempty(obj.YLabel),ylabel(obj.YLabel,obj.lfont{:}); end
            set(gca,'LineWidth', 1.5, obj.sfont{:})
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='ColorDict';
            params(1).validator=@(x)isa(x,'containers.Map');
            params(2).name='Marker';
            params(2).validator=@ischar;
            params(3).name='LineStyle';
            params(3).validator=@ischar;
            params(4).name='LineWidth';
            params(4).validator=@isscalar;
            params(5).name='XLabel';
            params(5).validator=@ischar;
            params(6).name='YLabel';
            params(6).validator=@ischar;
            params(7).name='sfont';
            params(7).validator=@iscell;
            params(8).name='lfont';
            params(8).validator=@iscell;
            params(9).name='MarkerSize';
            params(9).validator=@isscalar;
            params(10).name='ButtonDownFcn';
            params(10).validator=@(x) isempty(x) || isa(x, 'function_handle');
        end
        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end