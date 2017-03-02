classdef FASpotDisplay < MovieDataDisplay
    % Concrete display class for displaying FA detection points
    % Andrew R. Jamieson - Mar. 2017

    properties
        % FA_types = {'BA','NA','FC','FA'};
        ColorDict = containers.Map({'BA','NA','FC','FA'}, {'g','r', 'y', 'b'})
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
                        
        function obj=FASpotDisplay(varargin)
            
%             ip =inputParser;
%             ip.addRequired('obj');
%             ip.addParameter('ColorDict',ColorDict,@islogical);
%             ip.addParameter('output', outputList{1}, @(x) all(ismember(x,outputList)));
%             ip.parse(obj,iChan,varargin{:})
%             output = ip.Results.output;
%             iFrame = ip.Results.iFrame;
            
            
            obj = obj@MovieDataDisplay(varargin{:});
%             nVarargin = numel(varargin);
%             if nVarargin > 1 && mod(nVarargin,2)==0
%                 for i=1 : 2 : nVarargin-1
%                     obj.(varargin{i}) = varargin{i+1};
%                 end
%             end
        end
        function h=initDraw(obj,data,tag,varargin)

            if isempty(data.xCoord), h=[]; return; end            
            h = gobjects(numel(obj.ColorDict.keys),1);
                        
            index = 1;
            for FAtype = obj.ColorDict.keys  
                rows = data.state == FAtype{1} & data.pres == true;
                vars = {'xCoord','yCoord'};
                d = data{rows,vars};
                if ~isempty(d)
%                     FAtype
%                     obj.ColorDict(FAtype{1})
                    h(index) = plot(d(:,1) ,d(:,2), 'Color', obj.ColorDict(FAtype{1}), varargin{:});
                    set(h(index),'Tag',tag);                    
                    index = 1+index;
                end
            end

        end
        
        function updateDraw(obj,h,data)
            % Update handle xData and yData
            set(h,'XData',data{:,1},'YData', data{:,2});
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
    
        % function h=draw(obj,data,tag,varargin)
        %     % Template method to draw a movie data component
            
        %     % Check input
            
        %     ip =inputParser;
        %     ip.addRequired('obj',@(x) isa(x,'MovieDataDisplay'));
        %     ip.addRequired('data',obj.getDataValidator());
        %     ip.addRequired('tag',@ischar);
        %     ip.addParamValue('hAxes',gca,@ishandle);
        %     params = obj.getParamValidators;
        %     for i=1:numel(params)
        %         ip.addParamValue(params(i).name,obj.(params(i).name),params(i).validator);
        %     end
        %     ip.KeepUnmatched = true; % Allow unmatched arguments
        %     ip.parse(obj,data,tag,varargin{:});
        %     for i=1:numel(params)
        %         obj.(params(i).name)=ip.Results.(params(i).name);
        %     end
            
        %     % Retrieve the axes handle and call the create figure method 
        %     hAxes = ip.Results.hAxes;
        %     set(hAxes,'NextPlot','add');
            
        %     % Get the component handle and call the adapted draw function
        %     h = findobj(hAxes,'-regexp','Tag',['^' tag '$']);
        %     if ~isempty(h) && any(ishandle(h))
        %         obj.updateDraw(h,data);
        %     else
        %         h=obj.initDraw(data,tag,'Parent',hAxes);
        %     end
        % end


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
            f=@istable;
        end
    end    
end