classdef FASpotDisplay < MovieDataDisplay
    % Concrete display class for displaying FA detection points
    % Andrew R. Jamieson - Mar. 2017

    properties
        % FA_types = {'BA','NA','FC','FA'};
        ColorDict = containers.Map({'BA','NA','FC','FA'}, {'g','r', 'y', 'b'})
        Marker = 'o';
        MarkerSize = 5; 
        LineStyle = 'none';
        LineWidth = .5;
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
                rows = data.state == FAtype{1};
                vars = {'xCoord','yCoord'};
                d = data{rows,vars};
                if ~isempty(d)
                    h(index) = plot(d(:,1) ,d(:,2), 'Color',...
                        obj.ColorDict(FAtype{1}), varargin{:},...
                        'Linestyle', obj.LineStyle, 'LineWidth', obj.LineWidth,...
                        'MarkerSize', obj.MarkerSize, 'Marker',obj.Marker);
                    set(h(index),'Tag',tag);                    
                    index = 1+index;
                end
            end
            
            % Necessary?
            obj.setAxesProperties();

        end
        
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            
            delete(h(:));
            
            index = 1;
            for FAtype = obj.ColorDict.keys  
                
                rows = data.state == FAtype{1};
                vars = {'xCoord','yCoord'};
                d = data{rows,vars};
                if ~isempty(d)
                    h(index) = plot(d(:,1) ,d(:,2), 'Color',...
                        obj.ColorDict(FAtype{1}),... 
                        'Linestyle', obj.LineStyle, 'LineWidth', obj.LineWidth,...
                        'MarkerSize', obj.MarkerSize, 'Marker',obj.Marker);
                    set(h(index),'Tag',tag);                    
                    index = 1+index;
                end
            end            
%             obj.setLineProperties(h);
            obj.setAxesProperties();
        end
%         
%         function setLineProperties(obj, h)
%             set(h, 'MarkerSize', obj.MarkerSize,...
%                 'Color', obj.Color, 'Marker',obj.Marker,...
%                 'Linestyle', obj.LineStyle, 'LineWidth', obj.LineWidth,...
%                 'ButtonDownFcn', obj.ButtonDownFcn);
%         end
        
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
            f=@istable;
        end
    end    
end