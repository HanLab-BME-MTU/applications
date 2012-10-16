classdef PairsDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Linestyle='-';
        GapLinestyle='--';
        Color='r';
        showLabel=false;
    end
    methods
        function obj=PairsDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            nPairs = numel(data);
            h=-ones(nPairs,4);
            % List valid pairs
            validPairs = find(~arrayfun(@(x) all(isnan(x.coords1)),data));
            delta = -5;
            
            for i = validPairs'
                % Concatenate coordinates from positions 1 and 2
                coords = vertcat(data(i).coords1(1:3),data(i).coords2(1:3));

                % Plot pairs
                h(i,1)=plot3(coords(:,1),coords(:,2),coords(:,3),'-',...
                    'Color',obj.Color,'Linestyle',obj.Linestyle);
                
                % Draw EB intensities as circles
                if isfield(data(i),'kEBcoords1') && ~isnan(data(i).kEBcoords1(1))
                    h(i,2) =plot3(data(i).kEBcoords1(1,1),data(i).kEBcoords1(1,2),0,'o',...
                        'Color','r','MarkerSize',ceil(data(i).kEBamp1(1)*0.1));
                end
                if isfield(data(i),'kEBcoords2') && ~isnan(data(i).kEBcoords2(1))
                    h(i,3) =plot3(data(i).kEBcoords2(1,1),data(i).kEBcoords2(1,2),0,'o',...
                        'Color','r','MarkerSize',ceil(data(i).kEBamp2(1)*0.1));
                end
                
                % Display pair numbers if option is selected
                if obj.showLabel
                    h(i,4) = text(mean(coords(:,1)) + delta, mean(coords(:,2)) + delta, ...
                        mean(coords(:,3)) + delta, num2str(i), 'Color', obj.Color);
                end
            end
            
            set(h(ishandle(h)),'Tag',tag);
        end
        
        function updateDraw(obj,allh,data)
            tag=get(allh(1),'Tag');
            delete(allh);
            obj.initDraw(data,tag);
            return;
        end
    end
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@(x)ischar(x) ||isvector(x);
            params(2).name='Linestyle';
            params(2).validator=@ischar;
            params(3).name='GapLinestyle';
            params(3).validator=@ischar;
            params(4).name='showLabel';
            params(4).validator=@isscalar;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end
end