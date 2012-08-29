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
            h=-ones(nPairs,3);
            
            for i=1:nPairs
                coords = vertcat(data(i).coords1(1:3),data(i).coords2(1:3));
                h(i,1)=plot3(coords(:,1),coords(:,2),coords(:,3),'-',...
                    'Color',obj.Color,'Linestyle',obj.Linestyle);
                if isfield(data(i),'kMTcoords1') && ~isnan(data(i).kMTcoords1(1))
                   h(i,2) =plot3(data(i).kMTcoords1(1,1),data(i).kMTcoords1(1,2),0,'o',...
                    'Color','r','MarkerSize',ceil(data(i).kMTamp1(1)*10^3));
                end
                if isfield(data(i),'kMTcoords2') && ~isnan(data(i).kMTcoords2(1))
                   h(i,3) =plot3(data(i).kMTcoords2(1,1),data(i).kMTcoords2(1,2),0,'o',...
                    'Color','r','MarkerSize',ceil(data(i).kMTamp2(1)*10^3));
                end
            end
            set(h(ishandle(h)),'Tag',tag);            
        end

        function updateDraw(obj,allh,data)
            tag=get(allh(1),'Tag');
            delete(allh);
            obj.initDraw(data,tag);
            return;
            nTracks = numel(data);

            h=findobj(allh,'Type','line');
            delete(h(2*nTracks+1:end));
            h(2*nTracks+1:end)=[];
            hlinks=findobj(h,'LineStyle',obj.Linestyle);
            hgaps=findobj(h,'LineStyle',obj.GapLinestyle);
            
            % Update existing windows
            for i=1:min(numel(hlinks),nTracks) 
                xData= data.x{i}(max(1,end-obj.dragtailLength):end);
                yData= data.y{i}(max(1,end-obj.dragtailLength):end);
                set(hgaps(i),'Xdata',xData(~isnan(xData)),'YData',yData(~isnan(yData)));
                set(hlinks(i),'Xdata',xData,'YData',yData);
            end
            for i=min(numel(hlinks),nTracks)+1:nTracks
                xData= data.x{i}(max(1,end-obj.dragtailLength):end);
                yData= data.y{i}(max(1,end-obj.dragtailLength):end);
                hgaps(i)=plot(xData(~isnan(xData)),yData(~isnan(yData)),...
                    'Linestyle',obj.GapLinestyle','Color',obj.Color);
                hlinks(i)=plot(xData,yData,...
                    'Linestyle',obj.Linestyle,'Color',obj.Color);
            end
            set([hlinks hgaps],'Tag',tag);

            
            if isfield(data,'label') && obj.showLabel
                ht=findobj(allh,'Type','text');
                delete(ht(nTracks+1:end));
                ht(nTracks+1:end)=[];
                for i=1:min(numel(ht),nTracks)
                    set(ht(i),'Position',[data.x{i}(end),data.y{i}(end)],...
                        'String',data.label(i));
                end
                for i=min(numel(ht),nTracks)+1:nTracks
                    ht(i) = text(data.x{i}(end),data.y{i}(end),num2str(data.label(i)));
                end
                set(ht,'Tag',tag); 
            end
           
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