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
            % For classification, there is classLabel at the second cell.
            % Check with it
            numGroups=9;
            colors = distinguishable_colors(numGroups,'k');
            tempColor = colors(6,:);
            colors(6,:) = colors(9,:);
            colors(9,:) = tempColor;
            idGroupLabel=[];
            if iscell(data{1})
                idGroupLabel=data{2};
                data = data{1};
            end
            h = gobjects(size(data,1),1);
            
            if isempty(idGroupLabel)
                for i=1:numel(data) 
                    adhBoundary = data{i}; %(i).adhBoundary;
                    try
                        h(i) = line(adhBoundary(:,2), adhBoundary(:,1),...
                                    'Color',obj.Color,'LineStyle',obj.LineStyle,...
                                    varargin{:});
                    catch
                        disp('error');
                    end
                end
            else
                existingClasses=unique(idGroupLabel','sorted');
                j=0;
                for k = existingClasses
                    curIndices = find(idGroupLabel'==k);
                    for i=curIndices
                        j=j+1;
                        adhBoundary = data{i}; %(i).adhBoundary;
                        try
                            h(j) = line(adhBoundary(:,2), adhBoundary(:,1),...
                                        'Color',colors(k,:),'LineStyle',obj.LineStyle,...
                                        varargin{:});
                        catch
                            disp('error');
                        end
                    end
                end
            end
            set(h,'Tag',tag);
        end

        function updateDraw(obj, h, data, varargin)
            tag = get(h(1),'Tag');
            idGroupLabel=[];
            if ~isempty(data)
                if iscell(data{1})
                    idGroupLabel=data{2};
                    data = data{1};
                end
            end
            nTracks = numel(data); %size(data,1);
            
            % Delete tracks
            delete(h(nTracks+1:end));
            h(nTracks+1:end) = [];
            if nTracks==0, return; end
            
            existingTracks = false(nTracks,1);
            existingTracks(1:min(numel(h),nTracks)) = true;

            numGroups=9;
            colors = distinguishable_colors(numGroups,'k');
            tempColor = colors(6,:);
            colors(6,:) = colors(9,:);
            colors(9,:) = tempColor;
%             h = gobjects(size(data,1),1);
            
            if isempty(idGroupLabel)
                for j = 1:nTracks

                    adhBoundary = data{j}; %(j).adhBoundary;

                    if existingTracks(j)  
                        set(h(j),'XData',adhBoundary(:,2), 'YData', adhBoundary(:,1),...
                                 'Color',obj.Color,'LineStyle',obj.LineStyle,varargin{:});
                    else
                        h(j) = line(adhBoundary(:,2), adhBoundary(:,1),...
                                    'Color',obj.Color,'LineStyle',obj.LineStyle,...
                                    varargin{:});               
                    end
                end
            else
                % In case there is zero label, make it six (noise)
                idGroupLabel(idGroupLabel==0)=6;
                % if there is anything related to more than group 9, those
                % are assigned to be group 6 (noise)
                idGroupLabel(idGroupLabel>9)=6;
                existingClasses=unique(idGroupLabel','sorted');
                
                j=0;
                for k = existingClasses
                    curIndices = find(idGroupLabel'==k);
                    for i=curIndices
                        j=j+1;
                        adhBoundary = data{i}; %(i).adhBoundary;
                        if existingTracks(j)  
                            set(h(j),'XData',adhBoundary(:,2), 'YData', adhBoundary(:,1),...
                                     'Color',colors(k,:),'LineStyle',obj.LineStyle,varargin{:});
                        else
                            h(j) = line(adhBoundary(:,2), adhBoundary(:,1),...
                                        'Color',colors(k,:),'LineStyle',obj.LineStyle,...
                                        varargin{:});
                        end
                    end
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
            f=@iscell;
        end
    end    
end