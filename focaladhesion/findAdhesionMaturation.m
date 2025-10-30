function [indMature,indMatureNAtoFA,indMatureNAtoFC,indMatureFCtoFA,...
        indFail,indFailFC,indFailFA, indStableNA,indStableFC,indStableFA,...
        pNAtoFC,pNAtoFA,pFCtoFA,lifeTimeNAmaturing,lifeTimeNAfailing]=findAdhesionMaturation(tracksNA)

%function [indMature,indMatureNAtoFA,indMatureNAtoFC,indMatureFCtoFA,...
%indFail,indFailFC,indFailFA, indStableNA,indStableFC,indStableFA,...
%pNAtoFC,pNAtoFA,pFCtoFA,lifeTimeNAmaturing,lifeTimeNAfailing]=findAdhesionMaturation(tracksNA)
%classifies all the adhesions in tracksNA with very simple criteria of
%maturation. For example, if the emerging NA become FC, index
%indMatureNAtoFC is 1 at such adhesion. 
% 
% input:    tracksNA        tracksNA struct variable an output of
% AdhesionAnalysisProcess. You can use getTracksNAFromMD(MD) function to
% get this variable too.
% 
% output:   indMature       Nx1 logical vector where all NAs maturing to
% either FC or FA have true values
%           indMatureNAtoFA
%           indMatureNAtoFC
%           indFail
%           indFailFC
%           indFailFA
%           indStableNA
%           indStableFC
%           indStableFA
%           pNAtoFC
%           pNAtoFA
%           pFCtoFA
%           lifeTimeNAmaturing
%           lifeTimeNAfailing
% 
% Sangyoon Han December 2018

if ~isfield(tracksNA,'emerging')
    for k=1:numel(tracksNA)
        % look for tracks that had a state of 'BA' and become 'NA'
        firstNAidx = find(tracksNA(k).state==2,1,'first');
        % see if the state is 'BA' before 'NA' state
        if (~isempty(firstNAidx) && firstNAidx>1 && (tracksNA(k).state(firstNAidx-1)==1)) || (~isempty(firstNAidx) &&firstNAidx==1)
            tracksNA(k).emerging = true;
            tracksNA(k).emergingFrame = firstNAidx;
        else
            tracksNA(k).emerging = false;
        end        
    end
end

indMature = false(numel(tracksNA),1);
indMatureNAtoFA = false(numel(tracksNA),1);
indMatureNAtoFC = false(numel(tracksNA),1);
indMatureFCtoFA = false(numel(tracksNA),1);
indStableNA = false(numel(tracksNA),1);
indStableFC = false(numel(tracksNA),1);
indStableFA = false(numel(tracksNA),1);

indFail = false(numel(tracksNA),1);
indFailFC = false(numel(tracksNA),1);
indFailFA = false(numel(tracksNA),1);
pNAtoFC=0; q=0; qStable=0; qStableExisting=0;
pNAtoFA=0; pFCtoFA=0;
lifeTimeNAmaturing=[];
lifeTimeNAfailing=[];

for k=1:numel(tracksNA)
    if tracksNA(k).emerging 
        % maturing NAs up to FCs and FAs
        if (any(tracksNA(k).state(tracksNA(k).emergingFrame:end)==3) || ...
                any(tracksNA(k).state(tracksNA(k).emergingFrame:end)==4)) && ...
                sum(tracksNA(k).presence)>8

            tracksNA(k).maturing = true;
            indMature(k) = true;
            pNAtoFC=pNAtoFC+1;
            % lifetime until FC
            lifeTimeNAmaturing(pNAtoFC) = sum(tracksNA(k).state(tracksNA(k).emergingFrame:end)==2);
            if any(tracksNA(k).state(tracksNA(k).emergingFrame:end)==4)
                pNAtoFA = pNAtoFA+1;
                pFCtoFA=pFCtoFA+1;
                indMatureNAtoFA(k) = true;
            else
                indMatureNAtoFC(k) = true;
            end
            % it might be beneficial to store amplitude time series. But
            % this can be done later from trackNAmature
        elseif sum(tracksNA(k).presence)<61 && sum(tracksNA(k).presence)>6
        % failing NAs
            tracksNA(k).maturing = false;
            indFail(k) = true;
            q=q+1;
            % lifetime until FC
            lifeTimeNAfailing(q) = sum(tracksNA(k).state(tracksNA(k).emergingFrame:end)==2);
        else
            indStableNA(k) = true;
            % stable NAs
            qStable=qStable+1;
        end
    else %it means that adhesions started already with FC or FA or NA status
        % We check if these adhesion have occasions to convert into FAs
        if tracksNA(k).state(tracksNA(k).startingFrameExtra)==2
            if any(tracksNA(k).state(tracksNA(k).startingFrameExtra:end)==3) % NA becoming FC
                if any(tracksNA(k).state(tracksNA(k).startingFrameExtra:end)==4) % NA becoming FA
                    indMatureNAtoFA(k) = true;
                else % NA becoming only FC
                    indMatureNAtoFC(k) = true;
                end
            else
                % see if this NA ends with NA state
                if tracksNA(k).state(tracksNA(k).endingFrameExtra)==2
                    indStableNA(k) = true;
                else
                    indFail(k) = true;
                end
            end
        elseif tracksNA(k).state(tracksNA(k).startingFrameExtra)==3 %if it started with FC state
            if any(tracksNA(k).state(tracksNA(k).startingFrameExtra:end)==4) % FC becoming FA
                indMatureFCtoFA(k) = true;
            else
                % see if this ends with FC state
                if tracksNA(k).state(tracksNA(k).endingFrameExtra)==3
                    indStableFC(k) = true;
                else % turn-over FC
                    indFailFC(k) = true;
                end
            end
        elseif tracksNA(k).state(tracksNA(k).startingFrameExtra)==4 %if it started with FA state
            % see if this ends with FA state
            [~,areaGrowth]=regression(tracksNA(k).iFrame(tracksNA(k).startingFrame:tracksNA(k).endingFrame), tracksNA(k).area(tracksNA(k).startingFrame:tracksNA(k).endingFrame));
            if tracksNA(k).state(tracksNA(k).endingFrameExtra)==4
                indStableFA(k) = true;
            elseif areaGrowth<0 && ~isnan(tracksNA(k).disassemRate) && tracksNA(k).disassemRate>0 && ...
                    ismember(tracksNA(k).state(end), [2,3,5])% turn-over FA if it is decreasing size and intensity and ending with FC, NA then ANA
                indFailFA(k) = true;
            end
        end
%             if (any(tracksNA(k).state==2) || any(tracksNA(k).state==3)) 
%                 pFCtoFA=pFCtoFA+1;
%                 indMatureFCtoFA(k) = true;
%             else
%                 qStableExisting=qStableExisting+1;
%                 indStableNA(k) = true;
%             end
    end
end
