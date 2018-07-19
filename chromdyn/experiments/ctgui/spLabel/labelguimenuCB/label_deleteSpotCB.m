function label_deleteSpotCB(spid)
%called from popup-menu. Allows to delete either spot (forced recalc) or
%color (optional recalc)

%switch delete spot or color
what2del = questdlg('Do you want to delete this specific spot or do you want to delete this color in all frames?',...
    'Do you know what you''re doing?','delete spot','delete color','CANCEL','CANCEL');

switch what2del
    
    case {'CANCEL',''}
        hpop = findall(0,'Style','listbox','Tag','selcetlabel');
        if ~isempty(hpop)
            delete(hpop);
            %remove button
            delete(gcbo);
        end;
        return
        
    case 'delete color'
        really = questdlg('Really delete this color?','WARNING','yes','yes&recalc','NO','NO');
        switch really
            case {'NO',''}
                hpop = findall(0,'Style','listbox','Tag','selcetlabel');
                if ~isempty(hpop)
                    delete(hpop);
                    %remove button
                    delete(gcbo);
                end;
                return
            case 'yes'
                recalc = 0;
            otherwise
                recalc = 1;
        end %switch really
        
    case 'delete spot'
        recalc = 1;
        
end %switch what2del

cffig = openfig('labelgui','reuse');

imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');;
if isempty(imgFigureH)
    return;
end;

%clear popupmenu
hpop = findall(0,'Style','listbox','Tag','selcetlabel');
if ~isempty(hpop)
    delete(hpop);
end;
hpush = findall(0,'Style','pushbutton','Tag','deletespot');
if ~isempty(hpush)
    delete(hpush);
end;

% remove text
delete(get(gcbo,'UserData'));

%remove button
delete(gcbo);


timeslideH = findall(cffig,'Tag','slider3');
idlist = GetUserData(imgFigureH,'idlist');
tp = get(timeslideH,'Value');
rowIdx = find(idlist(tp).linklist(:,2)==spid);

switch what2del
    
    case 'delete color' %clear spot ( = clear color for all times)
        
        if length(rowIdx)>1
            h = warndlg('Sorry, You can''t delete multiple colors at the same time!','Warning');
            uiwait(h);
            return %end evaluation here
        end
        
        %find color to delete
        delColor = idlist(tp).linklist(rowIdx,4);
        
        %delete entry in list of labels
        labelIdx = log2(delColor)+1;
        idlist(1).stats.labelcolor(labelIdx)=[];
        
        t1 = 1;
        %delete color at all timepoints
        for t = 1:size(idlist,2)
            if ~isempty(idlist(t).linklist)
                %find row to delete
                delRow = find(idlist(t).linklist(:,4)==delColor);
                %spotNumber of color2delete
                delSpotNum = idlist(t).linklist(delRow,2);
                %intensity of color2delete
                delInt = idlist(t).linklist(delRow,8);
                
                    
                %delete row
                idlist(t).linklist(delRow,:) = [];
                
                %adjust all colors
                biggerColIdx = find(idlist(t).linklist(:,4)>delColor);
                idlist(t).linklist(biggerColIdx,4) = idlist(t).linklist(biggerColIdx,4)/2;
                
                %find if this deletes the whole spot
                spotRowIdx = find(idlist(t).linklist(:,2)==delSpotNum);
                if isempty(spotRowIdx) %only tag in spot -> adjust spotnumbers of all spots with higher number than the one in row spotRowIdx
                    biggerSpotNumIdx = find(idlist(t).linklist(:,2)>delSpotNum);
                    idlist(t).linklist(biggerSpotNumIdx,2) = idlist(t).linklist(biggerSpotNumIdx,2)-1;
                else %adjust intensity
                    idlist(t).linklist(spotRowIdx,3) = sum(idlist(t).linklist(spotRowIdx,4))*ones(size(spotRowIdx));
                    %intensity is distributed according to remainingInt/sumOfRemainingInts
                    sumInt = sum(idlist(t).linklist(spotRowIdx,8));
                    idlist(t).linklist(spotRowIdx,8) = idlist(t).linklist(spotRowIdx,8)*(1+delInt/sumInt);
                end
                
                
                %delete rows and cols in info-files
                if isfield(idlist,'info')
                    %we can just remove the tag from the Q-matrix: no tag
                    %changes its spot! (except for the one that moves to the no-spot...)
                    %info-matrices are 3-by-3 block matrices => rowIdx
                    %-> 3(rowIdx-1)+1:3rowIdx
                    tagNumber=log2(delColor)+1;
                    %delete rows
                    idlist(t).info.detectQ_Pix(3*(tagNumber-1)+1:3*tagNumber,:) = [];
                    % delete cols
                    idlist(t).info.detectQ_Pix(:,3*(tagNumber-1)+1:3*tagNumber) = [];
                    
                    %delete entry in info.noise (if-clause can be removed at a later date)
                    if isfield(idlist(tp).info,'noise')
                        idlist(t).info.noise = [];
                    end
                    
                    if ~isempty(idlist(tp).info.trackQ_Pix)
                        idlist(t).info.trackQ_Pix(3*(tagNumber-1)+1:3*tagNumber,:) = [];
                        idlist(t).info.trackQ_Pix(:,3*(tagNumber-1)+1:3*tagNumber) = [];
                    end
                end
                
                %redo spot color
                for i = 1:size(idlist(t).linklist,1)
                    sameIdx = find(idlist(t).linklist(:,2)==idlist(t).linklist(i,2));
                    idlist(t).linklist(i,3) = sum(idlist(t).linklist(sameIdx,4));
                end
                
                %adjust centroid
                idlist(t).centroid = mean(idlist(t).linklist(:,9:11),1);
                
                %sort idlist
                [idlist(t).linklist,sortIdx] = sortrows(idlist(t).linklist,4);
                
                %update Q-matrix
                detQ = [];
                nse = [];
                traQ = [];
                for j=1:size(idlist(t).linklist,1)
                    %the jth tag now was the ith tag before
                    i = sortIdx(j);
                    detQ = blkdiag(detQ,idlist(t).info.detectQ_Pix( (i-1)*3+1:i*3,(i-1)*3+1:i*3 ) );
                    %nse = [nse;idlist(curr_time).info.noise(i)];
                end
                if ~isempty(idlist(t).info.trackQ_Pix) 
                    for j=1:size(idlist(t).linklist,1)
                        i = sortIdx(j);
                        traQ = blkdiag(traQ,idlist(t).info.trackQ_Pix( (i-1)*3+1:i*3,(i-1)*3+1:i*3 ) );
                    end
                end
                idlist(t).info.detectQ_Pix=detQ;
                %idlist(curr_time).info.noise=nse;
                idlist(t).info.trackQ_Pix=traQ;
                
                %adjust linkup/linkdown
                if t>t1
                    %linkdown & update only if not is empty
                    if ~isempty(idlist(t1).linklist)
                        idlist(t1).linklist(:,7) = idlist(t).linklist(:,2);
                        %linkup
                        idlist(t).linklist(:,6) = idlist(t1).linklist(:,2);
                    end
                end  
                t1 = t;
            end
        end
        
        idlist(1).stats.maxColor = sum(idlist(tp).linklist(:,4))+1;
        
        %write idlist-status
        idlist(1).stats.status{end+1}=[date,': deleted color (in frame ',num2str(tp),')'];

        
    case 'delete spot' %delete spot spid
        %remember which tags were deleted
        tagNumber = log2(idlist(tp).linklist(rowIdx,4))+1; %list of tag numbers
        
        %delete linklist-entry
        idlist(tp).linklist(rowIdx,:) = [];
        %if there is no spot left at time tp, recalcIdlist will remove
        %the timepoint
        
        %adjust spot numbers if necessary
        biggerSpotNumIdx = find(idlist(tp).linklist(:,2)>spid);
        if ~isempty(biggerSpotNumIdx)
            idlist(tp).linklist(biggerSpotNumIdx,2) = idlist(tp).linklist(biggerSpotNumIdx,2)-1;
        end
        
        %delete row(s) in info-files (recalcIdlist does not read this per spot)
        %info-matrices are 3-by-3 block matrices => spid
        %-> 3(rowIdx-1)+1:3rowIdx
        if isfield(idlist,'info')
            %trackQ is per tag, not per spot, as is detectQ -> need to find tag "numbers"
            %first. There can be multiple tags deleted at the same time!
            %since this does not change where the other tags were 'sitting' 
            threeList = 1:3*max(tagNumber);
                threeMat = reshape(threeList,3,max(tagNumber)); %matrix [[1;2;3],[4;5;6]...]
                rmInfoIdx = threeMat(:,tagNumber); %all indices to remove rows and cols
                rmInfoIdx = rmInfoIdx(:); %to get a vector
                
            %delete rows
            idlist(tp).info.detectQ_Pix(rmInfoIdx,:) = [];
            % delete cols
            idlist(tp).info.detectQ_Pix(:,rmInfoIdx) = [];
            
            %delete info.noise (field is discontinued)
            if isfield(idlist(tp).info,'noise')
                idlist(tp).info.noise = [];
            end
            
            
            if ~isempty(idlist(tp).info.trackQ_Pix)
                
                idlist(tp).info.trackQ_Pix(rmInfoIdx,:) = [];
                idlist(tp).info.trackQ_Pix(:,rmInfoIdx) = [];
            end
        end
        
        %linkup/linkdown are not updated, idlist is passed to recalcIdlist
        
        %write idlist-status
        idlist(1).stats.status{end+1}=[date,': deleted spot ',num2str(spid),' (in frame ',num2str(tp),')'];
        
end %switch what2del




%if recalc: load dataProperties, recalc from beginning

if recalc
    dataProperties = GetUserData(imgFigureH,'dataProperties');
    idlist = recalcIdlist(idlist,1,[],dataProperties);
end



SetUserData(imgFigureH,idlist,1);

%get view3DH if exist and update data
view3DH = findall(0,'Tag','view3DGUI');
if ishandle(view3DH)
    view3D_generateHandles;
end

labelgui('refresh');





