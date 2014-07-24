function [fc1,fc2,cellInNv,cellAwNv]=assFcWithCells(constrForceField,frame,edgeId,fc,char)
edge=constrForceField{frame}.network.edge{edgeId};
pt1=round(edge.intf_internal(floor(end/2),:)+10*edge.nVec_internal);
pt2=round(edge.intf_internal(floor(end/2),:)-10*edge.nVec_internal);

nodes=edge.nodes;
% for cellId=nodes
%     figure(cellId)
%     imagesc(constrForceField{frame}.cell{cellId}.mask)
%     hold on
%     plot(pt1(1),pt1(2),'or')
%     plot(pt2(1),pt2(2),'xb')
% end

cellInNv=[];
cellAwNv=[];

% find the cells to which these points belong:
for cellId=nodes
    currMask=constrForceField{frame}.cell{cellId}.mask;
    [xMax,yMax]=size(currMask);
    if pt1(1)<=yMax && pt1(2)<=xMax && currMask(pt1(2),pt1(1))
        cellInNv=cellId;
    end
    
    if pt2(1)<=yMax && pt2(2)<=xMax && currMask(pt2(2),pt2(1))
        cellAwNv=cellId;
    end
end
% Cross validate with another method?!

% if we couldn't find the cell automatically, find it by hand. This happens
% very rarely:
if isempty(cellInNv) || isempty(cellAwNv)
    mask1=nodes(1)*constrForceField{frame}.cell{nodes(1)}.mask;
    mask2=nodes(2)*constrForceField{frame}.cell{nodes(2)}.mask;
    [rows1,cols1]=size(mask1);
    [rows2,cols2]=size(mask2);
    rowMax=max(rows1,rows2);
    colMax=max(cols1,cols2);
    if rows1<rowMax || cols1<colMax
        mask1(rowMax,colMax)=0;
    end
    if rows2<rowMax || cols2<colMax
        mask2(rowMax,colMax)=0;
    end    
    combMask=mask1+mask2;
    
    imagesc(combMask)
    colormap jet;
    colorbar;
    hold on;
    plot(pt1(1),pt1(2),'or');
    plot(pt2(1),pt2(2),'xb');
    hold off;
    center1=constrForceField{frame}.cell{nodes(1)}.center;
    center2=constrForceField{frame}.cell{nodes(2)}.center;
    text(center1(1),center1(2),num2str(nodes(1)));
    text(center2(1),center2(2),num2str(nodes(2)));
    title(['pt1=or; pt2=xb; frame= ',num2str(frame)]);
    if isempty(cellInNv)
        cellInNv=input('Enter cell ID of point pt1=or: ');
    end
    if isempty(cellAwNv)
        cellAwNv=input('Enter cell ID of point pt2=xb: ');
    end        
end



% check with the nodes connected to the edge:
if cellInNv==nodes(1) && cellAwNv==nodes(2)
    if char.status==-1
        fc1=-fc;
        fc2= fc;
    elseif char.status==1
        fc1= fc;
        fc2=-fc;
    elseif char.status==0
        % char.val =[ort_wghtd_mean ort_raw_mean];
        % the default is to take the weighted orientation:
        if char.val(1)<=0
            fc1=-fc;
            fc2= fc;
        else
            fc1= fc;
            fc2=-fc;
        end
    end
elseif cellInNv==nodes(2) && cellAwNv==nodes(1)
    if char.status==-1
        fc1= fc;
        fc2=-fc;
    elseif char.status==1
        fc1=-fc;
        fc2= fc;
    elseif char.status==0
        % char.val =[ort_wghtd_mean ort_raw_mean];
        % the default is to take the weighted orientation:
        if char.val(1)<=0
            fc1= fc;
            fc2=-fc;
        else
            fc1=-fc;
            fc2= fc;
        end
    end    
else
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    display('Couldnt associate cluster forces with cells')
    display('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    input('Report to Achim and hit [enter] to continue:...')
    cellInNv=NaN;
    cellAwNv=NaN;
    
    fc1= NaN;
    fc2= NaN;
end
    

 