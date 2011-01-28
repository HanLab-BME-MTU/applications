function [fc1,fc2,cellInNv,cellAwNv]=assFcWithCells(constrForceField,frame,edgeId,fc,char)
edge=constrForceField{frame}.network.edge{edgeId};
pt1=round(edge.intf_internal(floor(end/2),:)+10*edge.nVec_internal);
pt2=round(edge.intf_internal(floor(end/2),:)-10*edge.nVec_internal);

nodes=edge.nodes;
for cellId=nodes
    figure(cellId)
    imagesc(constrForceField{frame}.cell{cellId}.mask)
    hold on
    plot(pt1(1),pt1(2),'or')
    plot(pt2(1),pt2(2),'xb')
end


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
        if ort_wghtd_mean<=0
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
        if ort_wghtd_mean<=0
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
    

 