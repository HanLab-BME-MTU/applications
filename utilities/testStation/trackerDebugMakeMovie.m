function trackerDebugMakeMovie(sourceInfo,targetInfo)


    % make 3d images (arrange 6x2) of: Source, target, delta, dx, dy, dz
    [s1, s2, t1, t2, d1, d2, dx1, dx2, dy1, dy2, dz1, dz2] = deal(...
        repmat(NaN,sourceInfo(1).size));
    s1(sourceInfo(1).goodIdx) = sourceInfo(1).intList(sourceInfo(1).goodIdx);
    s2(sourceInfo(2).goodIdx) = sourceInfo(2).intList(sourceInfo(2).goodIdx);
    t1(sourceInfo(1).goodIdx) = targetInfo(1).intList(sourceInfo(1).goodIdx);
    t2(sourceInfo(1).goodIdx) = targetInfo(2).intList(sourceInfo(1).goodIdx);
    d1(sourceInfo(1).goodIdx) = targetInfo(1).deltaInt(sourceInfo(1).goodIdx);
    d2(sourceInfo(1).goodIdx) = targetInfo(2).deltaInt(sourceInfo(1).goodIdx);
    dx1(sourceInfo(1).goodIdx) = targetInfo(1).gradient(:,1);
    dx2(sourceInfo(1).goodIdx) = targetInfo(2).gradient(:,1);
    dy1(sourceInfo(1).goodIdx) = targetInfo(1).gradient(:,2);
    dy2(sourceInfo(1).goodIdx) = targetInfo(2).gradient(:,2);
    dz1(sourceInfo(1).goodIdx) = targetInfo(1).gradient(:,3);
    dz2(sourceInfo(1).goodIdx) = targetInfo(2).gradient(:,3);
    
    i1=cat(1,s1,t1,d1,dx1,dy1,dz1);
    i2=cat(1,s2,t2,d2,dx2,dy2,dz2);
    
    image = cat(2,i1,i2);
    ct = evalin('base','ct');
    if ~evalin('base','exist([''control'',num2str(ct)],''var'')')
        
        evalin('base',['control',num2str(ct),'={};']);
    end
    assignin('base','tmp',image);
    evalin('base',['control',num2str(ct),'{end+1}=tmp;']);
    
