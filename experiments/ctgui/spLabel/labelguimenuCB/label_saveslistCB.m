function label_saveslistCB
%save to project data or to file or to 

%load handles and data
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    error('no movie loaded.');
    return;
end;
labelguiH = findall(0,'Tag','labelgui');
runCtBatchIsActive = GetUserData(labelguiH,'runCtBatchIsActive',1);

dataFile = GetUserData(imgFigureH,'dataFile');
idlist = GetUserData(imgFigureH,'idlist');
idname = GetUserData(imgFigureH,'idname');

%make set idname to idlist/idlisttrack
if strcmp(idname(end-1:end),'_L')
    idname = idname(1:end-2);
end

if runCtBatchIsActive %runCtBatch is active
    
    %set correct filename
    eval([idname,' = idlist;']);
    
    %save new idlist to file (runCtBatch gets data directly from labelPanel)
    save([dataFile.path,idname,'-L-',nowString],idname);
    
    uiresume(labelguiH); %let runCtBatch continue
    return
    
else %normal save
    if ~isempty(dataFile) %data from datafile
        
        %save as idlist_L / idlisttrack_L
        if ~strcmp(idname(end-1:end),'_L')
            idSaveName = [idname,'-L-'];
            idname = [idname,'_L'];
        else
            idSaveName = [idname(1:end-2),'-L-'];
        end
        
        eval([idname,' = idlist;']); %set correct name (idlist/idlisttrack)
        save([dataFile.path,dataFile.name],idname,'-append');
        
        %update project properties
        load([dataFile.path,dataFile.name],'projProperties');
        prevStatus = bsum2bvec(projProperties.status);
        if strcmp(idname,'idlist_L')
            projProperties.status = sum(prevStatus(find(prevStatus<8)))+8;
        else %it's an idlisttrack
            projProperties.status = sum(prevStatus(find(prevStatus<32)))+32;
        end
        %save project properties
        save([dataFile.path,dataFile.name],'projProperties','-append');
        
        %update lastResult
        load([dataFile.path,dataFile.name],'lastResult');
        lastResult = idname;
        save([dataFile.path,dataFile.name],'lastResult','-append');
        
        %save outside dataFile
        save([dataFile.path,idSaveName,nowString],idname);
        
    else
        if ~strcmp(idname(end-1:end),'_L')
            idname_L = [idname,'_L'];
            eval([idname_L,' = ',idname]);
            idSaveName = [idname,'-L-'];
        else
            idSaveName = [idname(1:end-2),'-L-'];
            idname_L = idname;
        end
        [fname,pathName] = uiputfile([idSaveName,nowString],'save idlist');
        if(fname(1)==0)
            image = [];
            fname = [];
            return;
        end;
        save([pathName fname],idname_L);
    end
end




%export spots to file and standard output: does not work with idlist (yet)
%fstat =  savespots([path fname],slist,3);

