function [ output_args ] = GCAValidationReconstructDocumentErrorFinal(reconstructDir)
% 



x = dir(reconstructDir);
x = x(3:end); 
list = arrayfun(@(i) x(i).name,1:length(x),'uniformoutput',0); 

% will eventually have a filter here that will check if done and only load
% those not finished if user desires. 


idxInclude  = listSelectGUI(list,[],'move');


% listSelectGUI to select a file to work on 
toTest = list(idxInclude); 





for iProj = 1: numel(toTest) 
    % make documentation directory 
    overlayDirC = [reconstructDir filesep toTest{iProj} filesep  'ValidationOverlays']; 

    if ~isdir(overlayDirC)
    mkdir(overlayDirC) 
    end 
    
    
   listOfImages =  searchFiles('.png',[],[reconstructDir filesep toTest{iProj} filesep 'Reconstruct_Movie'],0,'all',1); 
    imgRaw = imread(listOfImages{1}); 
    imgOverlay = imread(listOfImages{end}); 
    imgLarge = [imgOverlay imgRaw]; 
    [nyLarge,nxLarge,~] = size(imgLarge); 

% move to next non-documented frame? 

% choose a file manually? 

% 
% load the chosen project directory 


% The first and last frames in the Troubleshoot Reconstruction Directory
% will be what we load. 
    [ny,nx,~] = size(imgRaw);
   
  
    setFigure(nxLarge,nyLarge,'on');
    imshow(imgLarge,[]);
    hold on
    
  
    
    %% Initiate Counter
    filoCount = 1;
    clickFilo = 1;
    coordsFN = [NaN,NaN];
    
   
    % initiate false false positive counter
    hText =  text(20,20, ['N false positives = 0 ']);
    % while asking the user to click filo
    while clickFilo == 1
        
        reply2 = questdlg('Document False Negative Filopodia?');
        % if yes ask the user to choose a point
        if strcmpi(reply2,'yes')
            %         imshow(img,[])
            %         hold on
            %         cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
            %         spy(poleMaskTotal,'r',50); % plot the poles thus far ;
            h(filoCount)=impoint;
            position = wait(h(filoCount));
            %idx = sub2ind(size(img),round(position(2)),round(position(1)));
            coordsFN(filoCount,2) = position(2);
            coordsFN(filoCount,1) = position(1);
            % setString(h,'False Negative');
            %  setColor(h,'r');
            filoCount = filoCount+1;
            % set(hText,'string','');
            %   total = nFilos + filoCount; % add the false negatives to the total
            % set(hText,'string',['Percent False Negatives = ' num2str((filoCount - 1)./total*100,3) ' for Frame ' ...
            %num2str(iFrame)]);
            set(hText,'string',['Number False Negatives = ' num2str(filoCount - 1) ]);
            
            % if no or cancel show the total
        elseif strcmpi(reply2,'no')  || strcmpi(reply2,'cancel')
            % if no plot the filo
            if ~isnan(coordsFN(1,1))
                for i = 1:size(h,2)
                    
                    %setString(h(i),'False Negative');
                    % tag it on opposite side
                    hfiloT(i) = text(coordsFN(i,1)+nx+2,coordsFN(i,2),'False Neg');
                    % setColor(h(i),'r');
                    hfilo(i) = scatter(coordsFN(i,1)+nx,coordsFN(i,2),10,'b','filled');
                end
                % ask if correct-
            end
            
            % Ask if correct
            reply3 = questdlg('Is This Final Number of False Negatives Correct?');
            % if yes end while move to next
            if strcmpi(reply3,'yes')
                clear h
                clickFilo = 0;
                % if no reset the values
            elseif strcmpi(reply3,'no') || strcmpi(reply3,'cancel')
                reply4 = questdlg('How Would You like to Fix It?','','Add a Filopodia','Remove a Filopodia','Start Over','Add a Filopodia');
                switch reply4
                    case 'Add a Filopodia';
                        
                        
                    case 'Remove a Filopodia'
                        
                    case 'Start Over'
                        
                        
                        delete(hText)
                        delete(hfilo)
                        delete(hfiloT)
                        delete(h)
                        % initiate false false positive counter
                        hText =  text(20,20, ['N false positives = 0']);
                        display(['Restarted False Negative Calculation']);
                        % if not or cancel restart counter
                        filoCount = 1;
                        coordsFN = [NaN,NaN];
                end
            end
            
            
        end % strcmpi
              
    end % while
    % save every frame
    %filoValidation.coordFN = coordsFN;
  
    save('coordsFN.mat','coordsFN'); 
    saveas(gcf,[overlayDirC filesep 'FalseNegativeOverlay.png']); 
%% End 


    
   
    setFigure(nxLarge,nyLarge,'on');
    % setFigure(nxLarge,nyLarge,'on');
    imshow(imgLarge,[]);
    hold on
    
    
    
    
   
    %%%%%%%%%%%%%%%%%% Initiate Counter  %%%%%%%%%%%%%%%%%%%
    filoCount = 1;
    clickFilo = 1;
    coordsFP = [NaN,NaN];
    
   
    
    % initiate false false positive counter
    hText =  text(20,20, ['N FALSE POSITIVES = 0 ']);
    % while asking the user to click filo
    while clickFilo == 1
        
        reply2 = questdlg('Document False Positive Filopodia?');
        % if yes ask the user to choose a point
        if strcmpi(reply2,'yes')
            %         imshow(img,[])
            %         hold on
            %         cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
            %         spy(poleMaskTotal,'r',50); % plot the poles thus far ;
            h(filoCount)=impoint;
            position = wait(h(filoCount));
            %idx = sub2ind(size(img),round(position(2)),round(position(1)));
            coordsFP(filoCount,2) = position(2);
            coordsFP(filoCount,1) = position(1);
            % setString(h,'False Negative');
            %  setColor(h,'r');
            filoCount = filoCount+1;
            % set(hText,'string','');
            %   total = nFilos + filoCount; % add the false negatives to the total
            % set(hText,'string',['Percent False Negatives = ' num2str((filoCount - 1)./total*100,3) ' for Frame ' ...
            %num2str(iFrame)]);
            set(hText,'string',['Number False Positives = ' num2str(filoCount - 1)] );
            
            % if no or cancel show the total
        elseif strcmpi(reply2,'no')  || strcmpi(reply2,'cancel')
            % if no plot the filo
            if ~isnan(coordsFP(1,1))
                for i = 1:length(h)
                    
                    %setString(h(i),'False Negative');
                    % tag it on opposite side
                    text(coordsFP(i,1)+nx+2,coordsFP(i,2),'False Positive');
                    % setColor(h(i),'r');
                    scatter(coordsFP(i,1)+nx,coordsFP(i,2),10,'r','filled')
                end
                % ask if correct-
            end
            
            % Ask if correct
            reply3 = questdlg('Is This Final Number of False Positives Correct?');
            % if yes end while move to next
            if strcmpi(reply3,'yes')
                delete(hText)
                clear h
                clickFilo = 0;
                % if no reset the values
            elseif strcmpi(reply3,'no') || strcmpi(reply3,'cancel')
                close gcf
                display('Restarted False Positive Calculation ');
                % if not or cancel restart counter
                filoCount = 1;
                coordsFP = [NaN,NaN];
            end
            
            
        end % strcmpi
        
        
        
        
        
    end % while
     


    




%%
%% Document Misconnections
%for iFrame = 1:movieData.nFrames_ -1
    
    
    imgRaw = imread(listOfFilesRaw{iFrame});
    imgOverlay =imread(listOfFilesOverlay{iFrame});
    [ny,nx,~] = size(imgRaw);
    imgLarge = [imgOverlay imgRaw];
    [nyLarge,nxLarge,~] = size(imgLarge);
    setFigure(nxLarge,nyLarge,'on');
    % setFigure(nxLarge,nyLarge,'on');
    imshow(imgLarge,[]);
    hold on
    
    %     filoInfoC = analInfo(iFrame).filoInfo;
    %
    %     % load the raw images
    %     img =  double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}]));
    %     [ny,nx] = size(img);
    %     %subplot(1,2,1);
    %
    %     imgLarge = [img img];
    %     [nyLarge,nxLarge] = size(imgLarge);
    %
    %     setFigure(nxLarge,nyLarge,'on');
    %     imshow(-imgLarge,[]);
    %     hold on
    %
    %     n = length(filoInfoC);
    %     c = linspecer(n);
    %     idxRand = randperm(n);
    %     c = c(idxRand,:);
    %     display(['Loading Overlays for Frame ' num2str(iFrame)]);
    %     for ifilo = 1:length(filoInfoC)
    %         filoInfoIdx = filoInfoC(ifilo);
    %         GCAVisualsMakeOverlaysFilopodia(filoInfoIdx,[ny,nx],1,1,c(ifilo,:),0);
    %         clear filoInfoIdx
    %     end
    %%%%%%%%%%%%%%%%%% Initiate Counter  %%%%%%%%%%%%%%%%%%%
    filoCount = 1;
    clickFilo = 1;
    coordsMis = [NaN,NaN];
    
    %load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_1' filesep ...
    %   'analInfoTestSave.mat']);
    % filoInfo = analInfo(iFrame).filoInfo;
    
    % initiate false false positive counter
    hText =  text(20,20, ['N Misconnections = 0 for Frame ' num2str(iFrame)]);
    % while asking the user to click filo
    while clickFilo == 1
        
        reply2 = questdlg('Document Misconnected Filopodia ?');
        % if yes ask the user to choose a point
        if strcmpi(reply2,'yes')
            %         imshow(img,[])
            %         hold on
            %         cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
            %         spy(poleMaskTotal,'r',50); % plot the poles thus far ;
            h(filoCount)=impoint;
            position = wait(h(filoCount));
            %idx = sub2ind(size(img),round(position(2)),round(position(1)));
            coordsMis(filoCount,2) = position(2);
            coordsMis(filoCount,1) = position(1);
            % setString(h,'False Negative');
            %  setColor(h,'r');
            filoCount = filoCount+1;
            % set(hText,'string','');
            %   total = nFilos + filoCount; % add the false negatives to the total
            % set(hText,'string',['Percent False Negatives = ' num2str((filoCount - 1)./total*100,3) ' for Frame ' ...
            %num2str(iFrame)]);
            set(hText,'string',['Number Filopodia Misconnected = ' num2str(filoCount - 1) ' for Frame ' ...
                num2str(iFrame)]);
            
            % if no or cancel show the total
        elseif strcmpi(reply2,'no')  || strcmpi(reply2,'cancel')
            % if no plot the filo
            if ~isnan(coordsMis(1,1))
                for i = 1:length(h)
                    
                    %setString(h(i),'False Negative');
                    % tag it on opposite side
                    text(coordsMis(i,1)+nx+2,coordsMis(i,2),'Misconnections');
                    % setColor(h(i),'r');
                    scatter(coordsMis(i,1)+nx,coordsMis(i,2),10,'g','filled')
                end
                % ask if correct-
            end
            
            % Ask if correct
            reply3 = questdlg('Is This Final Number of Misconnections Correct?');
            % if yes end while move to next
            if strcmpi(reply3,'yes')
                delete(hText)
                clear h
                clickFilo = 0;
                
                % if no reset the values
            elseif strcmpi(reply3,'no') || strcmpi(reply3,'cancel')
                close gcf
                display(['Restarted Misconnection Calculation for ' num2str(iFrame)]);
                % if not or cancel restart counter
                filoCount = 1;
                coordsMis = [NaN,NaN];
            end
            
            
        end % strcmpi
        
         %Record Misconnections in a structure
    filoValidation(iFrame).coordsMis = coordsMis;
    saveDir = [movieData.outputDirectory_ filesep 'ValidationDirectory' filesep 'coordsMis'];
    if ~isdir(saveDir) 
        mkdir(saveDir) 
    end 
    saveas(gcf,[saveDir filesep num2str(iFrame,'%03d') '.png']);
        
        
        
   
    
    
    end % while
    
    
    
    
%     reply4 = questdlg('Move To Next Frame?') ;
%     if strcmpi(reply4,'yes');
%         close gcf
%     else
%         close gcf
%         break
%     end
   
    
%end % iframe



%save([movieData.outputDirectory_ filesep 'filoValidation'],'filoValidation');

end % iProj 



