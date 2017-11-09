function [ coordsFN, number] = gcaValidationDocument(errorName,imSize,color)
% GCAValidationFiloBranchNetworkDocumentErrorTool
%

%% Initiate Counter
filoCount = 1;
clickFilo = 1;
coordsFN = [NaN,NaN];
ny = imSize(1); 
nx = imSize(2); 

% initiate false negative counter
hText =  text(20,ny-20, ['N ' errorName ' = 0 ']);
% while asking the user to click filo
while clickFilo == 1
    
    reply2 = questdlg(['Document ' errorName ' Segment?']);
    % if yes ask the user to choose a point
    if strcmpi(reply2,'yes')
        %         imshow(img,[])
        %         hold on
        %         cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
        %         spy(poleMaskTotal,'r',50); % plot the poles thus far ;
        h(filoCount)=impoint;
        position = wait(h(filoCount));
        
        if position(1) > nx || position(2) >ny
            clear position
            display('Please Click on the Left Image');
            delete(h(filoCount));
            h(filoCount) = impoint;
            position = wait(h(filoCount));
        end
        
        %idx = sub2ind(size(img),round(position(2)),round(position(1)));
        coordsFN(filoCount,2) = position(2);
        coordsFN(filoCount,1) = position(1);
        % setString(h,'False Negative');
        %  setColor(h,'r');
        
        % set(hText,'string','');
        %   total = nFilos + filoCount; % add the false negatives to the total
        % set(hText,'string',['Percent False Negatives = ' num2str((filoCount - 1)./total*100,3) ' for Frame ' ...
        %num2str(iFrame)]);
        set(hText,'string',['Number' errorName ' = ' num2str(filoCount) ]);
        %setString(h(i),'False Negative');
        % tag it on opposite side
        
        hfiloT(filoCount) = text(coordsFN(filoCount,1)+nx+2,coordsFN(filoCount,2),[errorName ' '  num2str(filoCount)],'FontSize',14);
        % setColor(h(i),'r');
        hfilo(filoCount) = scatter(coordsFN(filoCount,1)+nx,coordsFN(filoCount,2),30,color,'filled');
        filoCount = filoCount+1;
        % if no or cancel show the total
    elseif strcmpi(reply2,'no')
        %             % if no plot the filo
        %             if ~isnan(coordsFN(1,1))
        %                 for i = 1:size(h,2)
        %
        %
        %                 end
        %                 % ask if correct-
        %             end
        %
        % Ask if correct
        reply3 = questdlg(['Is This Final Number of ' errorName ' Correct?']);
        % if yes end while move to next
        if strcmpi(reply3,'yes')
            clear h
            clickFilo = 0;
            % if no reset the values
        elseif strcmpi(reply3,'no')
            reply4 = questdlg('How Would You like to Fix It?','','Add a Problem','Remove a Problem','Start Over','Add a Problem');
            switch reply4
                case 'Add a Problem';
                    
                    
                case 'Remove a Problem'
                    listFilo = num2cell(1:filoCount-1);
                    idxRemove  = listSelectGUI(listFilo,[],'move');
                    
                    
                    %                             delete(hfiloT(idxRemove));
                    %                             delete(hfilo(idxRemove));
                    
                    filoCount = filoCount-length(idxRemove)-1;
                    
                    test = zeros(size(coordsFN,1),1);
                    test(idxRemove) = 1;
                    test = logical(test);
                    % new coords
                    coordsFN = coordsFN(~test,:);
                    
                    % delete only the points on the first side
                    delete(h(idxRemove));
                    % delete everything else
                    delete(hfiloT);
                    delete(hfilo);
                    
                    if ~isempty(coordsFN)
                        % replot with new numbers
                        for i = 1:size(coordsFN,1)
                            hfiloT(i) = text(coordsFN(i,1)+nx+2,coordsFN(i,2),[errorName ' ' num2str(i)],'FontSize',14);
                            hfilo(i) = scatter(coordsFN(i,1)+nx,coordsFN(i,2),30,color,'filled');
                        end
                    end
                    set(hText,'string',['Number ' errorName ' = ' num2str(filoCount) ]);
                    filoCount = filoCount+1;
                    
                    
                case 'Start Over'
                    
                    
                    delete(hText)
                    if exist('hfilo','var');
                        delete(hfilo)
                    end
                    if exist('hfiloT','var'); 
                    delete(hfiloT)
                    end 
                    delete(h)
                    % initiate false false positive counter
                    hText =  text(20,ny-20, ['N ' errorName ' = 0']);
                    display(['Restarted ' errorName ' Negative Calculation']);
                    % if not or cancel restart counter
                    filoCount = 1;
                    coordsFN = [NaN,NaN];
                    
                    
            end
        elseif strcmpi(reply3,'cancel')
            display('Canceled By User: Exiting Without Saving');
            return
            
        end % strcmpi reply3
        
    elseif strcmpi(reply2,'cancel')
        display('Canceled By User: Exiting Without Saving');
        return
    end % strcmpi
    
end % while clickFilo
number = filoCount -1; 
if number>0 
      arrayfun(@(i) text(coordsFN(i,1)+2,coordsFN(i,2),[errorName ' ' num2str(i)],'FontSize',14),1:size(coordsFN(:,1)));
      arrayfun(@(i) scatter(coordsFN(i,1),coordsFN(i,2),30,color,'filled'),1:size(coordsFN(:,1)));
end 

end

