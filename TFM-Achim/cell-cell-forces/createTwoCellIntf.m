function constrForceField=createTwoCellIntf(constrForceField,frame,doPlot)
i=frame;
% first delete all existing twoCellIntf:
constrForceField{i}.twoCellIntf=[];

numItf=1;
for j=1:length(constrForceField{i}.interface)
    % intialization:
    m=1;
    ptsCurrIntf=length(constrForceField{i}.interface{j}.conct(:,1));
    pattern=constrForceField{i}.interface{j}.conct(m,:);
    
    % run through the whole length of the interface:
    while m<=ptsCurrIntf
        numPos=1;    
        % store the values for the first point, the starting vertex:
        twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
        twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);

        % determine which cells are connected through this interface. The
        % following procedure ensures that we get correct results if the
        % interface is composed of only two highly connected vertices:
        neigh1=find(constrForceField{i}.interface{j}.conct(m  ,:));
        neigh2=find(constrForceField{i}.interface{j}.conct(m+1,:));

        % The interface links the intersection of these two sets:
        twoCellIntf{numItf}.link=intersect(neigh1,neigh2);

        % now we can go to the next point of the interface:
        m=m+1;
        numPos=numPos+1;

        % while the connectivity of the points is <3 or we are not yet at the
        % end of this interface we go on:
        
        pattern=constrForceField{i}.interface{j}.conct(m,:);    
        while sum(pattern)<3 && m<=ptsCurrIntf            
            twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
            twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);
            m=m+1;
            numPos=numPos+1;
            if m<=ptsCurrIntf
                pattern=constrForceField{i}.interface{j}.conct(m,:);
            end
        end
        % check why the loop has terminated:
        if sum(pattern)>2 && m<ptsCurrIntf
            % we arrived at a higher order vertex within the interface:
            twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
            twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);
            
            % Don't increase m, the next interface should start at the same
            % vertex!
            numItf=numItf+1;
        elseif sum(pattern)>2 && m==ptsCurrIntf
            % we arrived at the end of the interface which is a higher 
            % order vertex. The position m=ptsCurrIntf has still to be
            % treated:
            twoCellIntf{numItf}.neigh{numPos,:}=find(pattern);
            twoCellIntf{numItf}.pos(numPos,:)  =constrForceField{i}.interface{j}.pos(m,:);
            
            % Now the interface has been run through completely, to 
            % terminate the loop:
            m=m+1;
            numItf=numItf+1;
        else
            % we arrived at the end of the interface which is a vertex of
            % connectivity<3. The interface has been run through
            % completely, to terminate the loop:
            m=m+1;
            numItf=numItf+1;
            % m is already at its maximum
        end        
    end
end
constrForceField{i}.twoCellIntf=twoCellIntf;

if nargin>2 && doPlot==1
    marker=['r','b','m','c','g','y','k'];
    figure(111)
    for p=1:length(twoCellIntf)
        plot(twoCellIntf{p}.pos(:,1),twoCellIntf{p}.pos(:,2),marker(mod(p,7)+1))
        hold on    
    end
    set(gca,'YDir','reverse')
    hold off
end