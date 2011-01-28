function stats=TFM_part_7_statistics(constrForceField)
% we only need the field .network for this function and .segmRes.hole for
% calculating the intf_internal_L!
display('The interfacial length inlcudes only the part of the interface that is within the segmented cell perimeter!')
% In order to calculate the interfacial length, we take only every 10th
% point:
dPts=10;

load('fileAndFolderNames.mat')

if ~strcmp(pwd,path_ProjFolder)
    display('Before running this script browse to the FSM project folder')
    return
end

if nargin < 1 || isempty(constrForceField)
    display('Loading cellCellForces.mat. This file is large and thus takes some while:...')
    tic;
    filestruct=load(path_cellCellForces);
    constrForceField=filestruct.constrForceField;
    toc;    
end

% first check which frames have been analyzed:
% first test if the cluster analysis has been performed on this frame. For
% example, there is only a single cell:
toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'clusterAnalysis')
        toDoList=horzcat(toDoList,frame);
    end
end
%toDoList=10:31;

%**************************************************************************
% 1) check that both network and cluster forces agree.                    *
%**************************************************************************

alpha_nVec_fc=[];
alpha_fc_fnet=[];

% run through all frames that have been analyzed.
forceNum=0;
for frame=toDoList
    % run through all edges determined in that frame:
    for edgeNum=1:length(constrForceField{frame}.network.edge)
        % Check if the force could be determined by the network analysis:
        if ~isempty(constrForceField{frame}.network.edge{edgeNum}.f1) && ~isempty(constrForceField{frame}.network.edge{edgeNum}.f2) 
            % fc should be there since cluster Analysis has been performed
            % on that frame, thus also on this edge. Calculte the mean from
            % both forces:
            forceNum=forceNum+1;
            % calculate the mean vector:
            fnet{forceNum}.vec =1/2*(constrForceField{frame}.network.edge{edgeNum}.f1-constrForceField{frame}.network.edge{edgeNum}.f2);
            fnet{forceNum}.mag=sqrt(sum(fnet{forceNum}.vec.^2));
            
            fc{forceNum}.vec=constrForceField{frame}.network.edge{edgeNum}.fc;
            fc{forceNum}.mag=sqrt(sum(constrForceField{frame}.network.edge{edgeNum}.fc.^2));
            
            nVec{forceNum}.vec=constrForceField{frame}.network.edge{edgeNum}.n_Vec;
            nVec{forceNum}.mag=sqrt(sum(nVec{forceNum}.vec.^2));
            
            % to make this check correctly we first have to test which
            % normal has been picked when calculating the surface
            % stress/force in the cluster analysis:           
            alpha_nVec_fnet(forceNum)=acos(dot(fnet{forceNum}.vec,nVec{forceNum}.vec)/(fnet{forceNum}.mag*nVec{forceNum}.mag));
            if alpha_nVec_fnet<=pi/2
                signToPick=1;
            else
                signToPick=-1;
            end
            % that, by definition should also be ~180 degress.
            alpha_nVec_fc(forceNum)=acos(dot(fc{forceNum}.vec,nVec{forceNum}.vec)/(fc{forceNum}.mag*nVec{forceNum}.mag));
                        
            alpha_fc_fnet(forceNum)=pi-acos(dot(signToPick*fnet{forceNum}.vec,fc{forceNum}.vec)/(fnet{forceNum}.mag*fc{forceNum}.mag));
        end
    end
end
% is this worth storing?
% stats.fnet=fnet;
% stats.fc=fc;
% stats.nVec=nVec;
% stats.alpha_fc_fnet=alpha_fc_fnet;

figure(1)
for j=1:forceNum
    plot(fc{j}.mag,fnet{j}.mag,'or')
    hold on
end
xlabel('force from cluster analysis [nN]')
ylabel('force from network analysis [nN]')
hold off

if ~isempty(alpha_fc_fnet)
    figure(2)
    hist(alpha_fc_fnet*360/(2*pi))
    title('This 180-angle(fcluster,fnet)')
    xlabel('angular deviation [deg]')
    ylabel('frequency')
end

% figure(3)
% hist(alpha_nVec_fnet*360/(2*pi))
% title('angle normal vector and network force')
% xlabel('angular deviation [deg]')
% ylabel('frequency')
% 
% figure(4)
% hist(alpha_nVec_fc*360/(2*pi))
% title('This, by definition should also be ~180 degress.')
% xlabel('angular deviation [deg]')
% ylabel('frequency')


%**************************************************************************
% 2) Determine the residual force and elastic energy per connectivity     *
%**************************************************************************

% first get the maximal degree of connectivity:
maxDeg=0;
for frame=toDoList %10:31% toDoList
    % run through all nodes (cells) determined in that frame:
    for nodeNum=1:length(constrForceField{frame}.network.node)
        currDeg=constrForceField{frame}.network.node{nodeNum}.deg;
        if currDeg>maxDeg
           maxDeg=currDeg;
        end
    end
end

resforcePerDeg{maxDeg}    =[];
elEnergyPerDeg{maxDeg}    =[];
sumIntForcePerDeg{maxDeg} =[];
sumIntLengthPerDeg{maxDeg}=[];
% run through all frames:
for frame=toDoList %10:31% toDoList
    % run through all nodes (cells) determined in that frame:
    for nodeNum=1:length(constrForceField{frame}.network.node)
        deg=constrForceField{frame}.network.node{nodeNum}.deg;
        
        %******************************************************************
        % 1) determine residual force per connectivity:                   *
        %******************************************************************
        if isempty(resforcePerDeg{deg}) % if it is the first value to be sorted in:
            resforcePerDeg{deg}.val  = constrForceField{frame}.network.node{nodeNum}.mag;
            resforcePerDeg{deg}.frame= frame;
            % is the node number also interesting? That might change for
            % a particular cell between cells!
        else
            resforcePerDeg{deg}.val  = horzcat(resforcePerDeg{deg}.val,constrForceField{frame}.network.node{nodeNum}.mag);
            resforcePerDeg{deg}.frame= horzcat(resforcePerDeg{deg}.frame,frame);
        end
        
        %******************************************************************
        % 2) determine elastic energy per connectivity:                   *
        %******************************************************************
        if isempty(elEnergyPerDeg{deg}) % if it is the first value to be sorted in:
            elEnergyPerDeg{deg}.val  = constrForceField{frame}.network.node{nodeNum}.elE;
            elEnergyPerDeg{deg}.frame= frame;
            % is the node number also interesting? That might change for
            % a particular cell between cells!
        else
            elEnergyPerDeg{deg}.val  = horzcat(elEnergyPerDeg{deg}.val,constrForceField{frame}.network.node{nodeNum}.elE);
            elEnergyPerDeg{deg}.frame= horzcat(elEnergyPerDeg{deg}.frame,frame);
        end
        
        %******************************************************************
        % 3) determine the sum of all interfacial forces of cells with a  *
        %    certain connectivity. Do the same for the length of the      *
        %    interfaces.                                                  *
        %******************************************************************
        % First get all edges/interfaces connected to the current node: 
        edges=constrForceField{frame}.network.node{nodeNum}.edges;
        
        % sum up the forces transmitted through these edges/interfaces,
        % as well as the length of the interfaces:
        sumIntForce =0;
        sumIntLength=0;
        for edgeNum=edges
%!!!        Here we take the force determined through the cluster analysis,
%           This can be changed also to the forces determined through the
%           network analysis:            
            
            % first do the sum of forces:            
            sumIntForce =sumIntForce + norm(constrForceField{frame}.network.edge{edgeNum}.fc);
                       
            % use only the internal part of the interface to calculate the length:            
            pixSize_mu=constrForceField{frame}.par.pixSize_mu;
            if ~isfield(constrForceField{frame}.segmRes,'hole')
                % display('This is an old data set!')
                constrForceField{frame}.segmRes.hole=[];
            end
            currentLength=pixSize_mu*calcCurveLength(constrForceField{frame}.network.edge{edgeNum}.intf_internal,dPts,constrForceField{frame}.segmRes.hole);                        
            sumIntLength=sumIntLength+currentLength;
            % store these value only temporally in the constrForceField
            % structure:
            constrForceField{frame}.network.edge{edgeNum}.intf_internal_L=currentLength;            
        end
        if isempty(sumIntForcePerDeg{deg}) % if it is the first value to be sorted in:
            % first do the sum of forces:
            sumIntForcePerDeg{deg}.val  = sumIntForce;
            sumIntForcePerDeg{deg}.frame= frame;
            
            % Now the some of interfacial length:
            sumIntLengthPerDeg{deg}.val  = sumIntLength;
            sumIntLengthPerDeg{deg}.frame= frame;
            
            % is the node number also interesting? That might change for
            % a particular cell between cells!
        else
            % first do the sum of forces:
            sumIntForcePerDeg{deg}.val  = horzcat(sumIntForcePerDeg{deg}.val  ,sumIntForce);
            sumIntForcePerDeg{deg}.frame= horzcat(sumIntForcePerDeg{deg}.frame,frame);
            
            % Now the some of interfacial length:
            sumIntLengthPerDeg{deg}.val  = horzcat(sumIntLengthPerDeg{deg}.val  ,sumIntLength);
            sumIntLengthPerDeg{deg}.frame= horzcat(sumIntLengthPerDeg{deg}.frame,frame);
        end
    end
end
% Store all values:
stats.resforcePerDeg     = resforcePerDeg;
stats.elEnergyPerDeg     = elEnergyPerDeg;
stats.sumIntForcePerDeg  = sumIntForcePerDeg;
stats.sumIntLengthPerDeg = sumIntLengthPerDeg;

% marker for plotting
marker=['r','b','m','c','g','y','k'];

% Plot the residual forces:
% That seems to be different if the first images are skipped. Then there
% seems to be a difference:
figure(20)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg
    if ~isempty(resforcePerDeg{deg})
        groupVecVal=horzcat(groupVecVal,resforcePerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(resforcePerDeg{deg}.val)));
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title(['Residual force for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Residual force [nN]')


figure(21)
if ~isempty(resforcePerDeg{1})
    plot(resforcePerDeg{1}.frame,resforcePerDeg{1}.val,'or')
end
hold on
if ~isempty(resforcePerDeg{2})
    plot(resforcePerDeg{2}.frame,resforcePerDeg{2}.val,'*b')
end
title('Residual forces of cells with deg: 1 2')
xlabel('frame number')
ylabel('Residual force [nN]')
% legend('deg 1 cell','deg 2 cell')
hold off

% Plot the elastic energy:
figure(23)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg
    if ~isempty(elEnergyPerDeg{deg})
        groupVecVal=horzcat(groupVecVal,elEnergyPerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(elEnergyPerDeg{deg}.val)));
    end
end
if ~isempty(groupVecVal)
    boxplot(groupVecVal,groupVecDeg,'notch','on')
end
title(['Elastic energy for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Elastic energy [pJ]')

figure(24)
if ~isempty(elEnergyPerDeg{1}) && ~isempty(elEnergyPerDeg{1}.val)
    plot(elEnergyPerDeg{1}.frame,elEnergyPerDeg{1}.val,'or')
end
hold on
if ~isempty(elEnergyPerDeg{2}) && ~isempty(elEnergyPerDeg{2}.val)
    plot(elEnergyPerDeg{2}.frame,elEnergyPerDeg{2}.val,'*b')
end
title('Elastic energy of cells with deg: 1 2')
xlabel('frame number')
ylabel('Elastic energy [pJ]')
legend('deg 1 cell','deg 2 cell')
hold off


% Correlate residual forces with the elastic energy per degree of the
% node:
for deg=1:maxDeg 
    if ~isempty(elEnergyPerDeg{deg}) && ~isempty(elEnergyPerDeg{deg}.val)
        figure(200+deg)
        plot(elEnergyPerDeg{deg}.val,resforcePerDeg{deg}.val,'or')
        title(['Correlation residual force with elastic energy for cells with connectivity: ',num2str(deg)])
        xlabel('Elastic energy [pJ]')
        ylabel('residual Force [nN]')
    end
end

figure(200+maxDeg+1)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg
    if ~isempty(resforcePerDeg{deg}) && (~isempty(elEnergyPerDeg{deg})  && ~isempty(elEnergyPerDeg{deg}.val))
        groupVecVal=horzcat(groupVecVal,resforcePerDeg{deg}.val./elEnergyPerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(elEnergyPerDeg{deg}.val)));
    end
end
if ~isempty(groupVecVal)
    boxplot(groupVecVal,groupVecDeg,'notch','on')
end
title(['Residual Force per elastic energy for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Residual force / elastic energy [nN/pJ]')

% That seems to be significantly different:
figure(30)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg
    if ~isempty(sumIntForcePerDeg{deg})
        groupVecVal=horzcat(groupVecVal,sumIntForcePerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(sumIntForcePerDeg{deg}.val)));
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title(['Sum of all interfacial forces (magnitude) of a cell with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Sum of interfacial forces [nN]')


figure(31)
if ~isempty(sumIntForcePerDeg{1})
    plot(sumIntForcePerDeg{1}.frame,sumIntForcePerDeg{1}.val,'or')
end
hold on
if ~isempty(sumIntForcePerDeg{2})
    plot(sumIntForcePerDeg{2}.frame,sumIntForcePerDeg{2}.val,'*b')
end
title('Sum of interfacial forces of cells with deg: 1 2')
xlabel('Frame number')
ylabel('Sum of interfacial forces [nN]')
legend('deg 1 cell','deg 2 cell')
hold off

% Correlate sum of interfacial forces with the elastic energy per degree of
% the node:
for deg=1:maxDeg
    if ~isempty(elEnergyPerDeg{deg}) && ~isempty(elEnergyPerDeg{deg}.val)
        figure(300+deg)
        plot(elEnergyPerDeg{deg}.val,sumIntForcePerDeg{deg}.val,'or')
        title(['Correlation sum of interfacial forces with elastic energy for cells with connectivity: ',num2str(deg)])
        xlabel('Elastic energy [pJ]')
        ylabel('sum of interfacial forces [nN]')
    end
end

figure(300+maxDeg+1)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg
    if ~isempty(sumIntForcePerDeg{deg}) && (~isempty(elEnergyPerDeg{deg}) && (~isempty(elEnergyPerDeg{deg}.val)))
        groupVecVal=horzcat(groupVecVal,sumIntForcePerDeg{deg}.val./elEnergyPerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(sumIntForcePerDeg{deg}.val)));
    end
end
if ~isempty(groupVecVal)
    boxplot(groupVecVal,groupVecDeg,'notch','on')
end
title(['Sum of interfacial forces per elastic energy for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Sum of interfacial forces / elastic energy [nN/pJ]')


% That seems to be significantly different:
figure(32)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg 
    if ~isempty(sumIntLengthPerDeg{deg})
        groupVecVal=horzcat(groupVecVal,sumIntLengthPerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(sumIntLengthPerDeg{deg}.val)));
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title(['Interfacial length of a cell with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Cumulative length of all interfaces: [um]')


figure(33)
if ~isempty(sumIntLengthPerDeg{1})
    plot(sumIntLengthPerDeg{1}.frame,sumIntLengthPerDeg{1}.val,'or')
end
hold on
if ~isempty(sumIntLengthPerDeg{2})
    plot(sumIntLengthPerDeg{2}.frame,sumIntLengthPerDeg{2}.val,'*b')
end
title('Interfacial length of cells with deg: 1 2')
xlabel('Frame number')
ylabel('Cumulative length of all interfaces: [um]')
legend('deg 1 cell','deg 2 cell')
hold off

% Correlate sum of interfacial forces with the cumulative length of the
% interface:
figure(320)
for deg=1:maxDeg
    if ~isempty(sumIntLengthPerDeg{deg})
        plot(sumIntLengthPerDeg{deg}.val,sumIntForcePerDeg{deg}.val,['o',marker(mod(deg,7)+1)])
    end
    hold on;
    title(['Correlation sum of interfacial forces with cumulative interface length of cells with connectivity: ',num2str(deg)])
    xlabel('Cumulative length of all interfaces: [um]')
    ylabel('sum of interfacial forces [nN]')    
end
legend(num2str(1:maxDeg))
hold off

figure(320+maxDeg+1)
groupVecDeg=[];
groupVecVal=[];
for deg=1:maxDeg
    if ~isempty(sumIntForcePerDeg{deg})
        groupVecVal=horzcat(groupVecVal,sumIntForcePerDeg{deg}.val./sumIntLengthPerDeg{deg}.val);
        groupVecDeg=horzcat(groupVecDeg,deg*ones(size(sumIntForcePerDeg{deg}.val)));
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title(['Sum of interfacial forces per cumulative interfacial length for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Sum of interfacial forces / cumulative interfacial length [nN/um]')


%**************************************************************************
% 4) Determine the interfacial force per connectivity pair of the nodes   *
%    connected by this interface.                                         *
%    Then do the same for the interfacial stress, that is force normalized*
%    by the interfacial length.                                           *
%**************************************************************************

intForcePerDegPair{maxDeg,maxDeg} =[];
intLengthPerDegPair{maxDeg,maxDeg}=[];
intStressPerDegPair{maxDeg,maxDeg}=[];
% run through all frames:
for frame=toDoList
    % run through all edges (interfaces) determined in that frame:
    for edgeNum=1:length(constrForceField{frame}.network.edge)
        % get the two nodes that are connected by that edge/interface:
        nodes=constrForceField{frame}.network.edge{edgeNum}.nodes;
        
        % determine the conncectivity of these nodes:
        deg1=constrForceField{frame}.network.node{nodes(1)}.deg;
        deg2=constrForceField{frame}.network.node{nodes(2)}.deg;
        
        % order it such that deg1<deg2:
        if deg1>deg2
            storedDeg=deg2;
            deg2=deg1;
            deg1=storedDeg;
        end
%!!!    Here we take the force determined through the cluster analysis,
%       This can be changed also to the forces determined through the
%       network analysis:
            
        if isempty(intForcePerDegPair{deg1,deg2}) % if it is the first value to be sorted in:
            % 1) get the interfacial force:
            intForcePerDegPair{deg1,deg2}.val   = norm(constrForceField{frame}.network.edge{edgeNum}.fc);
            intForcePerDegPair{deg1,deg2}.frame = frame;
            
            % 2) get the interfacial length:
            intLengthPerDegPair{deg1,deg2}.val   = constrForceField{frame}.network.edge{edgeNum}.intf_internal_L;
            intLengthPerDegPair{deg1,deg2}.frame = frame;
            
            % 2) get the interfacial stress:
            intStressPerDegPair{deg1,deg2}.val   = intForcePerDegPair{deg1,deg2}.val./intLengthPerDegPair{deg1,deg2}.val;
            intStressPerDegPair{deg1,deg2}.frame = frame;
            
            % is the node number also interesting? That might change for
            % a particular cell between cells!
        else
            % 1) get the interfacial force:
            intForcePerDegPair{deg1,deg2}.val   = horzcat(intForcePerDegPair{deg1,deg2}.val  ,norm(constrForceField{frame}.network.edge{edgeNum}.fc));
            intForcePerDegPair{deg1,deg2}.frame = horzcat(intForcePerDegPair{deg1,deg2}.frame,frame);
            
            % 2) get the interfacial length:
            intLengthPerDegPair{deg1,deg2}.val   = horzcat(intLengthPerDegPair{deg1,deg2}.val  ,constrForceField{frame}.network.edge{edgeNum}.intf_internal_L);
            intLengthPerDegPair{deg1,deg2}.frame = horzcat(intLengthPerDegPair{deg1,deg2}.frame,frame);
            
            % 2) get the interfacial stress:
            intStressPerDegPair{deg1,deg2}.val   = intForcePerDegPair{deg1,deg2}.val./intLengthPerDegPair{deg1,deg2}.val;
            intStressPerDegPair{deg1,deg2}.frame = horzcat(intStressPerDegPair{deg1,deg2}.frame,frame);
        end
    end
end
% Store the values in the data structure:
stats.intForcePerDegPair =intForcePerDegPair;
stats.intLengthPerDegPair=intLengthPerDegPair;
stats.intStressPerDegPair=intStressPerDegPair;

% 1) plot the interfacial forces:
figure(40)
groupVecDeg=[];
groupVecVal=[];
for deg1=1:maxDeg
    for deg2=1:maxDeg
        if ~isempty(intForcePerDegPair{deg1,deg2})
            groupVecVal=horzcat(groupVecVal,intForcePerDegPair{deg1,deg2}.val);
%!!!        The following will go wrong if deg2>9 (which is however, very unlikely):
            if deg2<10
                groupVecDeg=horzcat(groupVecDeg,10*deg1*ones(size(intForcePerDegPair{deg1,deg2}.val))+deg2*ones(size(intForcePerDegPair{deg1,deg2}.val)));
            else
                display('The sorting in the boxplot might be wrong since deg2>=10!')
                break;
            end
        end
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title('Interfacial forces at edges connecting cells with connectivity 1,1; 1,2; 2,2')
xlabel('Pairs of connectivities')
ylabel('Interfacial force [nN]')


figure(41)
if ~isempty(intForcePerDegPair{1,1})
    plot(intForcePerDegPair{1,1}.frame,intForcePerDegPair{1,1}.val,'or')
    hold on
end
if length(intForcePerDegPair)>1
    if ~isempty(intForcePerDegPair{1,2})
        plot(intForcePerDegPair{1,2}.frame,intForcePerDegPair{1,2}.val,'*b')
        hold on
    end
    if ~isempty(intForcePerDegPair{2,2})
        plot(intForcePerDegPair{2,2}.frame,intForcePerDegPair{2,2}.val,'sg')
    end
end
legend('edge 1,1','edge 1,2','edge 2,2')
title('Interfacial forces at edges connecting cells with connectivity 1,1; 1,2; 2,2')
xlabel('Frame number')
ylabel('Interfacial force [nN]')
hold off


% 2) plot the interfacial length:
figure(42)
groupVecDeg=[];
groupVecVal=[];
for deg1=1:maxDeg
    for deg2=1:maxDeg
        if ~isempty(intLengthPerDegPair{deg1,deg2})
            groupVecVal=horzcat(groupVecVal,intLengthPerDegPair{deg1,deg2}.val);
%!!!        The following will go wrong if deg2>9 (which is however, very unlikely):
            if deg2<10
                groupVecDeg=horzcat(groupVecDeg,10*deg1*ones(size(intLengthPerDegPair{deg1,deg2}.val))+deg2*ones(size(intLengthPerDegPair{deg1,deg2}.val)));
            else
                display('The sorting in the boxplot might be wrong since deg2>=10!')
                break;
            end
        end
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title('Interfacial length of interfaces connecting cells with connectivity 1,1; 1,2; 2,2')
xlabel('Pairs of connectivities')
ylabel('Interfacial length [um]')


figure(43)
if ~isempty(intLengthPerDegPair{1,1})
    plot(intLengthPerDegPair{1,1}.frame,intLengthPerDegPair{1,1}.val,'or')
    hold on
end
if length(intLengthPerDegPair)>1
    if ~isempty(intLengthPerDegPair{1,2})
        plot(intLengthPerDegPair{1,2}.frame,intLengthPerDegPair{1,2}.val,'*b')
        hold on
    end
    if ~isempty(intLengthPerDegPair{2,2})
        plot(intLengthPerDegPair{2,2}.frame,intLengthPerDegPair{2,2}.val,'sg')
    end
end
legend('edge 1,1','edge 1,2','edge 2,2')
title('Length of interfaces connecting cells with connectivity 1,1; 1,2; 2,2')
xlabel('Frame number')
ylabel('Interfacial length [um]')
hold off



figure(44)
groupVecDeg=[];
groupVecVal=[];
for deg1=1:maxDeg
    for deg2=1:maxDeg
        if ~isempty(intStressPerDegPair{deg1,deg2})
            groupVecVal=horzcat(groupVecVal,intStressPerDegPair{deg1,deg2}.val);
%!!!        The following will go wrong if deg2>9 (which is however, very unlikely):
            if deg2<10
                groupVecDeg=horzcat(groupVecDeg,10*deg1*ones(size(intStressPerDegPair{deg1,deg2}.val))+deg2*ones(size(intStressPerDegPair{deg1,deg2}.val)));
            else
                display('The sorting in the boxplot might be wrong since deg2>=10!')
                break;
            end
        end
    end
end
boxplot(groupVecVal,groupVecDeg,'notch','on')
title('Stress at interfaces connecting cells with connectivity 1,1; 1,2; 2,2')
xlabel('Pairs of connectivities')
ylabel('Interfacial stress [nN/um]')


figure(45)
if ~isempty(intStressPerDegPair{1,1})
    plot(intStressPerDegPair{1,1}.frame,intStressPerDegPair{1,1}.val,'or')
    hold on
end
if length(intStressPerDegPair)>1
    if ~isempty(intStressPerDegPair{1,2})
        plot(intStressPerDegPair{1,2}.frame,intStressPerDegPair{1,2}.val,'*b')
        hold on
    end
    if ~isempty(intStressPerDegPair{2,2})
        plot(intStressPerDegPair{2,2}.frame,intStressPerDegPair{2,2}.val,'sg')
    end
end
legend('edge 1,1','edge 1,2','edge 2,2')
title('Stress at interfaces connecting cells with connectivity 1,1; 1,2; 2,2')
xlabel('Frame number')
ylabel('Interfacial stress [nN/um]')
hold off

return;

%**************************************************************************
% To be done:                                                             *
%**************************************************************************

% 1) Correlation analysis:
     


% 8) Calculate these values for cells with the same connectivity but for
%    different network topologies.

% 9) How to normalize the found values to compare between different
%    clusters? By the c1 or c11 values or by the elEnergy per cell?

