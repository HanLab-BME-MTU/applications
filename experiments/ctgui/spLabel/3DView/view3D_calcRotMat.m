function view3D_calcRotMat
%calculates rotation matrices for align_axis

%load data
labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
view3DH = GetUserData(labelPanelH,'view3DGUIH');
view3D_handles = guidata(view3DH);
axisCoordinates = view3D_handles.axisCoordinates;
goodTime = view3D_handles.goodTime;

%init variables
tmax = size(axisCoordinates,3);
rotAngle = zeros(tmax,1);
rotVector = zeros(tmax,3);
rotMatrix = zeros(3,3,tmax);

%calculate axis vectors
axisVectors = axisCoordinates(2,:,:)-axisCoordinates(1,:,:);
axisVectors = squeeze(axisVectors)'; %makes it a tmax by 3 matrix

%norm axis vectors
[axisLength,n_axisVectors] = normList(axisVectors);

%Calculate rotMat for all timesteps with non-zero axis lengths
%rotate first one parallel to [1,0,0]
tStart = min(find(axisLength));
rotAngle(tStart) = acos(n_axisVectors(1,1)); %no need to calculate the dot-product with [1,0,0]
rotVector(tStart,:) = cross([1,0,0],n_axisVectors(tStart,:)); %first the direction we want to turn towards: rot3 rotates clockwise
rotMatrix(:,:,tStart) = rot3(rotVector(tStart,:),rotAngle(tStart),'rad');

t1 = tStart;
for t = tStart+1:max(goodTime) %no need to go looking further than that
    if axisLength(t)~=0
        %calculate rotMat
        rotAngle(t) = acos(n_axisVectors(t1,:)*n_axisVectors(t,:)');
        rotVector(t,:) = cross(n_axisVectors(t1,:),n_axisVectors(t,:));
        rotMatrix(:,:,t) = rotMatrix(:,:,t1)*rot3(rotVector(t,:),rotAngle(t),'rad'); %multiply from the right: first operation applied to the vector must be leftmost matrix
        %update counter
        t1 = t;
    end
end

% figure,arrow3(zeros(size(rotVector)),rotVector);
% title('rotation axes');

%fill in rotation matrices for good timepoints with zero-length axis

%find timepoints with zero-length axis
zeroLAtimes = goodTime(~axisLength(goodTime));
if ~isempty(zeroLAtimes)
    num2fillIn = length(zeroLAtimes);
    
    i = 1;
    time2fillIn = [];
    done = 0;
    
    while ~done
        time2fillIn = [time2fillIn;zeroLAtimes(i)];
        if i<num2fillIn 
            %check whether several consecutive times to fill in
            while any(zeroLAtimes==goodTime(find(zeroLAtimes(i)==goodTime)+1))
                i = i+1;
                time2fillIn = [time2fillIn;zeroLAtimes(i)];
            end
            if i==num2fillIn
                done = 1;
            end
        else
            done = 1;
        end
        
        LTFI = length(time2fillIn);
        
        %fill in rotation Matrices/vectors/axes
        switch (time2fillIn(1)==goodTime(1))+2*(time2fillIn(end)==goodTime(end))
            case 1 %no non-zero axis before: rotate like first non-zero axis (i.e. no rot-correction, just adjustment of coord-system)
                rT2 = goodTime(find(time2fillIn(end)==goodTime)+1); %first time with good rot-mat is one entry below in goodTime
                rotMatrix(:,:,time2fillIn) = repmat(rotMatrix(:,:,rT2),[1,1,LTFI]);
                rotVector(time2fillIn,:) = repmat(rotVector(rT2,:),[LTFI,1]);
                rotAngle(time2fillIn) = rotAngle(rT2);
            case 2 %no non-zero axis after: rotate like last non-zero axis
                rT1 = goodTime(find(time2fillIn(1)==goodTime)-1); %last time with good rot-mat is one entry above in goodTime
                rotMatrix(:,:,time2fillIn) = repmat(rotMatrix(:,:,rT1),[1,1,LTFI]);
                rotVector(time2fillIn,:) = repmat(rotVector(rT1,:),[LTFI,1]);
                rotAngle(time2fillIn) = rotAngle(rT1);
            case 0 %both non-zero axes before and after
                rT1 = goodTime(find(time2fillIn(1)==goodTime)-1); %last time with good rot-mat is one entry above in goodTime
                rT2 = goodTime(find(time2fillIn(end)==goodTime)+1); %first time with good rot-mat is one entry below in goodTime
                rotVector(time2fillIn,:) = repmat(rotVector(rT2,:),[LTFI,1]); %all have the same rot-vector (of t2)
                rotAnglePart = rotAngle(rT2)/(LTFI+1); %rotAngle is divided evenly among all timeSteps
                %calculate the first rotMat with rot3, others are just multiplications!
                rotAngle(time2fillIn(1)) = rotAnglePart;
                rotMatrix(:,:,time2fillIn(1)) = rot3(rotVector(time2fillIn(1),:),rotAnglePart,'rad')*rotMatrix(:,:,t1); 
                for j = 2:LTFI %for each: calulate rot-mat: r30*r30 = r60, r30*r60 = r90 etc
                    rotMatrix(:,:,time2fillIn(j)) = rotMatrix(:,:,time2fillIn(1))*rotMatrix(:,:,time2fillIn(j-1));
                    rotAngle(time2fillIn(1)) = rotAnglePart*j;
                end
            otherwise
                error('bad axis list (tags never seem to separate) or bug in view3D_calcRotMat or view3D_calcAxis')
        end
        i = i+1;
        time2fillIn = [];
    end %while-loop
end %if


%store data
view3D_handles.rotMatrix = rotMatrix;
view3D_handles.rotAngle = rotAngle;
view3D_handles.rotVector = rotVector;
view3D_handles.axisVectors = axisVectors;
view3D_handles.axisLength = axisLength;

guidata(view3DH,view3D_handles);