%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the raw displacement field.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rawDispV = cell(numDTimePts,1);
%rawDataP = cell(numDTimePts,1);
if strcmp(trackMethod,'speck')
   %Load the displacement field given in the multidimentional array, 'M' that is
   % produced by the speckle tracking software.
   load(dataFile);
end

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   if strcmp(trackMethod,'speck')
      frame1 = imgIndexOfDTimePts(jj)-firstImgIndex+1;
      frame2 = frame1+numAvgFrames-1;
      rawDispField.p = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),2:-1:1,frame1);
      %rawDataP{jj} = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),2:-1:1,frame1);

      if strcmp(dispType,'MFAverage')
         rawDispField.v = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),4:-1:3,frame1);
         %rawDispV{jj}  = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),:,frame1);
         for k = frame1+1:frame2-1
            rawDispField.v = [rawDispField.v; M(find(M(:,1,k)~=0 & M(:,3,k)~=0),4:-1:3,k)];
            %rawDispV{jj} = [rawDispV{jj}; ...
            %   M(find(M(:,1,k)~=0 & M(:,3,k)~=0),:,k)];
         end
      elseif strcmp(dispType,'MFTrack')
         MFT = mFrameTrajBuild(M(:,:,frame1:frame2-1),MFCorLen, ...
            [rawDispField.p(:,2) rawDispField.p(:,1)]);
         rawDispField.v = MFT(:,end:-1:end-1) - MFT(:,2:-1:1);
         %MFT          = mFrameTrajBuild(M(:,:,frame1:frame2-1));
         %MFT = mFrameTrajBuild(M(:,:,frame1:frame2-1),MFCorLen, ...
         %   [rawDataP{jj}(:,2) rawDataP{jj}(:,1)]);
         %rawDispV{jj} = MFT(:,[1 2 end-1 end]);
      elseif strcmp(dispType,'SFrame')
         rawDispField.v = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),4:-1:3,frame1);
         %rawDispV{jj} = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),:,frame1);
      end
   elseif strcmp(trackMethod,'corr')
      s = load(dataFile{jj});
      flowTrack = s.flowTrack;

      numInd = find(~isnan(flowTrack.v{1}(:,1)) & ~isnan(flowTrack.v{1}(:,2)));
      rawDispField.p = flowTrack.p{1}(numInd,:);
      rawDispField.v = flowTrack.v{1}(numInd,:);
      %rawDataP{jj} = flowTrack.p{1}(numInd,:);
      %rawDispV{jj} = [flowTrack.p{1}(numInd,2:-1:1) ...
      %    flowTrack.p{1}(numInd,2:-1:1)+flowTrack.v{1}(numInd,2:-1:1)];
   elseif strcmp(trackMethod,'dArray')
      s = load(dataFile{jj});
      disp_array = s.disp_array;

      numInd = find(~isnan(disp_array(:,3)) & ~isnan(disp_array(:,4)));
      rawDispField.p = disp_array(numInd,1:2);
      rawDispField.v = disp_array(numInd,3:4);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Exclude outlier by fitting an ellipse to i) the distribution cloud of the two
   % components of actin velocity; 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %[outlier inlier exEllipse] = idMultiRandVarOutlier([rawDispV{1}(:,3:4)-rawDispV{1}(:,1:2) rawDispV{1}(:,1:2)]);
   %[outlier inlier exEllipse] = idMultiRandVarOutlier(rawDispV{1}(:,4:-1:3)-rawDispV{1}(:,2:-1:1));
   %[outlier inlier exEllipse] = idTwoRandVarOutlier(rawDispV{1}(:,3)-rawDispV{1}(:,1),rawDispV{1}(:,4)-rawDispV{1}(:,2));
   %[outlier inlier exEllipse] = idMultiRandVarOutlier(rawDispV{1}(:,4:-1:3)-rawDispV{1}(:,2:-1:1));
   %[outlier inlier exEllipse] = idTwoRandVarOutlier(rawDispV{1}(:,3)-rawDispV{1}(:,1),rawDispV{1}(:,4)-rawDispV{1}(:,2));

   %[outlier inlier exEllipse] = idMultiRandVarOutlier(rawDispField.v);
   %rawDispField.outlier = outlier;
   %rawDispField.inlier  = inlier;
   rawDispField.outlier = find(isnan(rawDispField.v(:,1)) | isnan(rawDispField.v(:,2)));
   rawDispField.inlier  = [1:size(rawDispField.v,1)];
   rawDispField.inlier(rawDispField.outlier) = [];

   imgIndex = imgIndexOfDTimePts(jj);
   rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
   save(rawDispFieldFile,'rawDispField');
   %save([mechDir filesep 'rawDispField'],'rawDispV','rawDataP','rawOutlier','rawInlier');
end

%drawFieldBound;

