function [shift,transform,A,B] = driftMarkerRegistration(drift,drift_ref, ImSize)
% driftMarkerRegistration calculated interitively first through translation
% only (image based) and then uses procrustus to determine the rotation and translational
% shift. 
%
% Inputs: 
%         drift, n x 2 vector with locations of the drift markers to be
%                shifted.
%     drift_ref, m x 2 vector with drift marker locations in reference
%                image with m >= n
%        ImSize, image size of orginal image [npix,npix] must be square
%
% Outputs:
%         shift, the cumulative translation of points to be applied before
%                applying general rigid transform
%     Transform, transform structure from procrustus plus a term .shift
%                which is the best translation only alignment
%
% Written by Jeffrey Werbin
% HMS
%           
    
    s_d = size(drift);
    s_ref = size(drift_ref);
    
    if s_d(1) > s_ref(1)
        error('the reference must have an equal or larger number of drift markers');
    end
    
    %Create images of each point list with each point represented as a
    %gaussian with a sigma of 1 pix
    
    %shift data to remove negative numbers
    bump= min(vertcat(drift,drift_ref),[],1);
    bump(bump > 0) = 0;
    
    [img]=PointP_Plot2([drift-repmat(bump,[s_d(1),1]),ones(size(drift))],ImSize,1,10,'test',5);
    d_img = img(:,:,2);
    [img]=PointP_Plot2([drift_ref-repmat(bump,[s_ref(1),1]),ones(size(drift_ref))],ImSize,1,10,'test',5);
    ref_img = img(:,:,2);
    
    C = xcorr2(d_img,ref_img);
    mC = max(C(:));
    [x,y] = ind2sub(size(C),find(C == mC));
    shift = [x,y];
    
    shift = shift - size(ref_img);
    
    %apply shift to drift
    drift = drift - repmat(shift,[s_d(1),1]);
    
    %Calculate general transform
    %
    % first points in drift are matched to points in drift_ref
     
    dm = distMat2(drift,drift_ref);
    
    % Assigns drift markers in ref to markers in test
    [OneToTwo,TwoToOne]=lap(dm,-1,0,1);
    
    %Finds the transform that reconciles the two sets of points
    % is dependent on absor code
    tmp = OneToTwo(1:s_d(1));
    A = drift_ref(tmp,:);
    B = drift;
    
    %Z = TRANSFORM.b * Y * TRANSFORM.T + TRANSFORM.c.
    [d, Z, transform] = procrustes(A, B,'Reflection',false,'Scaling',false);
    trans = mean(A-B);
    
    transform.trans = trans;
    
end
    
    
%    shift = [0,0];    
%    % make a working copy of the inputs
%    drift_start = drift;     
% 
%     criteria = 1;
%     last_lap = zeros([s_d(1),1]);
%     Itt = 0;
%     
%     while criteria 
% 
%         dm = distMat2(drift,drift_ref);
% 
%         % Assigns drift markers in ref to markers in test
%         [OneToTwo,TwoToOne]=lap(dm,-1,0,1);
% 
%         %Finds the transform that reconciles the two sets of points
%         % is dependent on absor code
%         tmp = OneToTwo(1:s_d(1));
%         A = drift_ref(tmp,:);
%         B = drift;
%         trans = mean(A-B);
%         
%         %determine if another iteration is needed
%         criteria = any(tmp ~=last_lap);
%         last_lap = tmp;
%         shift = shift + trans;
%         
%         drift = drift + shift;
%         Itt = Itt+1;
%         
%     end
    
%    disp(['required ',num2str(Itt),' iterations']);