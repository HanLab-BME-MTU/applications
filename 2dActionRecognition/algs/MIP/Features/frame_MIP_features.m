function [points,frame_hist, MIP_for_display1, MIP_for_display2] = ...
    frame_MIP_features(f_id,Prev_F, Current_F, Next_F,val_map,MIP_params)

%% ================================================================================
% frame_MIP_features 
% -------------------------------------------------
% Calculate MIP features for one frame - one channel.
%
%
% Inputs:
% ~~~~~~
% f_id - frame id.
% Prev_F - The previous frame
% Current_F - The current frame
% Next_F - The next frame
% val_map - a map of binary encoded values with the size of the image.
% MIP_params - The parameters for MIP coding (params.MIP_params)
%
% Outputs:
% ~~~~~~~~
% points
% frame_hist
% MIP_for_display1
% MIP_for_display2
%
%
% Copyright 2012, Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
%
%% ================================================================================

hight = size(Current_F,1);
width = size(Current_F,2);
N = MIP_params.N;
M = MIP_params.M;
    
% Get MIP voting for moving pixels: (get MIP 8 trits per pixel)
MIP_mat = calculate_MIP(Prev_F, Current_F, Next_F, MIP_params); 
    
if (nargout>=3)
    % visualization
    MIP_for_display1 = sum(MIP_mat==1,3);
    MIP_for_display2 = sum(MIP_mat==-1,3);
end
    
% Lower bound on pixels in a patch (for data reduction):
patch_th = MIP_params.limit_patch_min.frac*N*M;% if have less pixels, remove the patch.
    
% Divide votes for positive and negatine maps:
[pos_mat,neg_mat] = MIP_val_map(MIP_mat,val_map);

% Adjust to valid in case valid_mask_en:
if(MIP_params.valid_mask_en)
    pos_mat = pos_mat.*MIP_params.valid_mask;
    neg_mat = neg_mat.*MIP_params.valid_mask;
    MIP_for_display1 =  MIP_for_display1.*MIP_params.valid_mask;
    MIP_for_display2 =  MIP_for_display2.*MIP_params.valid_mask;
end
    
% Vectorization : compose features vectors:
frame_hist = {};
f_ind = 1;
points = {};
for y = 6: N :hight - N
    for x = 6: M :width - M
        % find patch hist
        pos_B_mat = pos_mat( y: y + N -1, x: x + M -1,:);
        neg_B_mat = neg_mat( y: y + N -1, x: x + M -1,:);
        if (nnz(pos_B_mat(:)) ~= 0)
            pos_hist = histc(pos_B_mat(:),[0:255]);
            pos_hist(1) = 0;% 0 entry keep 0 always.
        else
            pos_hist = zeros(256,1,'single');
        end
        if(nnz(neg_B_mat(:)) ~= 0)
            neg_hist = histc(neg_B_mat(:),[0:255]);
            neg_hist(1) = 0;% 0 entry keep 0 always.
        else
            neg_hist = zeros(256,1,'single');
        end
        % Combine to one 512 vector:
        patch_hist = [pos_hist;neg_hist];
        % add non zero features only:
        if( (MIP_params.limit_patch_min.en && (sum(patch_hist) > patch_th))...
            || (~MIP_params.limit_patch_min.en && (nnz(patch_hist) > 0)) )
            frame_hist{f_ind} = patch_hist;
            points{f_ind} = [x+floor(N/2);y+floor(M/2);f_id];
            f_ind = f_ind+1;
        end
    end % on x
end % on y

% Convert to matrix:
frame_hist = cat(2,frame_hist{:});
points = cat(2,points{:}); 

end % of function: frame_MIP_features

%%
function [pos_mat,neg_mat] = MIP_val_map(MIP_B_mat,val_map)
    % fast implementation of map val
    pos_mat = (MIP_B_mat > 0);
    neg_mat = (MIP_B_mat < 0);
    pos_mat = sum(pos_mat.*val_map,3);
    neg_mat = sum(neg_mat.*val_map,3);
end % of function: MIP_val_map
