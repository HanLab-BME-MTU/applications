function [wrap_data] = frame_MIP_align(Prev_F, Current_F, Next_F,... % rescaled images
                                       Prev_F_org, Current_F_org, Next_F_org, ... % high res images
                                       rescale, MIP_params)
   
%% ================================================================================
% frame_MIP_align 
% -------------------------------------------------
% Find registratin between prev and next frame using a direct method:
% 1. find translation (T = T1*T2)
% 2. wrap next to prev & recode MIP.
% 3. recalculate affine transformation (H = H1*H2) using only the following pixels:
%  - Pixels with originally didn't have MIP
%  - Pixels with no MIP after translation T
%
% Inputs:
%  ~~~~~~
% Prev_F(_org) - The previous frame (rescaled and original)
% Current_F(_org) - The current frame (rescaled and original)
% Next_F(_org) - The next frame (rescaled and original)
% rescale - rescale factor 
% MIP_params - The parameters for MIP coding (params.MIP_params)
%
% Outputs:
% ~~~~~~~~
% wrap_data - struct that includes all wrapping information and results.
%
% Copyright 2012, Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
%
% ================================================================================

% Convert ORG images to Gray:
Prev_F_org = single(rgb2gray(Prev_F_org));
Current_F_org = single(rgb2gray(Current_F_org));
Next_F_org = single(rgb2gray(Next_F_org));
    
% Set default output struct:
wrap_data.MIP_H_p = eye(3);
wrap_data.MIP_H_n = eye(3);
wrap_data.Next_F_w = Next_F;
wrap_data.Prev_F_w = Prev_F;
wrap_data.valid_mask = zeros(size(Prev_F),'single');
    
% Get moving pixels for high res frames, use the silent ones for translation.
[MIP_mat] = calculate_MIP(Prev_F_org, Current_F_org, Next_F_org, MIP_params); 
org_silent_pix = sum(abs(MIP_mat),3) == 0;

% Parameters for align:
[w,h] = size(Current_F_org);
Lw = log2(w/MIP_params.align.min_size);
Lh = log2(h/MIP_params.align.min_size);
levels = floor(min(Lw,Lh))+1;

% 1) Find translation on all image pixels (no mask): 
% -------------------------------------------------------------
T1 = direct_align(Prev_F_org, Current_F_org, ... 
                                 levels, MIP_params.align.iterations, ...
                                 'trans',MIP_params.align.dbg_align_flag);                             
                             
T2 = direct_align(Current_F_org, Next_F_org, ...
                                 levels, MIP_params.align.iterations, ...
                                 'trans',MIP_params.align.dbg_align_flag);

if ( any(isnan(T1(:))) || any(isnan(T2(:))) )        
      % No valid transform, eliminate MIP of this frame.
      return;   
end
% 1) end ----------------------------------------------------------------------------------

% 2) Wrap next to prev with translation:
% ------------------------------------------------
T = inv(T2*T1);
if ( any(isnan(T(:))) || any(isinf(T(:))) )
     % No valid transform, eliminate MIP of this frame.
     return;        
end

% Wrap translation:
next_F_wrap = wrap_affine(Next_F_org, T);    
% Prepare valid mask (to ignore pixels on the margin)
prep_mask = zeros(size(Next_F_org),'single');
margin = MIP_params.align.valid_margin;
prep_mask(margin+1:end-margin-1,margin+1:end-margin-1) = 1;
valid_mask = wrap_affine(prep_mask, T);

% set output to this point:    
wrap_data.MIP_H_n = T;
wrap_data.Next_F_w = imresize(next_F_wrap, rescale); 

% 2) end ----------------------------------------------------------------------------------

% 3) Re-calculate MIP after translation to find trans_silent_pix 
% ----------------------------------------------------------------------------
% get MIP 8 trits per pixel        
[MIP_mat] = calculate_MIP(Prev_F_org, Current_F_org, next_F_wrap, MIP_params); 
trans_silent_pix = sum(abs(MIP_mat),3) == 0;
    
% Use one valid mask:
affine_valid_mask = org_silent_pix & trans_silent_pix & valid_mask;

% 3) end ----------------------------------------------------------------------------------

% 4) Find affine on non-silent pixels only
% -------------------------------------------------
% Find wrapping:
init_guess = T1;
H1 = direct_align( Prev_F_org, Current_F_org, ...
                  levels, MIP_params.align.iterations, ...
                  MIP_params.align.model_type, ...
                  MIP_params.align.dbg_align_flag, ...
                  init_guess, affine_valid_mask);

init_guess = T2;
H2 = direct_align( Current_F_org, Next_F_org, ...
                  levels, MIP_params.align.iterations, ...
                  MIP_params.align.model_type, ...
                  MIP_params.align.dbg_align_flag, ...
                  init_guess, affine_valid_mask);

% No valid transform:
% - take the translation T , set to struct.
% - next frame should already be wraped
% - update valid mask with moving pixels
if (any(isnan(H1(:))) || any(isnan(H2(:))))
    wrap_data.valid_mask = fix(imresize(valid_mask, rescale));
    wrap_data.valid_mask = imerode(wrap_data.valid_mask, ones(margin));
    return;
end
% 4) end ----------------------------------------------------------------------------------

% 5) Wrap next to prev
% ---------------------------
% Use one wrapping for next frame:
H = inv(H2*H1);

% No valid transform:
if ( any(isnan(H(:))) || any(isinf(H(:))) )
    wrap_data.valid_mask = fix(imresize(valid_mask, rescale));
    wrap_data.valid_mask = imerode(wrap_data.valid_mask, ones(margin));
    return;
end

% Warp the frame:
next_F_wrap = wrap_affine(Next_F_org, H);
next_F_wrap = imresize(next_F_wrap, rescale);

% Prepare valid mask with updated margins.
prep_mask = zeros(size(Next_F_org),'single');
margin = MIP_params.align.valid_margin;
prep_mask(margin+1:end-margin-1,margin+1:end-margin-1) = 1;
valid_mask = wrap_affine(prep_mask, H);
valid_mask = fix(imresize(valid_mask, rescale));
valid_mask = imerode(valid_mask,ones(margin));

% Update out Struct:
wrap_data.MIP_H_n = H;
wrap_data.Next_F_w = next_F_wrap;
wrap_data.valid_mask = valid_mask;
% 5) end ----------------------------------------------------------------------------------

end % of function: frame_MIP_align
