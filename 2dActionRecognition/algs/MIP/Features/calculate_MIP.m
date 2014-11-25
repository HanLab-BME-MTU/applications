function MIP_mat = calculate_MIP(Perv_F, Current_F, Next_F, MIP_params)

%% ================================================================================
% calculate_MIP 
% -------------------------------------------------
% Calculate MIP votes for one frame, and one channel.
%
%
% Inputs:
% ~~~~~~
% Prev_F - The previous frame
% Current_F - The current frame
% Next_F - The next frame
% MIP_params - The parameters for MIP coding (params.MIP_params)
%
% Outputs:
%  ~~~~~~~
% MIP_mat - map of the size MxNx8 with 8 TRITS for each pixel 
%
% sunroutine:
%  ~~~~~~~~~~
% calculate_one_patch - MIP votes for specific patch position, for the whole image at once
% 
%
% Copyright 2012, Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
%
%% ================================================================================

% Calculate the 8 trits per pixel
patches_pos = [0 4; 3 3; 4 0; 3 -3; 0 -4; -3 -3; -4 0; -3 3];  % rounded 4 pixel circle
prev_pos = patches_pos;
next_pos = circshift(patches_pos,-1*MIP_params.MIP_code.code_shift);
TA = MIP_params.T;

% Calculate current alpha votes:
for patch = 1:8
    MIP_mat(:,:,patch) = calculate_one_patch(Perv_F, Current_F, Next_F, prev_pos(patch,:), next_pos(patch,:), TA);
end

% Switch patches position for suppression:
prev_pos = next_pos;
next_pos = patches_pos;
for patch = 1:8
    MIP_mat_reverse(:,:,patch) = calculate_one_patch(Perv_F, Current_F, Next_F, prev_pos(patch,:), next_pos(patch,:), TA);
end

% Suppress votes where : different signs and both not 0:
MIP_mat(MIP_mat ~= 0 & MIP_mat_reverse ~= 0  & MIP_mat ~=  MIP_mat_reverse ) = 0;
    
end % of function: calculate_MIP

%%
function one_patch_mat = calculate_one_patch(Perv_F, Current_F, Next_F, patch_pos_prev, patch_pos_next, TA)
%  MIP votes for one patch .

% Calculate the SSD per pixel with the shift of the specific patch position
SSD1 = ( circshift(Perv_F,patch_pos_prev) - Current_F).^2;
SSD2 = ( circshift(Next_F,patch_pos_next) - Current_F).^2; 

% Sum up the SSD of 9 pixels around
PA_SSD1 = 0; 
PA_SSD2 = 0;
for i = -1:1                 
    for j = -1:1
        PA_SSD1 = PA_SSD1 + circshift(SSD1, [i j]);
        PA_SSD2 = PA_SSD2 + circshift(SSD2, [i j]);
    end
end
% SSD sum diff:
PA_SSD = PA_SSD1 - PA_SSD2;

% Calculate where ||SSD1 - SSD2|| > treshold
pos_mov = double(PA_SSD > TA);
neg_mov = double(PA_SSD < -TA);

% Combine the +1 and the -1 to one map -1, 0, +1
one_patch_mat = pos_mov - neg_mov;

end % of  function: calculate_one_patch
