% function score_maps = network_similarity_with_breaks(sm_current_model,MT_current_model,img_size,radius)
% function for calculation the similarity of two networks
% Liya Ding 06.2013.


SM_Break = filament_model_to_seg_bwim(break_list,img_size,[]);

%% Breaks to be required
all_SM_current_seg = (SM_current_seg + SM_Break)>0;

SM_current_seg = all_SM_current_seg;

%% one pixel allowance, fat
SM_Break_fat = imdilate(SM_Break,ones(3,3));

fat_SM_current_seg = imdilate(SM_current_seg,ones(3,3));
fat_MT_current_seg = imdilate(MT_current_seg,ones(3,3));

%% total from SM -- ground truth, and from MT -- Test

Total_Actual_Positive = length(find(SM_current_seg==1));
Total_Actual_Negative = length(find(SM_current_seg==0));

Total_Test_Positive = length(find(MT_current_seg==1));
Total_Test_Negative = length(find(MT_current_seg==0));


%% The True False Positive Negative

Positive_map = SM_current_seg==1 & MT_current_seg ==1;
positive_break_map = SM_current_seg==0 & SM_Break == 1 & MT_current_seg ==1;

True_Posi = length(find(SM_current_seg>0 & MT_current_seg>0));
True_Neg = length(find(SM_current_seg==0 & MT_current_seg ==0));
False_Posi = length(find(SM_current_seg==0 & MT_current_seg ==1));

%% "fat" one pixel allowance measurments

False_Posi_fat = length(find(fat_SM_current_seg==0 & MT_current_seg ==1));

True_Posi_fat = length(find(fat_SM_current_seg>0 & MT_current_seg>0));

True_Neg = length(find(SM_current_seg==0 & MT_current_seg ==0));

False_Neg_fat = length(find(SM_current_seg==1 & fat_MT_current_seg ==0 & SM_Break_fat == 0));

%% Breaks to be considered
%"Allowing" is to accept break pixels, but not to penalize missed breaks

False_Posi_break = length(find(SM_current_seg==0 & MT_current_seg ==1 & SM_Break == 1));
True_Posi_break = length(find((Positive_map+positive_break_map)>0));
True_Posi_break = length(find(SM_current_seg==0 & SM_Break == 1 & MT_current_seg ==1));

Allowing_Positive_map = fat_SM_current_seg==1 & MT_current_seg ==1;

Allowing_True_Posi_break = length(find((Allowing_Positive_map+positive_break_map)>0));

Allowing_False_Posi = Total_Test_Positive - Allowing_True_Posi_break;

Allowing_false_negative_map = SM_current_seg==1 & fat_MT_current_seg ==0 & SM_Break == 0;

allowing_false_negative = length(Allowing_false_negative_map>0);

allowing_true_negative = Total_Test_Negative - allowing_false_negative;

False_Neg = length(find(SM_current_seg==1 & MT_current_seg ==0 & SM_Break == 0));









%% Scores

Accuracy_measure(iFrame)  = (True_Posi+True_Posi_break+True_Neg)/(True_Posi+True_Neg+False_Posi+False_Neg);

Precision(iFrame)  = (True_Posi+True_Posi_break)/(True_Posi+True_Posi_break+False_Posi);

Sensitivity(iFrame)  = (True_Posi+True_Posi_break)/(True_Posi+True_Posi_break+False_Neg);

F1Score(iFrame) = 2*(True_Posi+True_Posi_break)/(2*(True_Posi+True_Posi_break)+False_Neg+False_Posi);


%% Scores with the one pixel allowance, "fat"

Precision_fat(iFrame)  = (True_Posi_fat)/(True_Posi_fat+False_Posi_fat);

Sensitivity_fat(iFrame)  = (True_Posi_fat)/(True_Posi_fat+False_Neg_fat);

F1Score_fat(iFrame) = 2*(True_Posi_fat)/(2*(True_Posi_fat)+False_Neg_fat+False_Posi_fat);

%% Scores with breaks allowance "Allowing"


Allowing_Accuracy_measure(iFrame)  = (Allowing_True_Posi_break+allowing_true_negative)...
    /numel(SM_current_seg(:));

Allowing_Precision(iFrame)  = (Allowing_True_Posi_break)/Total_Test_Positive;

%% TPR and FPR for ROC curve

TPR(iFrame) = (True_Posi+True_Posi_break)/Total_Actual_Positive;
FPR(iFrame) = (False_Posi-False_Posi_break)/Total_Actual_Negative;

TPR_fat(iFrame) = (True_Posi_fat)/Total_Actual_Positive;
FPR_fat(iFrame) = (False_Posi_fat)/Total_Actual_Negative;

Allowing_TPR(iFrame) = (Allowing_True_Posi_break)/Total_Actual_Positive;
Allowing_FPR(iFrame) = (Allowing_False_Posi)/Total_Actual_Negative;


%% build a matched map display, including index; and find out the linked rate by seeing if the two ends are in the same filament index

if(~isempty(MT_current_model))
    model_length = length(MT_current_model);
    
    MT_seg_MT_BY_MT =cell(1,model_length);
    
    for iFila = 1 : model_length
        
        MT_seg_MT_BY_MT{iFila} = zeros(img_size);
        line_i_x = MT_current_model{iFila}(:,1);
        line_i_y = MT_current_model{iFila}(:,2);
        
        MT_seg_MT_BY_MT{iFila}(sub2ind(img_size, round(line_i_y),round(line_i_x)))=1;
    end
    
    
    %  two_channel_seg(:,:,3)=SM_Break;
    %
    %     h2=figure(2);imagesc(two_channel_seg);axis equal;axis off;
    
    Matched_flag = zeros(1,length(break_list));
    
    for iBreak = 1 : length(break_list)
        Begin_point  = break_list{iBreak}(1,:);
        begin_x = Begin_point(1);
        begin_y = Begin_point(2);
        
        
        End_point  =   break_list{iBreak}(end,:);
        end_x = End_point(1);
        end_y = End_point(2);
        
        if(begin_x==end_x && begin_y==end_y)
            Matched_flag(iBreak) = -1;
            
            for iFila = 1 : model_length
                if(sum(sum(MT_seg_MT_BY_MT{iFila}(begin_y-1:begin_y+1,begin_x-1:begin_x+1)))>=3)
                    %                 if(sum(sum(MT_seg_MT_BY_MT{iFila}(end_y,end_x)))>0)
                    Matched_flag(iBreak) = iFila;
                    %                 end
                end
            end
        else
            for iFila = 1 : model_length
                if(sum(sum(MT_seg_MT_BY_MT{iFila}(begin_y-1:begin_y+1,begin_x-1:begin_x+1)))>0)
                    if(sum(sum(MT_seg_MT_BY_MT{iFila}(end_y-1:end_y+1,end_x-1:end_x+1)))>0)
                        Matched_flag(iBreak) = iFila;
                    end
                end
            end
        end
    end
    linked_breaks_rate(iFrame) = mean(double(Matched_flag>0));
else
    linked_breaks_rate(iFrame) = NaN;
end


%% Canny 5*5

for iP= 1 : 5
    for iT = 1 :5
        
        MT_current_seg = current_seg_canny_cell{iP,iT};
        
        fat_MT_current_seg = imdilate(MT_current_seg,ones(3,3));
        
        %% total from SM -- ground truth, and from MT -- Test
        
        Total_Actual_Positive = length(find(SM_current_seg==1));
        Total_Actual_Negative = length(find(SM_current_seg==0));
        
        Total_Test_Positive = length(find(MT_current_seg==1));
        Total_Test_Negative = length(find(MT_current_seg==0));
        
        
        %% The True False Positive Negative
        
        Positive_map = SM_current_seg==1 & MT_current_seg ==1;
        positive_break_map = SM_current_seg==0 & SM_Break == 1 & MT_current_seg ==1;
        
        True_Posi = length(find(SM_current_seg>0 & MT_current_seg>0));
        True_Neg = length(find(SM_current_seg==0 & MT_current_seg ==0));
        False_Posi = length(find(SM_current_seg==0 & MT_current_seg ==1));
        
        %% "fat" one pixel allowance measurments
        
        False_Posi_fat = length(find(fat_SM_current_seg==0 & MT_current_seg ==1));
        
        True_Posi_fat = length(find(fat_SM_current_seg>0 & MT_current_seg>0));
        
        True_Neg = length(find(SM_current_seg==0 & MT_current_seg ==0));
        
        False_Neg_fat = length(find(SM_current_seg==1 & fat_MT_current_seg ==0 & SM_Break_fat == 0));
        
        %% Breaks to be considered
        %"Allowing" is to accept break pixels, but not to penalize missed breaks
        
        False_Posi_break = length(find(SM_current_seg==0 & MT_current_seg ==1 & SM_Break == 1));
        True_Posi_break = length(find((Positive_map+positive_break_map)>0));
        True_Posi_break = length(find(SM_current_seg==0 & SM_Break == 1 & MT_current_seg ==1));
        
        Allowing_Positive_map = fat_SM_current_seg==1 & MT_current_seg ==1;
        
        Allowing_False_Posi = Total_Test_Positive - Allowing_True_Posi_break;
        
        Allowing_false_negative_map = SM_current_seg==1 & fat_MT_current_seg ==0 & SM_Break == 0;
        
        allowing_false_negative = length(Allowing_false_negative_map>0);
        
        allowing_true_negative = Total_Test_Negative - allowing_false_negative;
        
        False_Neg = length(find(SM_current_seg==1 & MT_current_seg ==0 & SM_Break == 0));
        
        %% Scores
        
        Canny_Results{iP,iT}.Accuracy_measure(iFrame)  = (True_Posi+True_Posi_break+True_Neg)/(True_Posi+True_Neg+False_Posi+False_Neg);
        
        Canny_Results{iP,iT}.Precision(iFrame)  = (True_Posi+True_Posi_break)/(True_Posi+True_Posi_break+False_Posi);
        
        Canny_Results{iP,iT}.Sensitivity(iFrame)  = (True_Posi+True_Posi_break)/(True_Posi+True_Posi_break+False_Neg);
        
        Canny_Results{iP,iT}.F1Score(iFrame) = 2*(True_Posi+True_Posi_break)/(2*(True_Posi+True_Posi_break)+False_Neg+False_Posi);
        
        
        %% Scores with the one pixel allowance, "fat"
        
        Canny_Results{iP,iT}.Precision_fat(iFrame)  = (True_Posi_fat)/(True_Posi_fat+False_Posi_fat);
        
        Canny_Results{iP,iT}.Sensitivity_fat(iFrame)  = (True_Posi_fat)/(True_Posi_fat+False_Neg_fat);
        
        Canny_Results{iP,iT}.F1Score_fat(iFrame) = 2*(True_Posi_fat)/(2*(True_Posi_fat)+False_Neg_fat+False_Posi_fat);
        
        %% Scores with breaks allowance "Allowing"
        
        Allowing_True_Posi_break = length(find((Allowing_Positive_map+positive_break_map)>0));
        
        Canny_Results{iP,iT}.Allowing_Accuracy_measure(iFrame)  = (Allowing_True_Posi_break+allowing_true_negative)...
            /numel(SM_current_seg(:));
        
        Canny_Results{iP,iT}.Allowing_Precision(iFrame)  = (Allowing_True_Posi_break)/Total_Test_Positive;
        
        %% TPR and FPR for ROC curve
        
        Canny_Results{iP,iT}.TPR(iFrame) = (True_Posi+True_Posi_break)/Total_Actual_Positive;
        Canny_Results{iP,iT}.FPR(iFrame) = (False_Posi-False_Posi_break)/Total_Actual_Negative;
        
        Canny_Results{iP,iT}.TPR_fat(iFrame) = (True_Posi_fat)/Total_Actual_Positive;
        Canny_Results{iP,iT}.FPR_fat(iFrame) = (False_Posi_fat)/Total_Actual_Negative;
        
        Canny_Results{iP,iT}.Allowing_TPR(iFrame) = (Allowing_True_Posi_break)/Total_Actual_Positive;
        Canny_Results{iP,iT}.Allowing_FPR(iFrame) = (Allowing_False_Posi)/Total_Actual_Negative;
    end
end
