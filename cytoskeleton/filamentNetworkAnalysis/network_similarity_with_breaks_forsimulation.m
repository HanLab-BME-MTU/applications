% function score_maps = network_similarity_with_breaks(sm_current_model,MT_current_model,img_size,radius)
% function for calculation the similarity of two networks
% Liya Ding 06.2013.


SM_Break = filament_model_to_seg_bwim(break_list,img_size,[]);

Total_Actual_Positive = length(find(SM_current_seg==1));
Total_Actual_Negative = length(find(SM_current_seg==0));

True_Posi = length(find(SM_current_seg==1 & MT_current_seg ==1));
True_Neg = length(find(SM_current_seg==0 & MT_current_seg ==0));
False_Posi = length(find(SM_current_seg==0 & MT_current_seg ==1));
False_Neg = length(find(SM_current_seg==1 & MT_current_seg ==0 & SM_Break == 0));


Accuracy_measure(iFrame)  = (True_Posi+True_Neg)/(True_Posi+True_Neg+False_Posi+False_Neg);

Precision(iFrame)  = (True_Posi)/(True_Posi+False_Posi);

TPR(iFrame) = True_Posi/Total_Actual_Positive;
FPR(iFrame) = False_Posi/Total_Actual_Negative;

if(isempty(MT_current_model))
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

