function output_feature_reorganized = organize_profiles(output_feature,feature_index)
                           
output_feature_reorganized = output_feature;
T_dis_perp =40;
T_center = 10;

if feature_index(19) == 1
       
    %19
    output_feature_reorganized.nmssum_ratio_pool = output_feature.profileAllCell.nmssum_pool(:);
    %20
    output_feature_reorganized.nmsmean_ratio_pool = output_feature.profileAllCell.nmsmean_pool(:);
    %21
    output_feature_reorganized.intsum_ratio_pool = output_feature.profileAllCell.intmean_pool(:);
    %22
    output_feature_reorganized.intmean_ratio_pool = output_feature.profileAllCell.intmean_pool(:);
    %23
    output_feature_reorganized.filasum_ratio_pool = output_feature.profileAllCell.filasum_pool(:);
    %24
    output_feature_reorganized.filamean_ratio_pool = output_feature.profileAllCell.filamean_pool(:);
    
    
    %% automatic determing Distance_periphery_center
    
    NoRegion = numel(output_feature.profileCell);
    profileIntSumPerpCenterRatio = nan(NoRegion,6);
    profileIntMeanPerpCenterRatio = nan(NoRegion,6);
    profileFilaSumPerpCenterRatio = nan(NoRegion,6);
    profileFilaMeanPerpCenterRatio = nan(NoRegion,6);
    profileStnmsSumPerpCenterRatio = nan(NoRegion,6);
    profileStnmsMeanPerpCenterRatio = nan(NoRegion,6);
        
        DistancePeriCenter = nan(NoRegion,6);
        
    for  iR = 1 : NoRegion
        
        
        profileIntSum    = output_feature.profileCell{1,iR}.profileIntSum;
        profileFilaSum   = output_feature.profileCell{1,iR}.profileFilaSum;
        profileStnmsSum  = output_feature.profileCell{1,iR}.profileStnmsSum;
        
        profilePixels    = output_feature.profileCell{1,iR}.profilePixels;
        
        profileIntMean   = output_feature.profileCell{1,iR}.profileIntSum./(output_feature.profileCell{1,iR}.profilePixels);
        profileFilaMean  = output_feature.profileCell{1,iR}.profileFilaSum./(output_feature.profileCell{1,iR}.profilePixels);
        profileStnmsMean = output_feature.profileCell{1,iR}.profileStnmsSum./(output_feature.profileCell{1,iR}.profilePixels);
        
        
        
        for iAngle =  1 : 6
            
            profileIntAve_array = profileIntMean(iAngle, :);
            profileSTAve_array = profileStnmsMean(iAngle, :);
            profileFilaAve_array = profileFilaMean(iAngle, :);
            
                        
            profileIntAve_array_smooth = imfilter(profileIntAve_array, fspecial('gaussian', 11, 4), 'replicate','same');
            profileSTAve_array_smooth = imfilter(profileSTAve_array, fspecial('gaussian', 11, 4), 'replicate','same');
            
            STOverInt_array = profileSTAve_array_smooth./profileIntAve_array_smooth;            
            STOverInt_array_smooth = imfilter(STOverInt_array, fspecial('gaussian', 11, 4), 'replicate','same');
            
            STOverInt_array_grad = STOverInt_array_smooth(2:end)-STOverInt_array_smooth(1:end-1);
            
            % this Dpc should not be really close to the nucleus
            STOverInt_array_grad(1:T_center) = nan;
            
            % this Dpc should not be really far away to the nucleus            
            if(numel(STOverInt_array_grad)>=200)
                STOverInt_array_grad(200:end) = nan;
            end
            
            MAX_grad = max(STOverInt_array_grad);
            
            if(MAX_grad<=0 || isnan(MAX_grad))
                ind_d = nan;
            else                
                ind_d = find(STOverInt_array_grad == MAX_grad);
            end
            
            if(isempty(ind_d))
                ind_d = nan;
            end
            
            
            ind_d = ind_d(1);
            
            DistancePeriCenter(iR, iAngle) = ind_d;
           
            if(~isnan(ind_d))
                
                profileIntSumPerpCenterRatio(iR,iAngle) = ...
                    sum(profileIntSum(iAngle,ind_d:end),2)./sum(profileIntSum(iAngle, 1:ind_d-1),2);                
                
                profileIntMeanPerpCenterRatio(iR,iAngle)  = ...
                    (sum(profileIntSum(iAngle,ind_d:end),2)./sum(profilePixels(iAngle,ind_d:end),2))...
                    ./(sum(profileIntSum(iAngle, 1:ind_d-1),2)./sum(profilePixels(iAngle, 1:ind_d-1),2));
                
                profileFilaSumPerpCenterRatio(iR,iAngle)  = ...
                    sum(profileFilaSum(iAngle,ind_d:end),2)...
                    ./sum(profileFilaSum(iAngle, 1:ind_d-1),2);
                
                profileFilaMeanPerpCenterRatio(iR,iAngle)  = ...
                    (sum(profileFilaSum(iAngle,ind_d:end),2)./sum(profilePixels(iAngle,ind_d:end),2))...
                    ./(sum(profileFilaSum(iAngle, 1:ind_d-1),2)./sum(profilePixels(iAngle, 1:ind_d-1),2));
                
                profileStnmsSumPerpCenterRatio(iR,iAngle)  = ...
                    sum(profileStnmsSum(iAngle,ind_d:end),2)...
                    ./sum(profileStnmsSum(iAngle, 1:ind_d-1),2);
                
                profileStnmsMeanPerpCenterRatio(iR,iAngle)  = ...
                    (sum(profileStnmsSum(iAngle,ind_d:end),2)./sum(profilePixels(iAngle,ind_d:end),2))...
                    ./(sum(profileStnmsSum(iAngle, 1:ind_d-1),2)./sum(profilePixels(iAngle, 1:ind_d-1),2));
                
            end
            
        end
    end
        
    
    %25
    output_feature_reorganized.nmssum_ratio_pool_autoDpc   = profileStnmsSumPerpCenterRatio(:);
    %26
    output_feature_reorganized.nmsmean_ratio_pool_autoDpc  = profileStnmsMeanPerpCenterRatio(:);
    %27
    output_feature_reorganized.intsum_ratio_pool_autoDpc   = profileIntSumPerpCenterRatio(:);
    %28
    output_feature_reorganized.intmean_ratio_pool_autoDpc  = profileIntMeanPerpCenterRatio(:);
    %29
    output_feature_reorganized.filasum_ratio_pool_autoDpc  = profileFilaSumPerpCenterRatio(:);
    %30
    output_feature_reorganized.filamean_ratio_pool_autoDpc = profileFilaMeanPerpCenterRatio(:);           

    
    
end



