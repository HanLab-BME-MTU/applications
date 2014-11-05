function O_filter = OrientationVote(input_Orientation,SteerabelRes_Segment, pace, radius)

DO = input_Orientation*2;
DO_filter = DO;

% radius = 5;
H = fspecial('gaussian',2*radius+1,5);
half_pace  = (pace-1)/2;

S1 = size(DO,1);
S2 = size(DO,2);

for i = 1 :pace: size(DO,1)
    for j = 1 :pace: size(DO,2)
        local_seg = SteerabelRes_Segment(max(1,i-radius):min(i+radius,S1), max(1,j-radius):min(j+radius,S2));
        local_DO_patch = DO(max(1,i-radius):min(i+radius,S1), max(1,j-radius):min(j+radius,S2));
        if(~isempty(find(local_seg>0)))
            loacal_DO_list = local_DO_patch(local_seg  >0);
            [h_bin,i_count] = rose(loacal_DO_list,8);
            h_bin_center = (h_bin(2:4:end) +h_bin(3:4:end))/2;
            i_count = i_count(2:4:end);
            [inda,indb] = find(i_count==max(i_count));
            gross_peak = h_bin_center(indb(1));
            
            loacal_DO_list_center0 = loacal_DO_list-gross_peak;
            [h_bin,i_count] = rose(loacal_DO_list_center0,16);
            h_bin_center = (h_bin(2:4:end) +h_bin(3:4:end))/2;
            i_count = i_count(2:4:end);            
            [inda,indb] = find(i_count==max(i_count));
            DO_filter(max(1,i-half_pace):min(S1,i+half_pace), max(1,j-half_pace):min(S2,j+half_pace)) = h_bin_center(indb(1)) - pi/16 + gross_peak;            
        else
            DO_filter(max(1,i-half_pace):min(S1,i+half_pace), max(1,j-half_pace):min(S2,j+half_pace)) = NaN;
        end
    end
end

DO_filter(find(DO_filter>pi)) = DO_filter(find(DO_filter>pi)) - 2*pi;
DO_filter(find(DO_filter<-pi)) = DO_filter(find(DO_filter<-pi)) + 2*pi;
DO_filter(find(DO_filter>pi)) = DO_filter(find(DO_filter>pi)) - 2*pi;
DO_filter(find(DO_filter<-pi)) = DO_filter(find(DO_filter<-pi)) + 2*pi;

O_filter = DO_filter/2;