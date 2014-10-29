function angle_out = orientaion_list_org(angle)

angle(find(angle>pi/2)) = angle(find(angle>pi/2)) -pi;
angle(find(angle<-pi/2)) = angle(find(angle<-pi/2)) +pi;

angle(find(angle>pi/2)) = angle(find(angle>pi/2)) -pi;
angle(find(angle<-pi/2)) = angle(find(angle<-pi/2)) +pi;

angle(find(angle>pi/2)) = angle(find(angle>pi/2)) -pi;
angle(find(angle<-pi/2)) = angle(find(angle<-pi/2)) +pi;

angle(find(angle>pi/2)) = angle(find(angle>pi/2)) -pi;
angle(find(angle<-pi/2)) = angle(find(angle<-pi/2)) +pi;

angle(find(angle>pi/2)) = angle(find(angle>pi/2)) -pi;
angle(find(angle<-pi/2)) = angle(find(angle<-pi/2)) +pi;

for i = 1 : length(angle)-1
    if angle(i+1) - angle(i) > pi*2/3
        angle(i+1:end) = angle(i+1:end) - pi;
    else
        if angle(i+1) - angle(i) < -pi*2/3
            angle(i+1:end) = angle(i+1:end) + pi;
        end
    end
end


angle_out = angle;