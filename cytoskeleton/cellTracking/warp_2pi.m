function angle_array = warp_2pi(AA)

angle_array = AA;

for i = 1 : length(AA)-1
    if(angle_array(i+1)>angle_array(i)+2*pi*0.8)
        angle_array(i+1) = angle_array(i+1) -2*pi;
    end
    if(angle_array(i+1)<angle_array(i)-2*pi*0.8)
        angle_array(i+1) = angle_array(i+1) +2*pi;
    end
    if(angle_array(i+1)>angle_array(i)+2*pi*0.8)
        angle_array(i+1) = angle_array(i+1) -2*pi;
    end
    if(angle_array(i+1)<angle_array(i)-2*pi*0.8)
        angle_array(i+1) = angle_array(i+1) +2*pi;
    end
end