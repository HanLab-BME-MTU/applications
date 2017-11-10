
sort_x = digital_x_section(1);
sort_y = digital_y_section(1);
sort_ang = [];
angle_pool = angles_section(1);

if(length(digital_x_section)>1)
    for j = 2 : length(digital_x_section)
        if digital_x_section(j)==sort_x(end) && digital_y_section(j)==sort_y(end)
            %same point as previous, put the angle to the pool
            angle_pool = [angle_pool; angles_section(j)];
        else
            % new digital point? add to the sorted x and put the angle of
            % last group to the
            sort_x = [sort_x;digital_x_section(j)];
            sort_y = [sort_y;digital_y_section(j)];
            sort_ang = [sort_ang; mean(angle_pool)];
            angle_pool = angles_section(j);
        end
    end
    sort_ang = [sort_ang; mean(angle_pool)];
else
    sort_ang = NaN;
end

scrable_digital_model{iFila_section} = [sort_x sort_y];
scrable_orientation_model{iFila_section} = sort_ang;

scrable_II = [scrable_II;iFila_section*(sort_x*0+1)];
scrable_XX = [scrable_XX;sort_x;];
scrable_YY = [scrable_YY;sort_y;];
scrable_OO = [scrable_OO;sort_ang;];