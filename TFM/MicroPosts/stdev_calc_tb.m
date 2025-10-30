% CAL 1/30/03 
% Revised CAL 2/5/03
%
% calculate standard deviation of edge posts

% calculate number of unoccupied posts

border_num=length(centr2_bot)-cell_num-length(find(degen_posts_x));
if border_num ~= 0
        sum_mag=0;
        total_border_x=0;
        total_border_y=0;

        % Add deflections of unoccupied posts as well as square of deflections

        for i=1:row;
           for j=1:col;
              sum_mag=sum_mag+border_posts_mag(i,j)^2;
              total_border_x=total_border_x+border_posts_x(i,j);
              total_border_y=total_border_y+border_posts_y(i,j);

           end
        end

        for i=1:row;
            for j=1:col;
                if and(degen_posts_x(i,j) ~= 0,degen_posts_y(i,j) ~= 0);
                    actual_x_bot_micron(i,j)=(sum(actual_x_bot_micron(i,:))-actual_x_bot_micron(i,j))/(col-1);
                    actual_y_bot_micron(i,j)=(sum(actual_y_bot_micron(:,j))-actual_y_bot_micron(i,j))/(row-1);
                end
            end
        end

        for i=1:row;
            for j=1:col;
                    std_bot_x(i)=std(actual_x_bot_micron(i,:));
                    std_bot_y(j)=std(actual_y_bot_micron(:,j));
            end
        end

        std_bot=(mean(std_bot_x)+mean(std_bot_y))/2;

        % Calculate x and y averages of error, as well as std dev

        avg_grid_err_x=total_border_x/border_num;
        avg_grid_err_y=total_border_y/border_num;
        stdev=sqrt(sum_mag/border_num);


else
     avg_grid_err_x=0;
     avg_grid_err_y=0;
     stdev=0;
     std_bot=0;
end

        
