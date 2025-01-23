function [VCB] = vcb(megaarray,delx,x,t_i,z_i)
ns = size(megaarray,3);
sum_zy = zeros(2,1);
sum_xyz = zeros(2,1);
    for i = 1:ns
        if t_i(i,1)>0
            y_t = megaarray(t_i(i,1),2,i)+(megaarray(t_i(i,1)+1,2,i)-megaarray(t_i(i,1),2,i))/(megaarray(t_i(i,1)+1,1,i)-megaarray(t_i(i,1),1,i))*(z_i-megaarray(t_i(i,1),1,i));
            for j = 1:t_i(i,1)-1
                    sum_zy(2) = sum_zy(2) + 0.5*(megaarray(j,1,i)+megaarray(j+1,1,i))^2*(megaarray(j+1,2,i)-megaarray(j,2,i));
                    sum_zy(1) = sum_zy(1) + (megaarray(j,1,i)+megaarray(j+1,1,i))*(megaarray(j+1,2,i)-megaarray(j,2,i));
                    if(i==23)
                        %disp(sum_zy);
                    end
            end
            sum_zy(2) = sum_zy(2) + 0.5*(megaarray(t_i(i,1),1,i)+z_i)^2*(y_t-megaarray(t_i(i,1),2,i));%multiply z_i (for the linear line)
            sum_zy(1) = sum_zy(1) + (megaarray(t_i(i,1),1,i)+z_i)*(y_t-megaarray(t_i(i,1),2,i));
                    % if(sum_zy(2) == 0)
                    %     disp(i);
                    % end
            if(sum_zy(1)==0)
                break;
            end
            sum_xyz(2) = sum_xyz(2) + delx(i)*(sum_zy(2)/sum_zy(1))*x(i);
            sum_xyz(1) = sum_xyz(1) + delx(i)*(sum_zy(2)/sum_zy(1));
            %disp(sum_zy)
            % if (sum_zy(2)==-inf)
            %     disp(i);
            % end
            sum_zy(2) = 0;
            sum_zy(1) = 0;
        end
    end
    VCB = sum_xyz(2)/sum_xyz(1);
end
