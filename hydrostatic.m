function [VCB,LCB,LCF,Disp,I,AWP] = hydrostatic(megaarray,delx,x,t_i,z_i)
ns = size(megaarray,3);
sum_zy = zeros(2,1);
sum_v = zeros(3,1);
sum_a = zeros(2,1);
flag = 0;
sum_zy_prev = zeros(2,1);

    for i = 1:ns     %ns
        if t_i(i,1)>0
            y_t = megaarray(t_i(i,1),2,i)+(megaarray(t_i(i,1)+1,2,i)-megaarray(t_i(i,1),2,i))/(megaarray(t_i(i,1)+1,1,i)-megaarray(t_i(i,1),1,i))*(z_i(i)-megaarray(t_i(i,1),1,i));
            for j = 1:t_i(i,1)-1
                    sum_zy(2) = sum_zy(2) + (0.5)^2*(megaarray(j,1,i)+megaarray(j+1,1,i))*(-megaarray(j,1,i)+megaarray(j+1,1,i))*(megaarray(j+1,2,i)+megaarray(j,2,i));
                    sum_zy(1) = sum_zy(1) + 0.5*(-megaarray(j,1,i)+megaarray(j+1,1,i))*(megaarray(j+1,2,i)+megaarray(j,2,i));
                    % if(i==23)
                    %     disp(sum_zy);
                    % end
            end
        end
            if (t_i(i,2)==0) %%error in condition
                sum_zy(2) = sum_zy(2) + (0.5)^2*(megaarray(t_i(i,1),1,i)+z_i(i))*(-megaarray(t_i(i,1),1,i)+z_i(i))*(y_t+megaarray(t_i(i,1),2,i));%multiply z_i (for the linear line)
                sum_zy(1) = sum_zy(1) + 0.5*(-megaarray(t_i(i,1),1,i)+z_i(i))*(y_t+megaarray(t_i(i,1),2,i));
            end
            if (t_i(i,1)>0)
                if flag==1
                    sum_v(3) = sum_v(3) + delx(i-1)*(sum_zy(1)+sum_zy_prev(1))*(0.5)^2*(x(i)+x(i-1));
                    sum_v(2) = sum_v(2) + delx(i-1)*(sum_zy(2)+sum_zy_prev(2))*0.5;
                    sum_v(1) = sum_v(1) + delx(i-1)*(sum_zy(1)+sum_zy_prev(1))*0.5;
                end
                sum_zy_prev(1) = sum_zy(1);
                sum_zy_prev(2) = sum_zy(2);
                y_t_prev = y_t;
                sum_zy(2) = 0;
                sum_zy(1) = 0;
            end
                    % if(sum_zy(2) == 0)
            % disp(t_i(i,1));
            % disp(y_t);
                    % end
            % if(sum_zy(1)==0)
            %     break;
            % end
            if flag==1            
            % disp(sum_zy(1));
            % disp(i);
                    if (t_i(i,2)==0)&&(t_i(i,1)>0)
                        % if (isnan(sum_a(2)/sum_a(1)))
                        %     disp('x=')
                        %     disp(x(i))
                        %     disp('y=')
                        %     disp(y_t)
                        %     disp('z=')
                        %     disp(z_i(i))
                        % end
                        sum_a(1) = sum_a(1) + 0.5*(y_t+y_t_prev)*sqrt(delx(i-1)^2+(z_i(i)-z_i(i-1))^2);
                        sum_a(2) = sum_a(2) + (0.5)^2*(x(i)+x(i+1))*(y_t+y_t_prev)*delx(i-1);
                    end
            end
            flag = 1;
            %disp(sum_zy)
            % if (sum_zy(2)==-inf)
            %     disp(i);
            % end
    end
    sum_a(:) = 2*sum_a(:);
    Disp = sum_v(1)*2;
    VCB = sum_v(2)*2/Disp;
    LCB = sum_v(3)*2/Disp;
    LCF = sum_a(2)/sum_a(1);
    AWP = sum_a(1,1);
    I = 0;
    if isnan(LCF)
        I=0;
    else
            flag = 0;
        for i=1:ns
            if t_i(i,1)>0
                y_t = megaarray(t_i(i,1),2,i)+(megaarray(t_i(i,1)+1,2,i)-megaarray(t_i(i,1),2,i))/(megaarray(t_i(i,1)+1,1,i)-megaarray(t_i(i,1),1,i))*(z_i(i)-megaarray(t_i(i,1),1,i));
                if flag == 1
                    I = I + (0.5)*((x(i)+x(i+1))/2-LCF)^2*(y_t+y_t_prev)*delx(i-1);
                end
                y_t_prev = y_t;
                flag = 1;
            end
        end
    end
end
    
    %disp(t_i(:,2));

