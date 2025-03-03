function [VCB,LCB,LCF,Disp,I] = hydrostat_properties(z_i)

load('sections.mat'); %code restricted to 100 datapoints per section

n = size(sections,1);
ns = sections(1,1); %no. of sections

megaarray = zeros(100,2,ns);
i=2;
j=0;
k = 1;


%z_i = input("z_i = ");
t_i = zeros(ns,2); 
delx = zeros(ns-1,1);
flag=0;
xarray = zeros(ns,1);
% z_i = zeros(ns,1);
% z = input('z_i =');
% for j=1:ns
%     z_i(j) = z;
% end

% while(i+j<=n)   
%     array = zeros(100,2);
%     j = sections(i,2);
% 
%     if(k~=ns)
%     x1 = sections(i,1);
%     x2 = sections(i+j+1,1);
%     end
%     if z_i>sections(i+j,1)
%         t_i(k,2) = 1;
%     else
%         t_i(k,2) = 0;
%     end
%     xarray(k) = sections(i,1);
%     % delx(k) = x2-x1;
%     for t=1:j-1
% 
%         array(t,1)=sections(t+i,1);
%         array(t,2)=sections(t+i,2);
%         array(t+1,1)=sections(t+i+1,1);
%         array(t+1,2)=sections(t+i+1,2);
%         if ((array(t+1,1)>=z_i)&&(array(t,1)<=z_i))
%             t_i(k,1) = t;
%             flag = 1;
%             break;
%         elseif (((array(t+1,1)<z_i)&&(array(t,1)<z_i)))
%             t_i(k,1) = t;                    
%             flag = 1;
%         end
%     end
%     if (flag == 0)
%             t_i(k,1) = 0;
%     end
%     flag = 0;
%     megaarray(:,:,k) = array;
%     if(k<ns)  
%         delx(k) = x2-x1;
%     end
%     k = k+1;
% 
%     i = i+j+1;
% end

while(i+j<=n)   
    array = zeros(100,2);
    j = sections(i,2);
    
    if(k~=ns)
    x1 = sections(i,1);
    x2 = sections(i+j+1,1);
    end
    
    %disp(sections(i+j,2));
    if z_i(k)>sections(i+j,2)
        t_i(k,2) = 1;
    else
        t_i(k,2) = 0;
    end
    
    xarray(k) = sections(i,1);
    % delx(k) = x2-x1;
    for t=1:j-1

        array(t,1)=sections(t+i,2);
        array(t,2)=sections(t+i,1);
        array(t+1,1)=sections(t+i+1,2);
        array(t+1,2)=sections(t+i+1,1);
        if ((array(t+1,1)>=z_i(k))&&(array(t,1)<=z_i(k)))
            t_i(k,1) = t;
            flag = 1;
            break;
        elseif (((array(t+1,1)<z_i(k))&&(array(t,1)<z_i(k))))
            t_i(k,1) = t;                    
            flag = 1;
        end
    end
    if (flag == 0)
            t_i(k,1) = 0;
    end
    flag = 0;
    megaarray(:,:,k) = array;
    if(k<ns)  
        delx(k) = x2-x1;
    end
 
    k = k+1;
    if k==107
        break;
    end
    i = i+j+1;
end

[VCB,LCB,LCF,Disp,I] = hydrostatic(megaarray,delx,xarray,t_i,z_i);

fprintf('Displacement = %d\n', Disp);
fprintf('VCB = %d\n', VCB);
fprintf('LCB = %d\n', LCB);
fprintf('LCF = %d\n', LCF);
fprintf('I = %d\n', I);





