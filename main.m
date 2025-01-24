clear all;
clc;

load('sections.mat'); %code restricted to 100 datapoints per section
n = size(sections,1);
megaarray = zeros(100,2,106);
i=2;
j=0;
k = 1;
ns = sections(1,1); %no. of sections
z_i = input("z_i = ");
t_i = zeros(ns,1);
delx = zeros(ns-1,1);
%disp(size(delx));
flag=0;
xarray = zeros(ns,1);
while(i+j<=n)
    %disp(i);   
    array = zeros(100,2);
    j = sections(i,2);
    
    if(k~=ns)
    x1 = sections(i,1);
    x2 = sections(i+j+1,1);
    end
    xarray(k) = sections(i,1);
    %disp(j);
    %disp(x2);
            if k == 107
            disp("Hell naw");
        end
    delx(k) = x2-x1;
    for t=1:j-1

        array(t,1)=sections(t+i,1);
        array(t,2)=sections(t+i,2);
        array(t+1,1)=sections(t+i+1,1);
        array(t+1,2)=sections(t+i+1,2);
        if ((array(t+1,1)>=z_i)&&(array(t,1)<=z_i))
            t_i(k) = t;
            flag = 1;
            % if((k==23))
            %     disp(array(t+1,2));
            %     disp(array(t+1,2));
            %     disp(z_i);
            %     disp(t);
            %     disp(j);
            %     disp(i);
            % end
            break;
        elseif (((array(t+1,1)<z_i)&&(array(t,1)<z_i)))
            t_i(k) = t;                    
            flag = 1;
        end
    end
    if (flag == 0)
            t_i(k) = 0;
    end
    flag = 0;
    megaarray(:,:,k) = array;
    if(k==ns-1)  
        break;
    end
    k = k+1;

    i = i+j+1;
    %disp(size(delx));
end


%disp(size(delx));

[VCB,LCB,LCF] = vcb(megaarray,delx,xarray,t_i,z_i);

fprintf('VCB = %d\n', VCB);
fprintf('LCB = %d\n', LCB);
fprintf('LCF = %d\n', LCF);
% mal1 = ma_l(file,1);
% LCF = mal2/mal1;
% mav2 = ma_v(file,2);
% mav1 = ma_v(file,1);
% VCF = mav2/mav1;
% LCB = 




