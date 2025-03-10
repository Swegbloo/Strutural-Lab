function [VCB,LCB,LCF,Disp,I,xarray,AWP] = hydrostat_properties(z_i)


load('sections.mat'); %code restricted to 100 datapoints per section

n = size(sections,1);
ns = sections(1,1); %no. of sections
% z_i = [12.8774
%    12.8773
%    12.8118
%    12.6808
%    12.6153
%    12.4844
%    12.4189
%    12.3263
%    12.1427
%    11.9590
%    11.9468
%    11.9461
%    11.9131
%    11.8728
%    11.8327
%    11.7523
%    11.5917
%    11.0874
%    11.0873
%    11.0863
%    10.9130
%    10.9123
%    10.2142
%     9.5714
%     8.9286
%     8.2858
%     7.6429
%     7.0001
%     6.3573
%     5.7145
%     5.0717
%     4.4288
%     3.7860
%     3.1432
%     2.5003
%     1.8575
%     1.2147
%     0.5719
%    -0.0709
%    -0.7138
%    -1.3566
%    -1.9994
%    -2.6422
%    -3.1932
%    -3.3190
%    -3.3191
%    -3.5605
%    -3.7442
%    -3.8429
%    -3.8430
%    -4.2952
%    -4.5248
%    -4.7084
%    -4.8462
%    -4.9380
%    -5.0069
%    -5.0413
%    -5.0586
%    -5.0672
%    -5.0693
%    -5.0715
%    -5.0758
%    -5.0786
%    -5.0815
%    -5.0873
%    -5.0987
%    -5.1217
%    -5.1676
%    -5.2135
%    -5.2365
%    -5.2594
%    -5.2824
%    -5.3054
%    -5.3226
%    -5.3340
%    -5.3426
%    -5.3470
%    -5.3484
%    -5.3491
%    -5.3513
%    -5.3570
%    -5.3627
%    -5.3656
%    -5.3685
%    -5.3713
%    -5.3742
%    -5.3771
%    -5.3800
%    -5.3814
%    -5.3828
%    -5.3843
%    -5.3857
%    -5.3871
%    -5.3885
%    -5.3893
%    -5.3900
%    -5.3907
%    -5.3914
%    -5.3921
%    -5.3929
%    -5.3936
%    -5.3943
%    -5.3950
%    -5.3957
%    -5.3965
%    -5.3972];
%z_i=2.62662599769878*ones(ns,1);
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

    if z_i(k)>sections(i+j,2) || z_i(k)<sections(i+1,2)
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
[VCB,LCB,LCF,Disp,I,AWP] = hydrostatic(megaarray,delx,xarray,t_i,z_i);

% fprintf('Displacement = %d\n', Disp);
% fprintf('VCB = %d\n', VCB);
% fprintf('LCB = %d\n', LCB);
% fprintf('LCF = %d\n', LCF);
% fprintf('I = %d\n', I);





