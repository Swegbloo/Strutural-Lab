clear all;
clc;

load('sections.mat');
n = size(sections,1);

i=2;
j=0;
while((i<=n)&&(i+j+1<=n))
    %disp(i);   
    array = zeros(50,2);
    disp(i+j+1);
    j = sections(i,2);
    %disp(j);
    x1 = sections(i,1);
    x2 = sections(i+j+1,1);
    for t=2:j+1
        array(1,1)=x2-x1;
        array(1,2)=0;
        array(t,1)=sections(t+i-1,1);
        array(t,2)=sections(t+i-1,2);
    end
    i = i+j+1;
end

% mal2 = ma_l(file,2); % write compiled file with degree
% mal1 = ma_l(file,1);
% LCF = mal2/mal1;
% mav2 = ma_v(file,2);
% mav1 = ma_v(file,1);
% VCF = mav2/mav1;
% LCB = 




