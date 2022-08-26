function [res_d,res_t]=averaging_conic(inlier_number,type)

res_d=zeros(10,11,100);
res_t=zeros(10,11,100);

i=1; 
cnt_1=0; %counter

for outlier_ratio=0:0.1:0.90 %first dimention, rows
ii=1;
    for noise_level=0:0.1:1  %second dimension, columns
        for iii=1:1:100
 
        cnt_1=cnt_1+1;
        disp(cnt_1)
        %type=randi([1 3],1);
        [XY_gt, XY, ~, CC]=data_generation_aver(inlier_number,1,1,type,noise_level,outlier_ratio,0);
        %CC_normalized=CC/CC(6);

 
        [C_it_est] = nonlin_fitting_conics(XY,0,0); %no diagnostic %no plotting

%C_it_est=L1_DLT_conic(XY,0,0); %no plotting

[A,B,C,D,E,F]=deal(C_it_est(1),C_it_est(2),C_it_est(3),C_it_est(4),C_it_est(5),C_it_est(6));
A_33=A*C-(B^2)/4;
if((A_33>-1.0000e-15)&&(A_33<1.0000e-15))
%parabola
est_type=3;
elseif(A_33<=-1.0000e-15)
%hyperbola
est_type=2;
elseif(A_33>=1.0000e-15)
%ellypse
est_type=1;
end

%different measures:
% 1) different type
% 2) thresholding parameters
% 3) different type + tresholding parameters

if(type~=est_type)
res_t(i,ii,iii)=1;
end

    
%C_it_est_normalized=C_it_est/C_it_est(6);
C_M=[A B/2 D/2; B/2 C E/2; D/2 E/2 F];
res_d(i,ii,iii)=sampson_distance_conics(C_M,XY_gt,'sum','yes')/size(XY_gt,2);


        end %end of averaging for
        
    ii=ii+1;
    end
i=i+1;
end
end


%if (norm(CC_normalized-C_it_est_normalized)>0.1)
%res_m(i,ii,iii)=1;
%end