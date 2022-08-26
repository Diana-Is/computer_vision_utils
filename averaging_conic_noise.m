function [res_d,res_t,time_m]=averaging_conic_noise(type,alg_type)

% % if (type==1)
% % load ./ellypse_insiem2.mat; 
% % elseif (type==2)
% % load ./hyperbola_insiem2.mat;
% % end


res_d=zeros(4,11,100); %4 -- levels of inliers,5,10,15,20
res_t=zeros(4,11,100); %then 11 levels of noise from 0.1 to 1, and averaged over 100 samples
time_m = zeros(4,11,100);

i=1;
k=1;
cnt=0;
for inlier_number=5:5:20 %first dimention, rows
ii=1;
k=1;
    for noise_level=0:0.1:1  %second dimension, columns
        for iii=1:1:100
        cnt=cnt+1;
        disp(cnt);
        %XY=LL(k:k+1,1:inlier_number,iii);
        
        [XY_gt, XY, C_Mat, CC]=data_generation_aver(inlier_number,0,1,type,noise_level,0,0);
        tic;
        if strcmp(alg_type,'LM')
        [C_it_est] = nonlin_fitting_conics(XY,0,0); %no diagnostic %no plotting
        elseif strcmp(alg_type,'subgrad')
        C_it_est=L1_DLT_conic(XY,0,0); %no plotting
        elseif strcmp(alg_type,'IRLS')
        C_it_est=L1_IRLS_conic(XY,0,0);
        elseif strcmp(alg_type,'ISTA')
            
        end
        time_m(i,ii,iii)=toc;
        
       
[A,B,C,D,E,F]=deal(C_it_est(1),C_it_est(2),C_it_est(3),C_it_est(4),C_it_est(5),C_it_est(6));
A_33=A*C-(B^2)/4;
if((A_33>-1.0000e-10)&&(A_33<1.0000e-10))
%parabola
est_type=3;
elseif(A_33<=-1.0000e-10)
%hyperbola
est_type=2;
elseif(A_33>=1.0000e-10)
%ellypse
est_type=1;
end

%different measures:
% 1) different type
% 2) thresholding parameters

if(type~=est_type)
res_t(i,ii,iii)=1;
end
    
%XY_gt=GG(k:k+1,1:inlier_number,iii);
C_M=[A B/2 D/2; B/2 C E/2; D/2 E/2 F];
res_d(i,ii,iii)=sampson_distance_conics(C_M,XY_gt,'sum','yes')/size(XY_gt,2);

            end %end of averaging 100 for        
ii=ii+1;
k=k+2;
    end%end of noise for
i=i+1;
end %end of inlier for

end %end of function