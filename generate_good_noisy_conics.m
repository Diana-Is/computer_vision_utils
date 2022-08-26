function [M,Mgt,CM]=generate_good_noisy_conics(num_of_inl,type)
M=zeros(22,num_of_inl,100);%2 rows*11rows = 22 rows, n columns because n inliers and 100 samples
Mgt=zeros(22,num_of_inl,100);%2 rows*11rows = 22 rows, n columns because n inliers and 100 samples
CM=zeros(11,6,100); % 11 rows 6 columns, 100 samples
for j=1:100
i=1;
k=1;
for noise_level=0:0.1:1
[XY_gt, XY, ~, CC]=data_generation_aver(num_of_inl,0,1,type,noise_level,0,0);
M(i:i+1,:,j)=XY;
Mgt(i:i+1,:,j)=XY_gt;
CM(k,:,j)=CC;
i=i+2;
k=k+1;
end
end
end