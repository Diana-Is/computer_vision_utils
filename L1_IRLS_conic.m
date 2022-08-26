function s=L1_IRLS_conic(XY,if_plot,if_diag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: The coordinates x,y, mu's initial value, eps -error tolerance
%Output: The estimated model s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1: Mapping correspondences into embeddings a i or b i to form the data matrix M.
XYs=[XY;ones(1,size(XY,2))];
%newXY=XY;
[newXY, T] = precondition(XYs);
xw=newXY(1,:)';
yw=newXY(2,:)';
a=my_IRLS(xw,yw,100,if_diag);
[A,B,C,D,E,F]=deal(a(2),a(3),a(1),a(5),a(4),1);
C_M_est=[A,B/2,D/2; B/2,C,E/2; D/2,E/2,F];

% denormalise C
C_M_est = T'*C_M_est*T;
s=[C_M_est(1,1),C_M_est(1,2)*2,C_M_est(2,2),C_M_est(1,3)*2,C_M_est(2,3)*2,C_M_est(3,3)];

if(if_plot==1)
plot_a_conic(C_M_est,XY,'L1 IRLS conic estimation');
end


end 