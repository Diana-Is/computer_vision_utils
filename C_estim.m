function [C_v_est,C_M_est]=C_estim(XY,scale,if_plot)
XYs=XY;
if(scale==1)
if(size(XY,1)==2)
XY=[XY;ones(1,size(XY,2))];
end
[newXY, T] = precond(XY);
xw=newXY(1,:)';
yw=newXY(2,:)';
else 
xw=XY(1,:);
yw=XY(2,:);
T=eye(3);
end

A=[xw.^2 xw.*yw yw.^2 xw yw ones(size(xw))];
[~,~,V] = svd(A);
a=V(:,end);

[A,B,C,D,E,F]=deal(a(1),a(2),a(3),a(4),a(5),a(6));
C_M_est=[A,B/2,D/2; B/2,C,E/2; D/2,E/2,F];

%a=my_IRLS(xw,yw);
%[A,B,C,D,E,F]=deal(a(2),a(3),a(1),a(5),a(4),1);
%C_M_est=[A,B/2,D/2; B/2,C,E/2; D/2,E/2,F];

% denormalise C
C_M_est = T'*C_M_est*T;
C_v_est=[C_M_est(1,1),C_M_est(1,2)*2,C_M_est(2,2),C_M_est(1,3)*2,C_M_est(2,3)*2,C_M_est(3,3)];

if(if_plot==1)
plot_a_conic(C_M_est,XYs,'least square solution (just SVD)');
end


end




