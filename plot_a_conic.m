function plot_a_conic(C_est,XY,my_title,XYc)

A=C_est(1,1); B=2*C_est(1,2); C=C_est(2,2); D=2*C_est(1,3); E=2*C_est(2,3); F=C_est(3,3);

x=XY(1,:);
y=XY(2,:);

lim_X_min=min(XY(1,:));
lim_X_max=max(XY(2,:));
lim_Y_min=min(XY(1,:));
lim_Y_max=max(XY(2,:));

if (nargin<4||isempty(XYc))
figure
fimplicit(@(s,t) A*s.^2+s.*t.*B+C*t.^2+s.*D+t.*E+F, [lim_X_min-20,lim_X_max+20,lim_X_min-20,lim_X_max+20])
hold on
plot(x,y,'*')
title(my_title)
else
xc=XYc(1,:);
yc=XYc(2,:);

figure
fimplicit(@(s,t) A*s.^2+s.*t.*B+C*t.^2+s.*D+t.*E+F, [lim_X_min-20,lim_X_max+20,lim_X_min-20,lim_X_max+20])
hold on
plot(x,y,'*')
hold on
plot(xc,yc,'*','color','k')
title(my_title)

end
end