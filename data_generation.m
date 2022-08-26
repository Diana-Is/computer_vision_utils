function [XY, C_Mat, CC]=data_generation(cntr,rot,type,noise,outliers,if_plot)
[x_n, y_n, C_v]=conic_5_points_generation(cntr,rot,type);

[A,B,C,D,E,F]=deal(C_v(1),C_v(2),C_v(3),C_v(4),C_v(5),C_v(6));
CC=[A,B,C,D,E,F];
CC=CC/norm(CC);
C_Mat=[A B/2 D/2; B/2 C E/2; D/2 E/2 F];
%C_Mat=C_Mat/norm(C_Mat);

if(noise>0)
nn=noise*randn(2,5);
x_n=x_n+nn(1,:);
y_n=y_n+nn(2,:);
end

if(outliers~=0)
outl=10*rand(2,10);
x_n=[x_n,outl(1,:)];
y_n=[y_n,outl(2,:)];
end

XY=[x_n;y_n];

if(if_plot==1)
plot_a_conic(C_Mat,XY,'ground truth, generated');
end

end



function [x, y, C_v]=conic_5_points_generation(cntr,rot,type)
%Q(x,y)=Ax^2+Bxy+Cy^2+Dx+Ey+F=0


t = unifrnd(-pi,pi,1,5);

if(type==1)%ellypse

params = unifrnd(1,10,1,4);
a=params(1);
b=params(2);
xc=params(3);
yc=params(4);

if(cntr==1)%if centered    
xc=0;
yc=0;
end

theta=0;
if(rot==1)
theta=unifrnd(0,pi/2);
end

ct=cos(theta);
st=sin(theta);

x=a*cos(t);
y=b*sin(t);

R=[ct,-st;st,ct];
tr=[xc;yc];
crd=[x;y];

crd=R*crd+tr;
x=crd(1,:);
y=crd(2,:);

C_v=abcd_to_ABCDEF(a,b,xc,yc,theta,type);



elseif(type==2)%hyperbola
  
params = unifrnd(1,10,1,4);
a=params(1);
b=params(2);
xc=params(3);
yc=params(4);

if(cntr==1)%if centered    
xc=0;
yc=0;
end

   
theta=0;
if(rot==1)
theta=unifrnd(0,pi/2);
end

x=a*sec(t);
y=b*tan(t);

ct=cos(theta);
st=sin(theta);

R=[ct,-st;st,ct];
tr=[xc;yc];
crd=[x;y];

crd=R*crd+tr;
x=crd(1,:);
y=crd(2,:);
    
C_v=abcd_to_ABCDEF(a,b,xc,yc,theta,type);


elseif(type==3)%parabola
 
params = unifrnd(1,10,1,3);
a=params(1);
xc=params(2);
yc=params(3);

if(cntr==1)%if centered    
xc=0;
yc=0;
end    
    

theta=0;
if(rot==1)
theta=unifrnd(0,pi/2);
end

x=0.5*a*t.^2;
y=a*t;


ct=cos(theta);
st=sin(theta);

R=[ct,-st;st,ct];
tr=[xc;yc];
crd=[x;y];

crd=R*crd+tr;
x=crd(1,:);
y=crd(2,:);
    
C_v=abcd_to_ABCDEF(a,0,xc,yc,theta,type);
    

end



end