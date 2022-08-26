function [XY_gt, XY, C_Mat, CC]=data_generation_aver(n_inliers,cntr,rot,type,noise,outliers_ratio,if_plot)
%function generating points(with/without noise or outliers) of a conic
%INPUT:
%-n_inliers -- number of inliers, points, lying on conic, possibly noisy
%-cntr --- if the conic should be necessarily centered(1/0-centered/not)
%-rot --- if the conic can be rotated (1/0 - rotated/not)
%-type --- (1-ellypse,2-hyperbola,3-parabola)
%-noise --- power (in [0,1]) of gaussian noise which is added to inlier points
%-outliers ratio --- (in [0,1]) ratio of added outliers wrt the inliers
%-if plot --- if plotting the result of generation, conic and point on it
%OUTPUT:
%XY_gt --- ground truth points lying on conic,no noise no outliers
%XY --- result of generated conic points, with possibly noise and outliers
%C_Mat --- matrix of algebraic conic parameters
%CC --- vector of algebraic conic parameters

[x_n, y_n, C_v, x_cn, y_cn]=conic_points_generation(n_inliers,cntr,rot,type);
outliers=[];


if(noise>0||outliers_ratio~=0)
[XY,outliers]=corrupt(x_n,y_n,noise,outliers_ratio,n_inliers,x_cn, y_cn);

%generated points need to be checked for collinearity
k=1;
flag=0;
flag_glag=0;
thrsh=0.75;
while k<=length(x_n)-2 %first index
      
for i=k+1:length(x_n)-1
    for j=i+1:length(x_n)
    test_set=[XY(1,k),XY(2,k);XY(1,i),XY(2,i);XY(1,j),XY(2,j)];
    %dts12=pdist([XY(1,k),XY(2,k);XY(1,i),XY(2,i)],'euclidean');
    %dts13=pdist([XY(1,k),XY(2,k);XY(1,j),XY(2,j)],'euclidean');
    %dts23=pdist([XY(1,i),XY(2,i);XY(1,j),XY(2,j)],'euclidean');
    distances=pdist(XY','euclidean');
    if(min(distances)<=thrsh)
    flag_glag=1;
    end
        if  ((if_collinear(test_set)==1)||flag_glag==1)
            [x_n, y_n, C_v, x_cn, y_cn]=conic_points_generation(n_inliers,cntr,rot,type);
            [XY,outliers]=corrupt(x_n,y_n,noise,outliers_ratio,n_inliers,x_cn, y_cn);
            flag=1;
            break;
        end
    end
    if(flag==1)
    flag=0;
    flag_glag=0;
    k=0;
    break;
    end
end
k=k+1;
end %end of testing collinearity while

elseif(noise==0&&outliers_ratio==0) %if no noise, no outliers
    
XY=[x_n;y_n];
k=1;
flag=0;
while k<=length(x_n)-2 %first index
for i=k+1:length(x_n)-1
    for j=i+1:length(x_n)
    test_set=[XY(1,k),XY(2,k);XY(1,i),XY(2,i);XY(1,j),XY(2,j)];
    dts12=pdist([XY(1,k),XY(2,k);XY(1,i),XY(2,i)],'euclidean');
    dts13=pdist([XY(1,k),XY(2,k);XY(1,j),XY(2,j)],'euclidean');
    dts23=pdist([XY(1,i),XY(2,i);XY(1,j),XY(2,j)],'euclidean');
    thrsh=0.75;
        if  ((if_collinear(test_set)==1)||dts12<=thrsh||dts13<=thrsh||dts23<=thrsh)
            [x_n, y_n, C_v,~,~]=conic_points_generation(n_inliers,cntr,rot,type);
            XY=[x_n;y_n];
            flag=1;
            break;
        end
    end
    if(flag==1)
    flag=0;
    k=0;
    break;
    end
end
k=k+1;
end %end of testing collinearity while

end

XY_gt=[x_n;y_n];
[A,B,C,D,E,F]=deal(C_v(1),C_v(2),C_v(3),C_v(4),C_v(5),C_v(6));
CC=[A,B,C,D,E,F];
CC=CC/norm(CC);
C_Mat=[A B/2 D/2; B/2 C E/2; D/2 E/2 F];

if(if_plot==1)
plot_a_conic(C_Mat,XY,'ground truth, generated',outliers);
end

end %end of main function

function [XYn,outl]=corrupt(xn,yn,noise,outliers_ratio,n_inliers,x_c,y_c)
range_x=abs(max(xn)-min(xn));
range_y=abs(max(yn)-min(yn));

nn=[0.05*range_x;0.05*range_y];

if(noise>0)
nn=noise*nn.*randn(2,n_inliers);
xn=xn+nn(1,:);
yn=yn+nn(2,:);
end

if(outliers_ratio~=0)
outl_temp=unifrnd(-max(range_x,range_y)/2,max(range_x,range_y)/2,2,round(outliers_ratio*n_inliers));
outl=outl_temp+[x_c;y_c];
xn=[xn,outl(1,:)];
yn=[yn,outl(2,:)];
else 
outl=[];
end

XYn=[xn;yn];
end




function tf = if_collinear(pts)
%function checking collinearity of each 3 points in the array
%of generated points
mat = [pts(1,1)-pts(3,1),pts(1,2)-pts(3,2); ...
      pts(2,1)-pts(3,1),pts(2,2)-pts(3,2)];
%  det(mat)
tf=(det(mat)>-0.05)&&(det(mat)<0.05);
end


function [x, y, C_v, xc, yc]=conic_points_generation(n_inliers,cntr,rot,type)
%Q(x,y)=Ax^2+Bxy+Cy^2+Dx+Ey+F=0

t = unifrnd(-pi,pi,1,n_inliers);

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

% 1     function tf = collinear7(p1,p2,p3)
% 2     mat = [p1(1)-p3(1) p1(2)-p3(2); ...
% 3            p2(1)-p3(1) p2(2)-p3(2)];
% 4     tf = det(mat) == 0;