function C6=abcd_to_ABCDEF(a,b,xc,yc,theta,type)


if(type==1)%ellypse

ct=cos(theta);
st=sin(theta);    
A=ct^2*b^2+st^2*a^2;
B=2*ct*st*(b^2-a^2);
C=st^2*b^2+ct^2*a^2;
D=-2*xc*(b^2*ct^2+a^2*st^2)+2*yc*ct*st*(-b^2+a^2);
E=-2*yc*(a^2*ct^2+b^2*st^2)+2*xc*ct*st*(-b^2+a^2);
F=b^2*xc^2*ct^2+2*xc*yc*ct*st*b^2+yc^2*b^2*st^2+a^2*yc^2*ct^2-2*xc*yc*ct*st*a^2+xc^2*a^2*st^2-a^2*b^2;
C6=[A,B,C,D,E,F];

elseif(type==2)%hyperbola
    
ct=cos(theta);
st=sin(theta);        

A=b^2*ct^2-a^2*st^2;
B=2*st*ct*(a^2+b^2);
C=b^2*st^2-a^2*ct^2;
D=2*xc*(-b^2*ct^2+a^2*st^2)-2*st*ct*yc*(b^2+a^2);
E=2*yc*(-b^2*st^2+a^2*ct^2)-2*st*ct*xc*(b^2+a^2);
F=xc^2*b^2*ct^2+2*ct*st*xc*yc*b^2+b^2*st^2*yc^2-a^2*ct^2*yc^2+2*a^2*ct*st*xc*yc-xc^2*a^2*st^2-a^2*b^2;
C6=[A,B,C,D,E,F];

elseif(type==3)%parabola

ct=cos(theta);
st=sin(theta);

R = xc*st - yc*ct;

A=st^2;
B=-2*ct*st;
C=ct^2;
D=2*(-R*st-a*ct);
E=2*(R*ct-a*st);
F=R^2+2*a*(xc*ct+yc*st);

C6=[A,B,C,D,E,F];
end




end