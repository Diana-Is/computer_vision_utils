function d=sampson_distance_conics(C,P,by_point,if_sqrt)

%Sampson distance or first order approximation of the geometric error
%INPUT: matrix of the conic C, matrix of input points P,
%by_point -- if output should be by point for thresholding or sum
%if_sqrt -- if from the resulting distance square root is taken
%OUTPUT: sum over all points of Sampson distance of each point to the conic
%
%
%                        (point'*C*point)^2
% d(point)=    --------------------------------------
%              (2*(C*point)(1))^2+(2*(C*point)(2))^2


if(size(P,1)==2)
P=[P;ones(1,size(P,2))];
end

eps=(alg_dist_conics(C,P))';  %nominator

J1=[];
J2=[];


%P(:,i) is a column vector as in the book
for i=1:size(P,2)
J1(i)=P(:,i)'*C*[1;0;0]; %denominator 1 term
J2(i)=P(:,i)'*C*[0;1;0]; %denominator 2 term
end

dd=[];

for i=1:size(P,2)

	if strcmp(if_sqrt,'yes')
	dd(i)=sqrt(0.25*eps(i)^2/(J1(i)^2+J2(i)^2));

	elseif strcmp(if_sqrt,'no')

	dd(i)=0.25*eps(i)^2/(J1(i)^2+J2(i)^2);
	else
	error ('Wrong if_sqrt parameter!');
	end

end


if strcmp(by_point,'by_point')
d=dd;
elseif strcmp(by_point,'sum')
d=sum(dd);
else
error ('Wrong by_point parameter!');
end

end