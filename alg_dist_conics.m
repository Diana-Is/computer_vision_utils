function alg_d=alg_dist_conics(C,P)

%function, computing algebraic distance, x'Cx, where x is the coordinates
%of each point
%C -- matrix of conic coefficients
%P -- matrix of points coordinates

if(size(P,1)==2)
P=[P;ones(1,size(P,2))];
end
alg_d=zeros(1,size(P,2));
for i=1:size(P,2)
alg_d(i)=P(:,i)'*C*P(:,i);
end

end