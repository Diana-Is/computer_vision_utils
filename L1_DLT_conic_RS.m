function s=L1_DLT_conic_RS(XY,if_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L1_DLT_conic_with RANSAC afterwards
%Input: The coordinates x,y, mu's initial value, eps -error tolerance
%Output: The estimated model s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1: Mapping correspondences into embeddings a i or b i to form the data matrix M.
XYs=[XY;ones(1,size(XY,2))];
[newXY, T] = precondition(XYs);
xw=newXY(1,:)';
yw=newXY(2,:)';
M=[xw.^2 xw.*yw yw.^2 xw yw ones(size(xw))];

% threshold value below which we consider an element to be zeros
delta = 1e-4;

%2: Initialize z as the right singular vector corresponding to the smallest singular value of M.
[~,~,V] = svd(M);
z=V(:,end);

k=0;%step counter
mu_init=0.5;
p=0.5;

f_k=norm(M*z);%first value in the array f_k(1)=M*z;
z_k(:)=z;%first value in the array z_k(1)=z;

while (k<3000) 
	g = M'*sign(M*z);%Compute sub-gradient: ;
	mu=mu_init*p^k; % Update the step size µ according to a certain rule.
	z = z - mu*g; %Sub-gradient descent
	z = z/norm(z,2); %Sphere projection: norm() returns the Euclidean norm of vector v. This norm is also called the 2-norm, vector magnitude, or Euclidean length.
	f_k(end+1) = norm(M*z);
  	z_k(end+1,:)=z;
	k=k+1;
	disp(k);
	if(abs(f_k(end)-f_k(end-1))<=delta || norm(z_k(end-1,:)-z_k(end,:))<=delta)%convergence criteria
		disp('stopped before 3000 iterations because fulfilled convergence criteria')
		break;
	end

end

[f_best,ind] = min(f_k);
s=z_k(ind,:);

[A,B,C,D,E,F]=deal(s(1),s(2),s(3),s(4),s(5),s(6));
C_M_est=[A,B/2,D/2; B/2,C,E/2; D/2,E/2,F];

% denormalise C
C_M_est = T'*C_M_est*T;

if(if_plot==1)
plot_a_conic(C_M_est,XY,'L1 DLT conic estimation');
end

%Compute the residuals of each correspondence with respect to z.
%so, compute the distance from the conic determined to the all other sets

%Determine the inlier set I by thresholding the residuals.

%Post-processing on I to obtain final estimation result s. -- RANSAC


end 