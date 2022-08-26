function iter_par_out = my_IRLS(xp,yp,iters,if_diag)
[A,~,y]=reshape_problem(xp,yp);

t=1;
W=eye(size(A,1));
delta=0.0001;
iter_par(:,1)=zeros(5,1);

while (t<iters)
    if(det(A'*W*A)<1.000e-10)
    disp('coord')
    %disp([xp,yp])
    %plot_a_conic(C_Mat,[xp';yp'],'ground truth, generated');
    break
    end
    iter_par(:,t+1) = (A'*W*A)^(-1)*A'*W*y; % weighted L_2 sol.
    %iter_par=iter_par/norm([iter_par;1]); %sphere projection %not necess because we are not doing M*z=0
    %and avoiding trivial solution

    t=t+1;
    
    if (if_diag==1)
        disp(t);
    end
    
    
    if(t==1)
    W=eye(size(A,1));
    else
        for i=1:size(A,1)
            W(i,i)=1/max(delta,abs(y(i)-A(i,:)*iter_par(:,t)));
        end
    end
    
    if(norm((A*iter_par-y),1)<=1e-24 || norm((iter_par(:,t-1)-iter_par(:,t)),1)<=1e-24)%convergence criteria
        if (if_diag==1)
        disp('stopped before 3000 iterations because fulfilled convergence criteria')
        end
        break;
    end

end
iter_par_out=iter_par(:,end);

end

function [A,par,rhs]=reshape_problem(x,y)
%Q(x,y)=Ax^2+Bxy+Cy^2+Dx+Ey+F=0

A=[sum(y.^4),sum((x.*y).^2),sum(x.*(y.^3)),sum(y.^3),sum(x.*(y.^2));
sum((x.*y).^2),sum(x.^4),sum(y.*(x.^3)),sum(y.*(x.^2)),sum(x.^3);
sum(x.*(y.^3)),sum(y.*(x.^3)),sum((x.*y).^2),sum(x.*(y.^2)),sum(y.*(x.^2));
sum(y.^3),sum(y.*(x.^2)),sum(x.*(y.^2)),sum(y.^2),sum(x.*y);
sum(x.*(y.^2)),sum(x.^3),sum(y.*(x.^2)),sum(x.*y),sum(x.^2)];

rhs=[-sum(y.^2);-sum(x.^2);-sum(x.*y);-sum(y);-sum(x)];

%first estimation of the parameter vector
par=pinv(A)*rhs;

end