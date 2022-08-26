function [C_it_est] = nonlin_fitting_conics(P,if_display,if_plot)

% The format of the xs is
% [x1 x2 x3 ... xn ; 
%  y1 y2 y3 ... yn ;
%  w1 w2 w3 ... wn]


[r,c] = size(P);

if ( r<2 || c<5)
 error ('Error:wrong input data')
end

if (size(P,1) == 2)
  PP=[P ; ones(1,size(P,2))]; %making homogenious coordinates
  %no preconditioning, it makes wrong results
end

%initial DLT conic fitting
[C_v_est,C_M_est]=C_estim(P,1,0);

if(if_display==1)
opt = optimset( optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','on', 'Display','on');
elseif(if_display==0)
opt = optimset( optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
end

% opt = optimset( optimset('lsqnonlin') , 'LargeScale','off', 'Diagnostics','off', 'Display','off');
C_it_est = lsqnonlin(@(C_v)lsq_func(C_v,PP),C_v_est,[],[],opt);
%[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(fun,x0,lb,ub,options)
[A,B,C,D,E,F]=deal(C_it_est(1),C_it_est(2),C_it_est(3),C_it_est(4),C_it_est(5),C_it_est(6));
C_M=[A B/2 D/2; B/2 C E/2; D/2 E/2 F];

%rms = sampson_distance_conics(C_M,PP,'sum','yes'); %as a potential output parameter

if(if_plot==1)
plot_a_conic(C_M,P,'HZ algorithm, L-M,Samson');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = lsq_func(C_v,PP)

  [A,B,C,D,E,F]=deal(C_v(1),C_v(2),C_v(3),C_v(4),C_v(5),C_v(6));
  C_M=[A B/2 D/2; B/2 C E/2; D/2 E/2 F];
  r = sampson_distance_conics(C_M,PP,'sum','yes');

end
  
