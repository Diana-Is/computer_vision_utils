function bestFit=conics_RANSAC(XY,outl_ratio,if_plot)
%data – A set of observations 
%model – A model to explain observed data points.
%n – Minimum number of data points required to estimate model parameters.
%k – Maximum number of iterations allowed in the algorithm.
%t – Threshold value to determine data points that are fit well by model.
%d – Number of close data points required to assert that a model fits well to data.

%Return:
%bestFit – model parameters which best fit the data (or null if no good model is found)

n=5;
p=0.99;
e=outl_ratio;
k=log(1-p)/log(1-(1-e)^5);
t=5;%treshold
tm=5;

iter = 0;
bestFit = [];
bestErr =  1e+24;%something really large
tic
while (iter < k) 
    maybeInliers = XY(:,randperm(size(XY,2),5));
    [~,maybeModel]=C_estim(maybeInliers,1,0);
    alsoInliers = [];
    
    d=sampson_distance_conics(maybeModel,XY,'by_point','no');

    for i=1:size(XY,2)
        if d(i)<5
            alsoInliers=[alsoInliers,XY(:,i)];
        end
    end 

    if (size(alsoInliers,2) >= tm) 
        [~,betterModel]=C_estim(alsoInliers,1,0);
        thisErr = sampson_distance_conics(maybeModel,alsoInliers,'sum','no'); 
        if (thisErr < bestErr)
            bestFit= betterModel;
            bestErr= thisErr;
        end 
    end 
    iter=iter+1;
end 
toc
if(if_plot==1)
plot_a_conic(bestFit,XY,'RANSAC for Conic');
end


end