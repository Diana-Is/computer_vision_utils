% % function plotting_averaging_results(data_matrix,my_title,type)
% % if strcmp(type,'t')
% % my_data=mean(data_matrix,3);
% % [X,Y] = meshgrid(0:0.1:1,0:0.1:0.90);
% % figure, surf(X,Y,my_data)
% % xlabel('noise level percent of standard'), ylabel('outlier ratio'), zlabel('failure rate')
% % xlim([0 1]);
% % ylim([0 1]);
% % zlim([0 1]);
% % title(my_title)
% % elseif strcmp(type,'d1')
% % my_data=mean(data_matrix,3);
% % [X,Y] = meshgrid(0:0.1:1,0:0.1:0.90);
% % figure, surf(X,Y,my_data)
% % xlabel('noise level percent of standard'), ylabel('outlier ratio'), zlabel('average sampson dist')
% % xlim([0 1]);
% % ylim([0 1]);
% % title(my_title)
% % elseif strcmp(type,'d2')
% % 
% % flrs=zeros(10,11,100);
% %     
% % for i=1:1:10
% %     for j=1:1:11
% %         for k=1:1:100
% %         if data_matrix(i,j,k)>0.2
% %             flrs(i,j,k)=1;
% %         end
% %         end
% %     end
% % 
% % end
% % 
% % failure_percent=mean(flrs,3);
% % 
% % [X,Y] = meshgrid(0:0.1:1,0:0.1:0.90);
% % figure, surf(X,Y,failure_percent)
% % xlabel('noise level percent of standard'), ylabel('outlier ratio'), zlabel('failure rate')
% % xlim([0 1]);
% % ylim([0 1]);
% % title(my_title)
% % 
% % end
% % end


function plotting_averaging_results(data_matrix,my_title,type)
if strcmp(type,'t')
my_data=mean(data_matrix,3);
[X,Y] = meshgrid(0:0.1:1,5:5:20);
figure, surf(X,Y,my_data)
xlabel('noise level percent of standard'), ylabel('inlier number'), zlabel('failure rate')
xlim([0 1]);
ylim([5 20]);
zlim([0 1]);
title(my_title)
elseif strcmp(type,'d1')
my_data=mean(data_matrix,3);
[X,Y] = meshgrid(0:0.1:1,5:5:20);
figure, surf(X,Y,my_data)
xlabel('noise level percent of standard'), ylabel('inlier number'), zlabel('average sampson dist from GT model')
xlim([0 1]);
ylim([5 20]);
zlim([0 1]);
title(my_title)

elseif strcmp(type,'d2')

flrs=zeros(10,11,100);
    
for i=1:1:10
    for j=1:1:11
        for k=1:1:100
        if data_matrix(i,j,k)>0.2
            flrs(i,j,k)=1;
        end
        end
    end

end

failure_percent=mean(flrs,3);

[X,Y] = meshgrid(0:0.1:1,0:0.1:0.90);
figure, surf(X,Y,failure_percent)
xlabel('noise level percent of standard'), ylabel('outlier ratio'), zlabel('failure rate')
xlim([0 1]);
ylim([0 1]);
title(my_title)

end
end