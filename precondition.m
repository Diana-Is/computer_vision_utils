function [newXY, T] = precondition(XY)

    if size(XY,1) ~= 3
        error('XY must be 3xN');
    end
    
    % For the finite points ensure homogeneous coords have scale of 1
    XY(1,:) = XY(1,:)./XY(3,:);
    XY(2,:) = XY(2,:)./XY(3,:);
    XY(3,:) = 1;
    
    centr = mean(XY(1:2,:)')';            % Centroid of finite points
    temp_pts(1,:) = XY(1,:)-centr(1); % Shift origin to centroid.
    temp_pts(2,:) = XY(2,:)-centr(2);
    
    dist = sqrt(temp_pts(1,:).^2 + temp_pts(2,:).^2);
    meandist = mean(dist(:));  
    
    scale = sqrt(2)/meandist;
    
    T = [scale   0   -scale*centr(1)
         0     scale -scale*centr(2)
         0       0      1      ];
    
    newXY = T*XY;

end

% % tx = mean(X(1, :));
% % ty = mean(X(2, :));
% % 
% % s = mean(std(X, [], 2));
% % 
% % T = [1/s, 0, -tx/s; 0, 1/s, -ty/s; 0, 0, 1];
% % 
% % XCond = T*X;


