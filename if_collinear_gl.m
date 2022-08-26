function tf = if_collinear_gl(pts)
mat = [pts(1,1)-pts(3,1),pts(1,2)-pts(3,2); ...
      pts(2,1)-pts(3,1),pts(2,2)-pts(3,2)];
det(mat)
tf=(det(mat)>-0.1)&&(det(mat)<0.1);
end
