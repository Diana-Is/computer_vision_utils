function plotting_using_C(C)

im = zeros(500, 500);
for ii = 1 : 500
    for jj = 1 : 500
        im(ii, jj) = [jj, ii, 1] * C * [jj; ii; 1];
    end
end

figure 
imshow(im);

end