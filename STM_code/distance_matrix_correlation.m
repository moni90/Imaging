function dist_corr = distance_matrix_correlation(matrix1, matrix2)

%extract only half matrix
[mat1_half_2d, mat1_half_1d, index_pair] = lower_half(matrix1);
[mat2_half_2d, mat2_half_1d, index_pair] = lower_half(matrix2);
dist_corr = corr(mat1_half_1d(:), mat2_half_1d(:));
end


    

