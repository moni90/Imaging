function mat_sh = shuffle_symmetric_matrix(mat0)

[n_rows, n_cols] = size(matrix);
if n_rows ~= n_cols
    error('Distance matrix should be a square matrix!')
else
    mat_sh = zeros(size(mat0));
    [half_mat, half_mat_flatten, index_pair] = lower_half(mat0);
    half_mat_flatten_sh = half_mat_flatten(randperm(length(half_mat_flatten)));
    for i = 1:size(index_pair,1)
        mat_sh(index_pair(i,1),index_pair(i,2)) = half_mat_flatten_sh(i);
        if index_pair(i,1) ~= index_pair(i,2)
            mat_sh(index_pair(i,2),index_pair(i,1)) = half_mat_flatten_sh(i);
        end
    end
end
    

