function [half_mat, half_mat_flatten, index_pair] = lower_half(matrix)

[n_rows, n_cols] = size(matrix);
if n_rows ~= n_cols
    error('Distance matrix should be a square matrix!')
else
    half_mat = matrix + tril(NaN*ones(n_rows,n_rows),-1);
    index_id = find(1-isnan(half_mat));
    [index_rows, index_cols] = ind2sub([n_rows, n_rows], index_id);
    index_pair = [index_rows(:), index_cols(:)];
    half_mat_flatten = half_mat;
    half_mat_flatten(isnan(half_mat_flatten)) = [];
end

end