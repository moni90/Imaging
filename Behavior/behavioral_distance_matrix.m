function [p_no_go_matrix, distance_mat] = behavioral_distance_matrix(data)
%This function computes the behavioral distance matrix using the definition
%in the biorXiv paper

n_odors = 24;
count_matrix = zeros(n_odors, n_odors);
no_go_matrix = zeros(n_odors, n_odors);

for i = 1:length(data)
    row = data(i).Odor1_id;
    col = data(i).Odor2_id;
    if row == col
        count_matrix(row, col) = count_matrix(row ,col)+1;
        no_go_matrix(row, col) = no_go_matrix(row, col) + (1-data(i).go_trial);
    else
        count_matrix(row, col) = count_matrix(row, col)+1;
        count_matrix(col, row) = count_matrix(col, row)+1;
        no_go_matrix(row, col) = no_go_matrix(row, col) + (1-data(i).go_trial);
        no_go_matrix(col, row) = no_go_matrix(col, row) + (1-data(i).go_trial);
    end
end

p_no_go_matrix = no_go_matrix./count_matrix;
p_no_go_AA = repmat(diag(p_no_go_matrix),1,n_odors);
p_no_go_BB = repmat(diag(p_no_go_matrix)',n_odors,1);
distance_mat = 1 - (1 - p_no_go_matrix)./sqrt((1-p_no_go_AA).*(1-p_no_go_BB));
distance_mat = distance_mat.*(distance_mat>=0);