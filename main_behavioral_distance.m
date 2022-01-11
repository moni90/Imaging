data_path = 'D:\mmoroni\Hiro data\Behavior\Odor_set1';
code_path = 'D:\mmoroni\github_repos\Imaging';
addpath(genpath(code_path));
save_path = 'D:\mmoroni\Hiro_project\analyses\Behavior\Odor_set1\';

folder_list = dir(data_path);

for id_folder = 3:length(folder_list)
    if exist(fullfile(save_path, [folder_list(id_folder).name, '_behavioral_data.mat']))
        load(fullfile(save_path, [folder_list(id_folder).name, '_behavioral_data.mat']))
    else
        cd(fullfile(data_path, folder_list(id_folder).name));
        filename_list = dir('**\*.h5');
        cd(code_path);
        
        for i = 1:length(filename_list)
            disp(i)
            filename = fullfile(filename_list(i).folder, filename_list(i).name);
            data_temp = read_behavioral_data(filename);
            if i == 1
                data = data_temp;
            else
                data = [data, data_temp];
            end
        end
        save(fullfile(save_path, [folder_list(id_folder).name, '_behavioral_data.mat']), 'data');
        
    end
    if exist(fullfile(save_path, [folder_list(id_folder).name, '_distance_matrix.mat']))
        load(fullfile(save_path, [folder_list(id_folder).name, '_distance_matrix.mat']))
    else
        [p_no_go_matrix_temp, distance_mat_temp] = behavioral_distance_matrix(data);
        save(fullfile(save_path, [folder_list(id_folder).name, '_distance_matrix.mat']), 'p_no_go_matrix_temp', 'distance_mat_temp');
    end
    figure;
    subplot(1,2,1);imagesc(p_no_go_matrix_temp); title(['Mouse ',folder_list(id_folder).name, '. P(no-go)']); colormap(gray); caxis([0, 1]); colorbar;
    subplot(1,2,2); imagesc(distance_mat_temp); title(['Mouse ', folder_list(id_folder).name, '. Behavioral distance']); colormap(gray); caxis([0, 1]); colorbar;
    
    if id_folder == 3
        data_pooled = data;
    else
        data_pooled = [data_pooled, data];
    end
end

[p_no_go_matrix, distance_mat] = behavioral_distance_matrix(data_pooled);
figure;
subplot(1,2,1);imagesc(p_no_go_matrix); title('P(no-go)'); colormap(gray); caxis([0, 1]); colorbar;
subplot(1,2,2); imagesc(distance_mat); title('Behavioral distance'); colormap(gray); caxis([0, 1]); colorbar;

save(fullfile(save_path, 'pooled_behavioral_data.mat'), 'data_pooled');
save(fullfile(save_path, 'pooled_distance_matrix.mat'), 'p_no_go_matrix', 'distance_mat');