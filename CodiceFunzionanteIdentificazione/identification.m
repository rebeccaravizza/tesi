clear; clc; close all;

registry_file = "execution_registry.txt";
% Registry file check
if exist(registry_file, 'file') == 0
    fileID = fopen(registry_file, 'w');
    fclose(fileID);
    disp('Execution registry file created.');
end

registry = fileread(registry_file);

% parent folder is the path to the folder containing preprocessed data
parentFolder = 'preprocessedData';
subfolders = dir(parentFolder);

% save path is the path to the folder where you want to save data
savePath = 'IdentificationData';

mkdir(savePath);

% Iterate over each item in the parent folder
for i = 1:numel(subfolders)
    if subfolders(i).isdir && startsWith(subfolders(i).name, 'Test')
        subfolderPath = fullfile(parentFolder, subfolders(i).name);
        
        % Create save folder
        saveFolder = fullfile(savePath, subfolders(i).name);
        mkdir(saveFolder);
        
        % Retrieve a list of files in the subfolder
        files = dir(subfolderPath);
        
        
        % Iterate over each file in the subfolder
        for j = 1:numel(files)
            if ~files(j).isdir && startsWith(files(j).name, 'data')
                filePath = fullfile(subfolderPath, files(j).name);
                if contains(registry, files(j).name)
                    disp(["File", files(j).name, "already processed"]);
                else
                    outputPath = fullfile(saveFolder, strrep(files(j).name, '.mat', ''));
                    mkdir(outputPath)
                    % Process the file here
                    identification_fn(filePath, outputPath, 60000);
                    fileID = fopen(registry_file, 'a');
                    fprintf(fileID, '%s\n', files(j).name);
                    fclose(fileID);
                end
            end
        end
    end
end