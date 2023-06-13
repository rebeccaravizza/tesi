parentFolder = 'RawData';
subfolders = dir(parentFolder);
savePath = 'PreprocessedData';

mkdir 'PreprocessedData';

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
            if ~files(j).isdir && startsWith(files(j).name, 'Test')
                filePath = fullfile(subfolderPath, files(j).name);
                outputPath = fullfile(saveFolder, strrep(files(j).name, 'Test', 'data_filt'));
                % Process the file here
                preprocess_fn(filePath, outputPath)
            end
        end
    end
end