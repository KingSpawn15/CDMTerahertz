function data = read_refractive_index(filename)
    % Open the file
    fileID = fopen(filename, 'r');

    % Read the header line
    header = fgetl(fileID);

    % Initialize variables
    lambda_nm = [];
    n = [];
    k = [];

    % Read the data
    while ~feof(fileID)
        line = fgetl(fileID);
        C = strsplit(line, '\t'); % Modify separator if necessary

        % Extract the values
        lambda_nm = [lambda_nm; str2double(C{1})];
        n = [n; str2double(C{2})];
        k = [k; str2double(C{3})];
    end

    % Close the file
    fclose(fileID);

    % Store the data in a struct
    data = struct('lambda_nm', lambda_nm, 'n', n, 'k', k);
end