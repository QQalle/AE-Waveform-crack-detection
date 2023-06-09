    %Define
disp(append('    Loading Data: ',num2str(experimentNo)));
ASCIIOutPut = importdata(append('Data/EXP', experimentNo, '.txt'));
CSVRawData = readtable(append('Data/Specimen_RawData_',...
    experimentNo,'.csv'),"Delimiter",';');
CSVData = varfun(@(x) str2double(replace(x, ",", ".")), CSVRawData);
CSVData = CSVData(2:end,:);

    %Import waveform
ASCIIWaveforms = append('Data/EXP', experimentNo);
filePattern = fullfile(ASCIIWaveforms, '*.txt');
TheFiles = dir(filePattern);

    %Sort the files
C = transpose({TheFiles.name});
% Extract indices from each string using regular expressions
filesindices = cellfun(@(x) regexp(x, '(?<=_)\d+(?=_\d+\.txt)',...
    'match'), C, 'UniformOutput', false);
% Convert indices from cell array of strings to numeric array
filesindices = cellfun(@(x) str2double(x), filesindices, ...
    'UniformOutput', false);
% Sort the indices
[sorted_indices, order] = sort([filesindices{:}]);
% Reorder the original cell array using the sorted indices
TheFiles = TheFiles(order);

baseFN = TheFiles(end).name;
%baseFileName = FilesSorted(end);
fullFN = fullfile(TheFiles(end).folder, baseFN);
EndTime = str2double(regexp(fullFN,...
        '(?<=_)\d+(?=\.txt)', 'match','once'))/1000000; %Find time of hit
EndTime = ceil(EndTime);