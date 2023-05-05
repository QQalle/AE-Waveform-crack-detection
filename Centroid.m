close all
clear

experimentNo = '1003'; %Specify which experiment to analize
ASCIIOutPut = importdata(append('Data/EXP', experimentNo, '.txt'));
ASCIIWaveforms = append('Data/EXP', experimentNo);

%Import waveform
filePattern = fullfile(ASCIIWaveforms, '*.txt');
TheFiles = dir(filePattern);

fileNameArray = append({TheFiles.folder}, {TheFiles.name});

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

FilePaths = transpose(fullfile({TheFiles.folder}, {TheFiles.name}));

Signals = cellfun(@(x) load(x), FilePaths, 'UniformOutput', false);
disp("Signals Done!")

L = numel(Signals);

FFTs = cellfun(...
    @(s) (fft(s - mean(s)/L)), Signals, 'UniformOutput', false);
disp("FFT Done!")

Fs = 10*10^6;
Fn = Fs/2;
N = numel(FFTs);
Fv = linspace(0, 1, fix(L/2)+1)*Fn;
Iv = 1:numel(Fv);

Centroids = cellfun(@(F) (spectralCentroid(abs(F(Iv))*2, Fv)), FFTs, 'UniformOutput', false);
disp("Centroids Done!");

% plot(cell2mat(Centroids));
% waterfall(abs(FFTs{1}))



