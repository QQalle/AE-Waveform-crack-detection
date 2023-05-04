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
%Comments
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

% load the signals from file paths into a cell array
Signals = cellfun(@(x) load(x), FilePaths, 'UniformOutput', false);

% perform wavelet transform on each signal in the cell array
wt = cell(size(Signals)); % create an empty cell array to store the wavelet coefficients
f = cell(size(Signals)); % create an empty cell array to store the frequencies
for i = 1:length(Signals) % loop over the signals
    [wt{i},f{i}] = cwt(Signals{i},'amor',5e6,'FrequencyLimits',[200e3 1200e3]); % perform wavelet transform on each signal
end

%% 
maxAmp = [];
maxInd = [];
maxCol = [];
maxFreq = [];
for j = 1:length(Signals)
    [maxAmp(j), indtemp] = max(wt{j}, [], "all");
    [maxInd(j), maxCol(j)] = ind2sub(size(wt{j}), indtemp);
    maxFreq(j) = f{j}(maxInd(j));
end

plot3(maxFreq, maxCol, abs(maxAmp), ".");
xlabel("Freq.");
ylabel("index");
zlabel("Amp.");

% plot the wavelet coefficients in a waterfall plot
%figure; % create a new figure
%hold on; % hold on to plot multiple surfaces
%for i = 1:length(wt) % loop over the wavelet coefficients
%    waterfall(f{i},(1:length(Signals{i}))',abs(wt{i})'); % plot each surface
%end
%hold off; % release the hold