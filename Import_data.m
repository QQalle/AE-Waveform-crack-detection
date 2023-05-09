    %Define
disp(append('Loading Data: ',num2str(experimentNo)));
ASCIIOutPut = importdata(append('Data/EXP', experimentNo, '.txt'));
ASCIIWaveforms = append('Data/EXP', experimentNo);
CSVRawData = readtable(append('Data/Specimen_RawData_',...
    experimentNo,'.csv'),"Delimiter",';');
CSVData = varfun(@(x) str2double(replace(x, ",", ".")), CSVRawData);
CSVData = CSVData(2:end,:);
