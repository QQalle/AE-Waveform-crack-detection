
CheckAvgFreq = round(mean(CheckVariableTable.PeakFrequency));
CheckAvgAmp = round(mean(CheckVariableTable.Amplitude));
CheckAvgDur = round(mean(CheckVariableTable.Duration));
CheckAvgEner = round(mean(CheckVariableTable.Energy));
CheckAvgCount = round(mean(CheckVariableTable.Counts));
CheckAvgRise = round(mean(CheckVariableTable.RiseTime));

CheckSTDFreq = round(3*std(CheckVariableTable.PeakFrequency));
CheckSTDAmp = round(3*std(CheckVariableTable.Amplitude));
CheckSTDDur = round(3*std(CheckVariableTable.Duration));
CheckSTDEner = round(3*std(CheckVariableTable.Energy));
CheckSTDCount = round(3*std(CheckVariableTable.Counts));
CheckSTDRise = round(3*std(CheckVariableTable.RiseTime));

disp("CHECKED STATISTICS MEAN ± 3σ");
disp("Avg peak frequency: " + num2str(CheckAvgFreq) + "kHz ± " + num2str(CheckSTDFreq)); 
disp("Avg amplitude: " + num2str(CheckAvgAmp) + "dB ± " + num2str(CheckSTDAmp)); 
disp("Avg duration: " + num2str(CheckAvgDur) + "μs ± " + num2str(CheckSTDDur)); 
disp("Avg energy: " + num2str(CheckAvgEner) + "aJ ± " + num2str(CheckSTDEner)); 
disp("Avg counts: " + num2str(CheckAvgCount) + " ± " + num2str(CheckSTDCount)); 
disp("Avg rise time: " + num2str(CheckAvgRise) + "s ± " + num2str(CheckSTDRise)); 

Stats = table();
Stats(1,:) = table([CheckAvgFreq CheckSTDFreq],[CheckAvgAmp CheckSTDAmp],...
    [CheckAvgDur CheckSTDDur],[CheckAvgEner CheckSTDEner],...
    [CheckAvgCount CheckSTDCount],[CheckAvgRise CheckSTDRise]);

if istable(Matrixcracks)
        %Compute statistics
    MCAvgFreq = round(mean(Matrixcracks.PeakFrequency));
    MCAvgAmp = round(mean(Matrixcracks.Amplitude));
    MCAvgDur = round(mean(Matrixcracks.Duration));
    MCAvgEner = round(mean(Matrixcracks.Energy));
    MCAvgCount = round(mean(Matrixcracks.Counts));
    MCAvgRise = round(mean(Matrixcracks.RiseTime));

    MCSTDFreq = round(3*std(Matrixcracks.PeakFrequency));
    MCSTDAmp = round(3*std(Matrixcracks.Amplitude));
    MCSTDDur = round(3*std(Matrixcracks.Duration));
    MCSTDEner = round(3*std(Matrixcracks.Energy));
    MCSTDCount = round(3*std(Matrixcracks.Counts));
    MCSTDRise = round(3*std(Matrixcracks.RiseTime));
    
    disp("MATRIX CRACK STATISTICS MEAN ± 3σ");
    disp("Avg peak frequency: " + num2str(MCAvgFreq) + "kHz ± " + num2str(MCSTDFreq)); 
    disp("Avg amplitude: " + num2str(MCAvgAmp) + "dB ± " + num2str(MCSTDAmp)); 
    disp("Avg duration: " + num2str(MCAvgDur) + "μs ± " + num2str(MCSTDDur)); 
    disp("Avg energy: " + num2str(MCAvgEner) + "aJ ± " + num2str(MCSTDEner)); 
    
    Stats(end+1,:) = table([MCAvgFreq MCSTDFreq],[MCAvgAmp MCSTDAmp],...
    [MCAvgDur MCSTDDur],[MCAvgEner MCSTDEner],...
    [MCAvgCount MCSTDCount],[MCSTDRise MCSTDRise]);
end
if istable(Debondings)
    DBAvgFreq = round(mean(Debondings.PeakFrequency));
    DBAvgAmp = round(mean(Debondings.Amplitude));
    DBAvgDur = round(mean(Debondings.Duration));
    DBAvgEner = round(mean(Debondings.Energy));
    DBAvgCount = round(mean(Debondings.Counts));
    DBAvgRise = round(mean(Debondings.RiseTime));
    
    DBSTDFreq = round(3*std(Debondings.PeakFrequency));
    DBSTDAmp = round(3*std(Debondings.Amplitude));
    DBSTDDur = round(3*std(Debondings.Duration));
    DBSTDEner = round(3*std(Debondings.Energy));
    DBSTDCount = round(3*std(Debondings.Counts));
    DBSTDRise = round(3*std(Debondings.RiseTime));
    
    disp("DEBONDING STATISTICS MEAN ± 3σ");
    disp("Avg peak frequency: " + num2str(DBAvgFreq) + "kHz ± " + num2str(DBSTDFreq)); 
    disp("Avg amplitude: " + num2str(DBAvgAmp) + "dB ± " + num2str(DBSTDAmp)); 
    disp("Avg duration: " + num2str(DBAvgDur) + "μs ± " + num2str(DBSTDDur)); 
    disp("Avg energy: " + num2str(DBAvgEner) + "aJ ± " + num2str(DBSTDEner)); 
    
    Stats(end+1,:) = table([DBAvgFreq DBSTDFreq],[DBAvgAmp DBSTDAmp],...
    [DBAvgDur DBSTDDur],[DBAvgEner DBSTDEner],...
    [DBAvgCount DBSTDCount],[DBSTDRise DBSTDRise]);
end
Stats.Properties.VariableNames = {'MEAN | Frequency | 3σ' 'Amplitude' 'Duration'...
    'Energy' 'Counts' 'Rise time'};

if istable(Matrixcracks) %Check what rows to name
    if istable(Debondings)
        Stats.Properties.RowNames = {'Check' 'Matrix cracks' 'Debondings'};
    else
        Stats.Properties.RowNames = {'Check' 'Matrix cracks'};
    end
elseif istable(Debondings)
    Stats.Properties.RowNames = {'Check' 'Debondings'};
else
    Stats.Properties.RowNames = {'Check'};
end
