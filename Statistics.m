
if istable(Matrixcracks)
        %Compute statistics
    MCAvgFreq = round(mean(Matrixcracks.PeakFrequency),1);
    MCAvgAmp = round(mean(Matrixcracks.Amplitude),1);
    MCAvgDur = round(mean(Matrixcracks.Duration),1);
    MCAvgEner = round(mean(Matrixcracks.Energy),1);

    MCDevFreq = max(Matrixcracks.PeakFrequency) - round(MCAvgFreq,1);
    MCDevAmp = max(Matrixcracks.Amplitude) - round(MCAvgAmp,1);
    MCDevDur = max(Matrixcracks.Duration) - round(MCAvgDur,1);
    MCDevEner = max(Matrixcracks.Energy) - round(MCAvgEner,1);
    
    disp("MATRIX CRACK STATISTICS");
    disp("Avg peak frequency: " + num2str(MCAvgFreq) + "kHz ± " + num2str(MCDevFreq)); 
    disp("Avg amplitude: " + num2str(MCAvgAmp) + "dB ± " + num2str(MCDevAmp)); 
    disp("Avg duration: " + num2str(MCAvgDur) + "μs ± " + num2str(MCDevDur)); 
    disp("Avg energy: " + num2str(MCAvgEner) + "aJ ± " + num2str(MCDevEner)); 
end
if istable(Debondings)
    DBAvgFreq = round(mean(Debondings.PeakFrequency),1);
    DBAvgAmp = round(mean(Debondings.Amplitude),1);
    DBAvgDur = round(mean(Debondings.Duration),1);
    DBAvgEner = round(mean(Debondings.Energy),1);
    
    DBDevFreq = max(Debondings.PeakFrequency) - round(DBAvgFreq,1);
    DBDevAmp = max(Debondings.Amplitude) - round(DBAvgAmp,1);
    DBDevDur = max(Debondings.Duration) - round(DBAvgDur,1);
    DBDevEner = max(Debondings.Energy) - round(DBAvgEner,1);
    
    disp("DEBONDING STATISTICS");
    disp("Avg peak frequency: " + num2str(DBAvgFreq) + "kHz ± " + num2str(DBDevFreq)); 
    disp("Avg amplitude: " + num2str(DBAvgAmp) + "dB ± " + num2str(DBDevAmp)); 
    disp("Avg duration: " + num2str(DBAvgDur) + "μs ± " + num2str(DBDevDur)); 
    disp("Avg energy: " + num2str(DBAvgEner) + "aJ ± " + num2str(DBDevEner)); 
end

CheckAvgFreq = round(mean(CheckVariableTable.PeakFrequency),1);
CheckAvgAmp = round(mean(CheckVariableTable.Amplitude),1);
CheckAvgDur = round(mean(CheckVariableTable.Duration),1);
CheckAvgEner = round(mean(CheckVariableTable.Energy),1);
CheckAvgCount = round(mean(CheckVariableTable.Counts),1);
CheckAvgRise = round(mean(CheckVariableTable.RiseTime),1);

CheckDevFreq = max(CheckVariableTable.PeakFrequency)...
    - round(CheckAvgFreq,1);
CheckDevAmp = max(CheckVariableTable.Amplitude)...
    - round(CheckAvgAmp,1);
CheckDevDur = max(CheckVariableTable.Duration)...
    - round(CheckAvgDur,1);
CheckDevEner = max(CheckVariableTable.Energy)...
    - round(CheckAvgEner,1);
CheckDevCount = max(CheckVariableTable.Counts)...
    - round(CheckAvgCount,1);
CheckDevRise = max(CheckVariableTable.RiseTime)...
    - round(CheckAvgRise,1);

disp("CHECKED STATISTICS");
disp("Avg peak frequency: " + num2str(CheckAvgFreq) + "kHz ± " + num2str(CheckDevFreq)); 
disp("Avg amplitude: " + num2str(CheckAvgAmp) + "dB ± " + num2str(CheckDevAmp)); 
disp("Avg duration: " + num2str(CheckAvgDur) + "μs ± " + num2str(CheckDevDur)); 
disp("Avg energy: " + num2str(CheckAvgEner) + "aJ ± " + num2str(CheckDevEner)); 
disp("Avg counts: " + num2str(CheckAvgCount) + " ± " + num2str(CheckDevCount)); 
disp("Avg rise time: " + num2str(CheckAvgRise) + "s ± " + num2str(CheckDevRise)); 

