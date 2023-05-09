figure('name', 'Amplitude-Frequency','Position',[60,60,1400,700])
hold on
plot(PFreqList./1000, HAFImpAmpList, '.');
xlim([floor(min(PFreqList./1000)/10)*10 ...
    ceil(max(PFreqList./1000)/10)*10]); %nearest 10-number of max sample
ylim([floor(min(HAFImpAmpList)/10)*10 ...
    ceil(max(HAFImpAmpList)/10)*10]);
title('Amplitude-Frequency');
xlabel('Peak frequency [kHz]');
ylabel('Amplitude [dB]');
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PeakFrequency, Matrixcracks.Amplitude, 'o');
    legend(plus("Hits: ",num2str(K)),...
        plus("Matrix cracks: ",...
        num2str(length(Matrixcracks.HitIndex))),...
        'location','south outside',"AutoUpdate","off");
else
    legend(plus("Hits: ",num2str(K)),...
        plus("Matrix cracks: ",...
        num2str(length(Matrixcracks))),'location','south outside',"AutoUpdate","off");
end
patch([MCminFreq/1000 MCminFreq/1000 MCmaxFreq/1000 MCmaxFreq/1000],...
    [MCminAmp MCmaxAmp MCmaxAmp MCminAmp],'r','FaceAlpha',0,...
    'EdgeColor','r');
plot(CheckVariableTable.PeakFrequency, CheckVariableTable.Amplitude,'*');
hold off