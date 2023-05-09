figure('name', 'Energy-Frequency Spectrum','Position',[60,60,1400,700])
hold on
plot((1:length(StackEner)),StackEner);
xlim([0 1000]); 
ylim([0 max(StackEner)]);
title('Energy-Frequency Spectrum');
xlabel('Peak frequency [kHz]');
ylabel('Energy [aJ]');
xl1 = xline(MCminFreq/1000);
xl1.Color = 'r';
xl2 = xline(MCmaxFreq/1000);
xl2.Color = 'r';
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PeakFrequency, Matrixcracks.Energy, 'o');
end
hold off