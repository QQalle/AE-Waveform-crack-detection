figure('name','Amplitude vs Parameter 1','Position',[60,60,1400,700])
hold on
scatter(HAFImpPARA1, HAFImpAmpList, 60, '.');
title('Amplitude vs Parameter 1');
ylim([floor(min(HAFImpAmpList))*0.9 ceil(max(HAFImpAmpList))*1.1]);
xlabel('PARA1 [?]');
ylabel("Amplitude [dB]");
if istable(Matrixcracks) %If true = there are matrixcracks
    plot(Matrixcracks.PARA1, Matrixcracks.Amplitude, 'o');
end
hold off