figure('name', 'Frequency_Duration_Energy','Position',[60,60,1400,700])
hold on
scatter(PFreqList./1000, HAFImpDurList, 30, HAFImpEnerList, '.');
colorbar
title('Frequency Duration Energy');
xlabel('Peak frequency [kHz]');
ylabel('Duration [μs]');