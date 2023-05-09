figure('name', 'Stress vs Time','Position',[60,60,1400,700])
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_TensileStress)
title('Stress vs Time');
xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel('Stress [MPa]');