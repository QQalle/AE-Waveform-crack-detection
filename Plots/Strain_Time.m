figure('name', 'Strain vs Time','Position',[60,60,1400,700])
plot(CSVDataOffs.Fun_Time, CSVDataOffs.Fun_TensileStrain*100)
title('Strain vs Time');
xlim([0 TimeEnd]);
xlabel('Time [s]');
ylabel('Strain [%]');