figure('name', 'Estimate cracks','Position',[60,60,1400,700])
hold on
for k = 1 : exp
    plot(SV.time2{k}, smooth(SV.EstCracksTotal{k},200));
    p = plot(SV.time2{k}, SV.EstCracksTotal{k},'--');
end
title(append('Estimated Matrix Cracks vs Time',captext));
xlabel('Time [s]');
ylabel('Matrix cracks');
legend(Legendtext,'location','east outside');
grid on
% yticks((0:1:15)*10^7)
% xticks(0:50:800)
set(gca,'FontSize',14)
hold off