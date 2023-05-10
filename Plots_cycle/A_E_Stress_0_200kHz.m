figure('name', 'A-E vs Stress 0-200kHz','Position',...
    [60,60,1400,700])
hold on
for k = 1 : exp
%     yyaxis left
%     plot(SV.time2{k}/SV.Resolution(k), SV.SpreadEner{k});
    plot(SV.Fun_TensileStress{k}, SV.SpreadEner1{k});
%     yyaxis right
%     plot(SV.Fun_Time{k}, SV.Fun_TensileStress{k})
end
title('Cumulative Acoustic Energy vs Stress 0-200kHz');
% xlim([0 max(SV.Fun_Time)]);
% xlabel('Time [s]');
xlabel('Stress [MPa]');
ylabel('Energy [aJ]');
legend(num2str(round(SV.Table.MPa)));
hold off