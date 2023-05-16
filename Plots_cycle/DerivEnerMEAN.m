figure('name', 'Energy-time derivative','Position',[60,60,1400,700])
disp("calculating Energy-time derivative...")
hold on
% avg = zeros(1,length(SV.DerivEner{1}));
% for k = 1 : exp
%     for i = 1 : length(SV.DerivEner{k})
% %         if i <= length(SV.DerivEner{1})
%         avg(i) = avg(i) + SV.DerivEner{k}(i);
% %         end
%     end
% end
% 
% avg2 = zeros(1,length(SV.Fun_TensileStress{1}));
% for k = 1 : exp
%     for i = 1 : length(SV.Fun_TensileStress{k})
%         avg2(i) = avg2(i) + SV.Fun_TensileStress{k}(i);
%     end
% end
% 
% flipped = flip(SV.Table{:,1});
% e = [2,4,6,8,9,10,11,12,13,14];
% k = [14,12,10,8,6,5,4,3,2,1];
% for i = 1 : length(e)
%     if i == 1
%         firstvalue = 0;
%     else
%         firstvalue = round(flipped(e(i-1))*2.4857);
%     end
%     nextvalue = round(flipped(e(i))*2.4857);
%     if i == 10
%         avg(firstvalue+1:end) = avg(firstvalue+1:end)./k(i);
%         avg2(firstvalue+1:end) = avg2(firstvalue+1:end)./k(i);
%     else
%         avg(firstvalue+1:nextvalue) = avg(firstvalue+1:nextvalue)./k(i);
%         avg2(firstvalue+1:nextvalue) = avg2(firstvalue+1:nextvalue)./k(i);
%     end
% end

% avg(1:(round(200*2.4857))) = avg(1:(round(200*2.4857)))./14;
% avg((round(200*2.4857)):(round(250*2.4857))) = avg((round(200*2.4857)):(round(250*2.4857)))./12;
% avg((round(250*2.4857)):(round(300*2.4857))) = avg((round(250*2.4857)):(round(300*2.4857)))./10;
% avg((round(300*2.4857)):(round(350*2.4857))) = avg((round(300*2.4857)):(round(350*2.4857)))./8;
% avg((round(350*2.4857)):(round(400*2.4857))) = avg((round(350*2.4857)):(round(400*2.4857)))./6;
% avg((round(400*2.4857)):(round(500*2.4857))) = avg((round(400*2.4857)):(round(500*2.4857)))./5;
% avg((round(500*2.4857)):(round(550*2.4857))) = avg((round(500*2.4857)):(round(550*2.4857)))./4;
% avg((round(550*2.4857)):(round(600*2.4857))) = avg((round(550*2.4857)):(round(600*2.4857)))./3;
% avg((round(600*2.4857)):(round(650*2.4857))) = avg((round(600*2.4857)):(round(650*2.4857)))./2;
% % avg((round(650*2.4857)):(round(700*2.4857))) = avg((round(650*2.4857)):(round(700*2.4857)))./1;
% avg = smooth(avg,10);

for k = 1 : exp
    PullStopInd = find(SV.Fun_TensileStress{k} >= SV.Table{k,1}, 1);
    plot(SV.Fun_TensileStress{k}, SV.DerivEner{k},Markers(k),...
        'MarkerIndices',PullStopInd,'MarkerSize',10);
%     disp(k);
%     disp(append('Peaks: ',num2str(EnerPeaksNo)));
%     disp(append('Peaks divided: ',num2str(EnerPeaksNo2)));
end
% stress = SV.Fun_TensileStress{1};
% stress = sortrows(stress);
% plot(stress, avg,'-','LineWidth',1,'Color','k');

title('Energy Derivative vs Stress');
legend([Legendtext;"MEAN"],'location','east outside');
xlabel('Stress [MPa]');
ylabel('Energy derivative');
set(gca,'FontSize',14)
hold off