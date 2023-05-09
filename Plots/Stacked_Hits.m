figure('name', 'Frequency vs Stacked hits','Position',[60,60,1400,700])
disp("Stacking hits...");
hold on
multiplier = 1;
Spreadf = cell(1,multiplier*ceil(max(PFreqList./1000)));
Spreadf2 = [];
for i = 1 : length(PFreqList)
    Spreadf{1,round(PFreqList(i)./1000*multiplier)} = ...
        [Spreadf{1,round(PFreqList(i)./1000*multiplier)};...
        i*TimeEnd/HighestIndex];
end
% Spreadf = Spreadf(1,1:(1000*multiplier)); %Cap the frequency to 1MHz
for i = 1 : length(Spreadf)
    len = length(Spreadf{1,i});
    if len ~= 0
        if len == 1
            Spreadf2(1,i) = 1;
        else
            Spreadf2(1:len,i) = Spreadf{1,i};
        end
    end
end
imagesc(Spreadf2/multiplier);
cc = colorbar;
title(cc,'Time [s]')
shading interp
mycolormap = colormap(hsv);
% mycolormap(1, :) = [1,1,1]; %White background
mycolormap(1, :) = [0, 0, 0]; %Black background
bound = MCstrboundtime*height(mycolormap)/80;
mycolormap(1:round(bound), :) = zeros(round(bound), 3); %Reduced lower values
colormap(mycolormap);
title('Frequency vs Stacked hits');
xlabel('Peak frequency [kHz]');
ylabel('Hits');