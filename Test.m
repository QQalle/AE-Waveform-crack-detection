


s = Signal;                                         % Signal Vector
Ts = 1/Fs;                                          % Sampling Interval
Fs = 1/Ts;                                          % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
L = numel(s);                                       % Signal Length (Number Of Samples)
FTs = fft(s - mean(s))/L;                           % Fourier Transform
% Fv = linspace(0, 1, fix(L/2)+1)*Fn;                 % Frequency Vector
% Iv = 1:numel(Fv);      % Index Vector

cFFT = FFTMat(:, 500);
N = length(cFFT);
Fv = linspace(0, 1, length(ccFFT));
Iv = 1:numel(Fv);

centriod = spectralCentroid(abs(FTs(Iv))*2 , Fv);

figure
plot(fVals, abs(FFTf)*2)
grid
xlabel('Frequency')
ylabel('Amplitude')
% refl = xline(centriod);
% refl.Color = 'r';


