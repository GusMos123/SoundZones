clc
clear
close all

% Parameters
fs = 16000;         % Sampling frequency
f_sinusoid = 1000;  % Frequency of the sinusoid
duration = 1;       % Duration of the signal in seconds

% Time vector
t = 0:1/fs:duration-(1/fs);

% Create sinusoidal signal
x = sin(2*pi*f_sinusoid*t);

% Perform FFT
nfft = 2048;  % FFT size
X = fft(x, nfft);

% Compute the frequency axis
frequencyAxis = linspace(-fs/2, fs/2, nfft);

% Shift the FFT result to center around 0 Hz
X_shifted = fftshift(X);

% Plot the magnitude spectrum
figure;
plot(frequencyAxis, abs(X_shifted));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum (Centered)');


% Perform inverse FFT
x_reconstructed = ifft(X, nfft);

% Plot original and reconstructed signals
figure;
subplot(2,1,1);
plot(t(1:200), x(1:200));
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Signal');

subplot(2,1,2);
t_reconstructed = (0:length(x_reconstructed)-1) / fs;
plot(t_reconstructed(1:200), real(x_reconstructed(1:200)));  % Plot real part due to numerical artifacts
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed Signal');
