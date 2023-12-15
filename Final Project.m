close all; clear all; clc;

% Read audio file
[audio, Fs] = audioread('IMG_3525.wav');

% Normalize audio signal
audio = audio / max(abs(audio));

% Play the original audio
% sound(audio, Fs);

%original audio snr
snrOri = snr(audio);

% Plot the noisy signal
figure(1);
subplot(2, 1, 1);
t = (0:length(audio)-1) / Fs;
plot(t, audio);
title(sprintf('Original Recording (SNR: %.2f dB)', snrOri));
xlabel('Time (s)');
ylabel('Amplitude');

% Spectral analysis of the signal
subplot(2, 1, 2);
N = length(audio);
f = (-Fs/2 : Fs/N : Fs/2 - Fs/N);
freqDomain = fftshift(abs(fft(audio, N)));
plot(f, freqDomain);
title('Spectrum of the Original Recording');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Uniform condition for filter testing
o = 8; % Filter order
wn = [600]/(Fs/2); % Frequency cutoff (adjust as needed)

%% Butterworth
% Design Butterworth Low Pass IIR Filter and apply on the sequence
[b, a] = butter(o, wn, 'low');

% See frequency response of the filter
[h, w] = freqz(b, a, 1024, Fs);
figure(2);
plot(w, 20*log10(abs(h)));
title('IIR Filter (Butterworth Lowpass)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Filter the signal
filteredAudio = filter(b, a, audio);

% Play the filtered audio
%sound(filteredAudio, Fs);

% Save the filtered audio as a separate audio file
filteredFilename = 'butterworthFiltered.wav';  % Specify the desired filename
audiowrite(filteredFilename, filteredAudio, Fs);

% filtered signal snr
snrButterworth = snr(filteredAudio);

% Plot the filtered signal
figure(3);
subplot(2, 1, 1);
plot(t, filteredAudio);
title(sprintf('Filtered Recording with Butterworth (SNR: %.2f dB)', snrButterworth));
xlabel('Time (s)');
ylabel('Amplitude');

% Spectral analysis of the filtered signal
subplot(2, 1, 2);
filteredFreqDomain = fftshift(abs(fft(filteredAudio, N)));
plot(f, filteredFreqDomain);
title('Spectrum of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% Chebysev
% Design Chebyshev Type 1 Low Pass IIR Filter and apply on the sequence
[b, a] = cheby1(o, 1, wn, 'low');

% See frequency response of the filter
[h, w] = freqz(b, a, 1024, Fs);
figure(4);
plot(w, 20*log10(abs(h)));
title('IIR Filter (Chebyshev Lowpass)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Filter the signal
filteredCheby = filter(b, a, audio);

% Play the filtered audio
sound(filteredCheby, Fs);

% Save the filtered audio as a separate audio file
filteredFilename = 'chebyshevFiltered.wav';  % Specify the desired filename
audiowrite(filteredFilename, filteredCheby, Fs);

% filtered signal snr
snrCheby = snr(filteredCheby);

% Plot the filtered signal
figure(5);
subplot(2, 1, 1);
plot(t, filteredCheby);
title(sprintf('Filtered Recording with Chebyshev (SNR: %.2f dB)', snrCheby));
xlabel('Time (s)');
ylabel('Amplitude');

% Spectral analysis of the filtered signal
subplot(2, 1, 2);
filteredFreqDomain = fftshift(abs(fft(filteredCheby, N)));
plot(f, filteredFreqDomain);
title('Spectrum of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% Blackman
nyquistFreq = Fs/2;
windowLength = o + 1; % Window length is one greater than the filter order
window = blackman(windowLength);
b = fir1(o, wn, 'low', window);

% See frequency response of the filter
[h, w] = freqz(b, 1, 1024, Fs);
figure(6);
plot(w, 20*log10(abs(h)));
title('FIR Filter (Blackman Lowpass)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Filter the signal
filteredAudio = filter(b, 1, audio);

% Play the filtered audio
% sound(filteredAudio, Fs);

% Save the filtered audio as a separate audio file
filteredFilename = 'blackmanFiltered.wav';  % Specify the desired filename
audiowrite(filteredFilename, filteredAudio, Fs);

% Filtered signal SNR
snrBlackman = snr(filteredAudio);

% Plot the filtered signal
figure(7);
subplot(2, 1, 1);
plot(t, filteredAudio);
title(sprintf('Filtered Recording with Blackman (SNR: %.2f dB)', snrBlackman));
xlabel('Time (s)');
ylabel('Amplitude');

% Spectral analysis of the filtered signal
subplot(2, 1, 2);
filteredFreqDomain = fftshift(abs(fft(filteredAudio, N)));
plot(f, filteredFreqDomain);
title('Spectrum of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% Hamming
windowLength = o + 1; % Window length is one greater than the filter order
window = hamming(windowLength);
b = fir1(o, wn, 'low', window);

% See frequency response of the filter
[h, w] = freqz(b, 1, 1024, Fs);
figure(2);
plot(w, 20*log10(abs(h)));
title('FIR Filter (Hamming Lowpass)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Filter the signal
filteredAudio = filter(b, 1, audio);

% Play the filtered audio
% sound(filteredAudio, Fs);

% Save the filtered audio as a separate audio file
filteredFilename = 'hammingWindowFiltered.wav';  % Specify the desired filename
audiowrite(filteredFilename, filteredAudio, Fs);

% Filtered signal SNR
snrHamming = snr(filteredAudio);

% Plot the filtered signal
figure(3);
subplot(2, 1, 1);
plot(t, filteredAudio);
title(sprintf('Filtered Recording with Hamming (SNR: %.2f dB)', snrHamming));
xlabel('Time (s)');
ylabel('Amplitude');

% Spectral analysis of the filtered signal
subplot(2, 1, 2);
filteredFreqDomain = fftshift(abs(fft(filteredAudio, N)));
plot(f, filteredFreqDomain);
title('Spectrum of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');