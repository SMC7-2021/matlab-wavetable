clear; close all;

addpath('../helpers');
%[x, Fs] = audioread('wavetable_square_a0.wav');
%[x, Fs] = audioread('./surge_square_C8_1.aif');
%[x, Fs] = audioread('Wavetable_A0.wav');
%[x, Fs] = audioread('retro_square_c7.aif');
%[x, Fs] = audioread('retro_square_A4.aif');
%[x, Fs] = audioread('vital_square_A4.aif');

x = x(1:44100, 1);
x = normalize(x)/2;

a0 = 27.5; a4 = 440; c6 = 1046; c8 = 4186;

%this value is used to find the SNR 
f0 = c8;

f1 = figure('Position', [500, 300, 1000, 400]); tfPlotToPrint(x(20000:end, 1), Fs, .01);
figure; spectrogram(x(:,1), 512, 64, 512, Fs, 'yaxis'); 

pitchEstimate = mean(pitch(x, Fs, 'Range', [20, 5000], 'Method', 'LHS'));
specCentroid = mean(spectralCentroid(x,Fs));
rmsValues = rms(x);
snrValues = snr(x, Fs, round(Fs/2/f0));

%fileN = sprintf('./SotaPlots/tfPlot abletonArbitaryC6 SOTA.jpg');
%fileN2 = sprintf()
%exportgraphics(f1,fileN, 'Resolution','400');
sound(x,Fs);

%NO! if you're still in doubt about applying f0- not for sota
%TODO: for a fourier square, perhaps...
%f0 = [440, 1046.5, 4086];
%for f=1:length(f0)
    
    