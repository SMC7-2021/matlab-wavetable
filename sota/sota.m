clear; close all; clc;

addpath('../helpers');
[x, Fs] = audioread('./surge_sine_c5_long.aif');
% [x, Fs] = audioread('./retro_square_c7.aif');
tfPlot(x(20000:end, 1), Fs, .01);
figure; spectrogram(x(:,1), 512, 64, 512, Fs, 'yaxis');