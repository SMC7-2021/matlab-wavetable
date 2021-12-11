clear; close all; clc;

addpath('../helpers');
[x, Fs] = audioread('wavetable_square_c6.wav');
% [x, Fs] = audioread('./retro_square_c7.aif');
x = x(1:44100, 1);
x = normalize(x)/2;
tfPlot(x(20000:end, 1), Fs, .01);
figure; spectrogram(x(:,1), 512, 64, 512, Fs, 'yaxis'); 