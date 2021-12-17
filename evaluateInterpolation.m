%Evaluation script
%- Mean values: Pitch, SpectralCentroid, RMS
%- Plots: SNR

addpath('./helpers');
addpath('./wavetables');

clear, close all;

doPlots = false;            %!!!BEFORE RUNNING, PLEASE CONSIDER DOING SOME CLEANING!!!

Fs = 44100;
dur = 1;
f0 = 440; %c6 = 1046.5

f0vector = [27.5, 440, 1046.5, 4186.01];
intTypes = {'truncate', 'linear', 'cubic', 'sinc'};
arbiWave = audioread('./wavetables/braids2.wav');
waveTypes = {"square", "sawtooth", arbiWave, "noise"};


%temp increment value (for mipmaps)
N = 4;

%wavetable vector
testTable = zeros(N, 1);

%mean value vectors 
%for now: row will be wavetable, column mipmapPerOctave
pitchEstimates = zeros(N, 1);
specCentroids = zeros(N, 1);
rmsValues = zeros(N, 1);
snrValues = zeros(N, 1);

%loop for running f0 through variables:
%   1 - wavetable Type: i = 1:4
%   2 - MipmapsPerOctave: (n = 1:N)-1

for n=1:N
    for i=1:4
        %define wavetable 
        testTable = wavetable(Fs, 1, f0vector(n), 'InterpolationType', intTypes{i}, ...
            'Wavetables', waveTypes{1}, 'Oversample', 2, 'MipmapsPerOctave', n-1, 'WavetableLength',2048);
        
        %get values of test wavetables
        pitchEstimates(n,i) = mean(pitch(testTable, Fs, 'Range', [20, 5000], 'Method', 'LHS'));
        specCentroids(n,i) = mean(spectralCentroid(testTable,Fs));
        rmsValues(n,i) = rms(testTable);
        snrValues(n,i) = snr(testTable, Fs, round(Fs/2/f0vector(n)));
        
        if doPlots
            %plot
            f1 = figure, snr(testTable, Fs, round(Fs/2/f0vector(n)));
            f2 = figure('Position', [500, 300, 1000, 400]), tfPlotToPrint(testTable, Fs);
            %figure, spectrogram(testTable, 512, 64, 512, Fs, 'yaxis');
            
            %save the plots
            
%             if i == 3
%                 fileName = sprintf('./plots/Figure wavetable %s, Mipmaps%d.jpg', "arbitrary", n-1);
%                 file2Name = sprintf('./plots/tfPlot wavetable %s, Mipmaps%d.jpg', "arbitrary", n-1);
%             else
                
%                 fileName = sprintf('./plots/SNR wavetable %s, Mipmaps%d.jpg', waveTypes{i}, n-1);
%                 file2Name = sprintf('./plots/tfPlot wavetable %s, Mipmaps%d.jpg', waveTypes{i}, n-1);
            
            %fileName = sprintf('./plots/SNR oversample factor %d, mipmaps %d.jpg', i, n-1);
            %file2Name = sprintf('./plots/tfPlot oversample factor %d, mipmaps %d.jpg', i, n-1);
                
            %end
            
            %exportgraphics(f1,fileName, 'Resolution','400');
            %exportgraphics(f2,file2Name, 'Resolution','400');
            
        end
    end
end