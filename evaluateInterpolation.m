%Evaluation script on interpolation methods
%- Mean values: Pitch, SpectralCentroid, RMS
%- Plots: SNR

clear, close all;

Fs = 44100;
dur = 1;
f0 = 440; %c6 = 1046.5
intTypes = {'truncate', 'linear', 'cubic', 'sinc'};
%waveTypes = {"square", "sawtooth", "sine", "noise"};

%temp increment value (for mipmaps)
N = 5;

%wavetable vector
testTable = zeros(N, 1);

%mean value vectors 
%for now: row will be interpolation, column mipmapPerOctave
pitchEstimates = zeros(N, 1);
specCentroids = zeros(N, 1);
rmsValues = zeros(N, 1);
snrValues = zeros(N, 1);

%loop for running f0 through variables:
%   1 - InterpolationType: i = 1:4
%   2 - MipmapsPerOctave: (n = 1:N)-1 (number of mipmaps 0:4)

for n=1:N
    for i=1:4
        %define wavetable with different interpolation methods
        testTable = wavetable(Fs, 1, f0, 'InterpolationType', intTypes{i}, ...
            "square", 'Oversampling', 2, 'MipmapsPerOctave', n-1);
        
        %get values of test wavetable
        pitchEstimates(n,i) = mean(pitch(testTable, Fs, 'Range', [50, 5000], 'Method', 'LHS')); 
        specCentroids(n,i) = mean(spectralCentroid(testTable,Fs));
        rmsValues(n,i) = rms(testTable);
        snrValues(n,i) = snr(testTable, Fs, round(Fs/2/f0));
        
        %plot
        f1 = figure, snr(testTable, Fs, round(Fs/2/f0));
        %figure, spectrogram(testTable, 512, 64, 512, Fs, 'yaxis');
        
        %save the plots 
        fileName = sprintf('./plots/Figure intType%d, Mipmaps%d.jpg', i, n-1);
        exportgraphics(f1,fileName)
  
    end
end