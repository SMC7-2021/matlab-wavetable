close all;
Fs = 44100;
duration = 0.5;
Lt = 2^11;
outAmp = 1;
oversample = 1;
F0 = 440;
interpolationType = 'truncate';
wavetableType = "square";

    fontName = 'Times';
    fontSize = 12;
    

% Mipmap perameters
numOctaves = 9;
mipmapsPerOctave = 1;
% Adjust sampling rate according to oversampling factor.
Fs = Fs * oversample;

    % Make F0 a vector.
    switch length(F0)
        case 1
            F0 = linspace(F0, F0, Fs*duration)';
        otherwise
            F0 = linspace(F0(1), F0(2), Fs*duration)';
    end
    
    sourceWavetables = cell(length(wavetableType), 1);
    

sourceWavetables{i} = square(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');

    % Create a placeholder for the wavetable mipmaps.
    wavetables = cell(1, length(sourceWavetables));

        % Create a vector of frequencies to serve as the basis for the wavetable 
    % mipmaps.
    % Start from the 'fundamental' frequency of a wavetable of length Lt at
    % sampling rate Fs ((44100*oversample)/2048 = oversample*21.5332…) and 
    % proceed according to mipmapsPerOctave.
    mipmapFreqRatio = 2^((12/mipmapsPerOctave)/12);
    basisF0s = (Fs/Lt) * mipmapFreqRatio.^(0:numOctaves*mipmapsPerOctave);

 % Set up wavetable mipmaps.
    for i=1:length(sourceWavetables)
        if mipmapsPerOctave > 0
            wt = sourceWavetables{i};
            x = zeros(Lt, length(basisF0s));

            % For each fundamental frequency, compute a mipmap of the wavetable 
            % bandlimited to Nyqvist.
            for f = 1:length(basisF0s)
                x(:, f) = computeMipmap(wt, Fs, oversample, basisF0s(f), mipmapsPerOctave);
            end
        else
            x = sourceWavetables{i};
        end

        wavetables{i} = x;
    end
%%
    plottable = nan;
t=linspace(1,1,length(x))';
   figure 
    for nn=1:width(x)
   
        plottable = [plottable x(:,nn)'];

    end
    plot(plottable)
    hold off
        xlim([0 2048*(width(x))]), ...
        set(gca,'fontsize',fontSize,'fontname',fontName), ...
        xlabel('Octave'), ...
        zlabel('Numerical value'), ...
        grid on;    
        set(gcf, 'color', 'none');  
        set(gca,'XTick',[1 2048 2048*2 2048*3 2048*4 2048*5 2048*6 2048*7 2048*8 2048*9]), ...
        set(gca,'XTickLabel',[{'1st'}, '2nd' '3rd' '4th' '5th' '6th' '7th' '8th' '9th' '10th']), ...
export_fig mipmapillu.png -m4
        