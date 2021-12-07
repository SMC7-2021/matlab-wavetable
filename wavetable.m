function y = wavetable(Fs, duration, F0, varargin)
    Lt = 2^11;
    outAmp = 1;
    oversample = 1;
    interpolationType = 'truncate';
    wavetableType = 'square';
    
    % Mipmap perameters
    numOctaves = 9;
    mipmapsPerOctave = 0;
    
    for i = 1:nargin-3
        switch varargin{i}
            case 'Oversample'
                if floor(varargin{i + 1}) ~= varargin{i + 1} || varargin{i + 1} < 1
                    error('Oversampling factor must be a positive integer.');
                end
                oversample = varargin{i + 1};
            case 'InterpolationType'
                if ~ismember(varargin{i + 1}, ['truncate', 'linear', 'cubic', 'sinc'])
                    error("InterpolationType must be one of 'zeroth', 'linear', 'cubic', or 'sinc'.");
                end
                interpolationType = varargin{i + 1};
            case 'WavetableType'
                wavetableType = varargin{i + 1};
            case 'OutputAmplitude'
                if ~isnumeric(varargin{i + 1})
                    error('OutputAmplitude must be a number.');
                end
                outAmp = varargin{i + 1};
            case 'MipmapsPerOctave'
                if floor(varargin{i + 1}) ~= varargin{i + 1} || varargin{i + 1} < 0
                    error('MipmapsPerOctave must be a positive integer or zero.');
                end
                mipmapsPerOctave = varargin{i + 1};
        end
    end
    
    % Adjust sampling rate according to oversampling factor.
    Fs = Fs * oversample;
    
    % Make F0 a vector.
    switch length(F0)
        case 1
            F0 = linspace(F0, F0, Fs*duration)';
        otherwise
            F0 = linspace(F0(1), F0(2), Fs*duration)';
    end
    
    % Create an array of source wavetables.
    sourceWavetables = {
        % A sampled wavetable.
%         resample(x1, Lt, length(x1));
        % A sine wavetable.
%         sin(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');
        % A discontinuous sine wavetable.
%         sin(linspace(0, 2.2*pi, Lt)');
        % A square wavetable.
        square(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');
        % Another sampled wavetable.
%         resample(x2, Lt, length(x2));
        % A sawtooth wavetable.
%         sawtooth(linspace(0, (2*pi) - (2*pi)/Lt, Lt)');
        % A noise wavetable
%         (rand(Lt, 1) * 2) - 1;
    };
    
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
                x(:, f) = computeMipmap(wt, Fs, basisF0s(f), mipmapsPerOctave);
            end
        else
            x = sourceWavetables{i};
        end

        wavetables{i} = x;
    end

    % Placeholder for the wavetable morph samples.
    wt = zeros(2, 1);
    % Initialize the wavetable sample index.
    phase = 1;
    % Create output placeholder.
    y = zeros(Fs * duration, 1);
    
    % Generate output.
    for n=1:Fs*duration
        % Calculate how far through the output we are.
        outputDelta = n / (Fs * duration);
        % Multiply this by the number of regions between wavetables (e.g. 5 
        % wavetables, 4 regions) to work out how far we are between the current 
        % wavetable pair.
        wtDelta = outputDelta * (length(wavetables) - 1);
        % The integer part is the index of the previous wavetable.
        wtIndex = floor(wtDelta) + 1;
        % The fractional part is the morph amount between the pair of wavetables.
        morph = rem(wtDelta, 1);
        % Calculate the index of the next wavetable.
        nextWtIndex = wtIndex + 1;
        % At the very end of output, wtDelta == length(wavetables) so nextWtIndex 
        % will exceed the length of the wavetable array; so check for this.
        if nextWtIndex > length(wavetables)
            nextWtIndex = length(wavetables) - 1;
        end
        if nextWtIndex == 0
            nextWtIndex = 1;
        end

        % Determine, based on the current output frequency, which mipmap to use.
        % Default to mipmap 1 (for F0 below Fs/Lt)
        mipmap = 1;
        for f=length(basisF0s):-1:1
            if F0(n) > basisF0s(f)
                mipmap = f;
                break
            end
        end

        % Calculate the transitional wavetable samples.
        switch interpolationType
            case 'truncate'
                wt(1) = interpolateZeroth(wavetables{wtIndex}(:, mipmap), phase);
                wt(2) = interpolateZeroth(wavetables{nextWtIndex}(:, mipmap), phase);
            case 'linear'
                wt(1) = interpolateLinear(wavetables{wtIndex}(:, mipmap), phase);
                wt(2) = interpolateLinear(wavetables{nextWtIndex}(:, mipmap), phase);
            case 'sinc'
                wt(1) = interpolateSinc(wavetables{wtIndex}(:, mipmap), phase);
                wt(2) = interpolateSinc(wavetables{nextWtIndex}(:, mipmap), phase);
            case 'cubic'
                wt(1) = interpolateCubic(wavetables{wtIndex}(:, mipmap), phase);
                wt(2) = interpolateCubic(wavetables{nextWtIndex}(:, mipmap), phase);
        end

        % Calculate the output sample.
        y(n) = outAmp * (...
            (1 - morph) * wt(1) + ...
            morph * wt(2) ...
        );

        % Calculate the number of samples per period of the wavetable to produce the
        % current frequency.
        sampsPerPeriod = Fs / F0(n);
        phaseIncrement = Lt / sampsPerPeriod;

        % Upadate the wavetable sample index for the next iteration.
        phase = mod(phase + phaseIncrement, Lt);
    end
    
    if oversample > 1
        y = decimate(y, oversample);
    end
end

