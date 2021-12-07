# Wavetable synthesis

## Syntax

```
y = wavetable(Fs, duration, F0, 'OptionalParam', 'optionalValue', ... );
```

## Parameters

`Fs` (integer) sampling rate.

`duration` (double) output duration in seconds.

`F0` (double or vector) if double, wavetable will produce output of fixed 
frequency. If vector, wavetable will use the first and second elements in the 
vector to create a linear sweep from `F0(1)` to `F0(2)`.

## Optional parameters
Specified as key/value pairs.

`'OutputAmplitude'` (double) default `1`.

`'InterpolationType'` (string) one of `truncate` (default), `linear`, `cubic`, 
or `sinc` to be applied to the wavetable when handling fractional phase values.

`'Wavetables'` (string or vector) default `"square"`. Either a single string
wrapped in *double quotes*, one of `"square"`, `"sine"`, `"sawtooth"`, `"noise"`
or `"sineBroken"`, or an array of those. Alternatively, specify a vector of 
audio samples; this will be split into chunks of 256 samples and the chunks
resampled to 2048-sample wavetables.

`'MipmapsPerOctave'` (integer) default `0`. If `0`, no mipmapping will be
applied, i.e. the same untreated wavetable will be used for all frequencies.

`'Oversample'` (integer, `>= 1`) default `1`. If `1`, no oversampling will take
place. NB, mipmaps will be calculated with respect to `Fs * Oversample`.

## Returns
`y` (vector) the output signal.