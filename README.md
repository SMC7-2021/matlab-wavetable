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

`'WavetableType'` (string) default `square`. Only `square` supported for now.

`'MipmapsPerOctave'` (integer) default `0`. If `0`, no mipmapping will be
applied, i.e. the same untreated wavetable will be used for all frequencies.

`'Oversample'` (integer, `>= 1`) default `1`. If `1`, no oversampling will take
place. NB, mipmaps will be calculated with respect to `Fs * Oversample`.

## Returns
`y` (vector) the output signal.