function rampedSound = onsetOffsetRamp(wav,ms,Fs)

% Adds onset and offset volume ramping to WAV
% wav = WAV file
% ms = desired length of ramp (e.g., 20 ms)
% Fs = sample rate

if nargin < 3
    error('Needs wav, ms, and Fs input arguments')
end

if isrow(wav)
    wav  = wav';
end

ramp = round(ms/1000*Fs);
rampVector = linspace(0,1,ramp)';
wav(1:ramp) = wav(1:ramp).*rampVector;
wav(end-ramp+1:end) = wav(end-ramp+1:end).*rampVector(end:-1:1);
rampedSound = wav;
end
