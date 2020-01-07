function rampedSound = onsetOffsetRamp(wav,ms,Fs)

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