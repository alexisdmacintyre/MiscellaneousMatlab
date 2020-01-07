function pinkNoise = genPinkNoise(ms,Fs,vol)

% pinkNoise = genPinkNoise(ms,Fs)
% ms = milliseconds
% Fs = sampling frequency
% vol = higher values are quieter, defaults to 50

% https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html

if nargin < 2 
    error('ms and Fs required')
elseif nargin < 3
    vol = 50;
elseif nargin > 3
    error('too many inputs')
end
 

Nx = round((ms/1000)*Fs);  % number of samples to synthesize
B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
A = [1 -2.494956002   2.017265875  -0.522189400];
nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
x = filter(B,A,v);    % Apply 1/F roll-off to PSD
pinkNoise = (x(nT60+1:end)./vol)';

end