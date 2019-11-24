function vectorOut = timeStamps2Vector(x,time,fs_in,fs_out)

% Accepts a vector of timestamps (e.g., drum hits, vowel onsets) and returns
% a one-hot coded vector, to be used, for example, as a feature in ML
% applications. Specify optional 4th argument fs_out to resample the
% timestamps to desired scale.

% Usage:

% x = 1xN vector of timestamps (e.g., [400 804 860 978 1001])
% time = desired temporal space of output vector in seconds (e.g., 5.5)
% fs_in = sampling rate of timestamps in Hz (e.g., 1000)
% fs_out = optional (new) output sampling rate for one-hot coded vector

switch nargin
    case 3
        resampleYesNo = 0;
        fs_out = fs_in;
    case 4
        resampleYesNo = 1;
end


if resampleYesNo == 1
    x = round(x*(fs_out/fs_in));
end

if x(1) == 0
    x(1) = 1;
end

vectorOut = zeros(1,(round(fs_out*time)));
vectorOut(x) = 1;

end