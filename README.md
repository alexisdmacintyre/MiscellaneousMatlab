# MiscellaneousMatlab
Odds and ends that may be of use...

## calculateAsynchrony.m
MATLAB function for calculating asynchrony (euclidean distance) between vectors of datapoints (e.g., time stamps, drum hits, linguistic annotations). Plots your pairings and produces a 1d series of absolute distances between paired points, as determined by serial closest two-way matches. Optionally returns indices of unpaired values from each input vector.

![Example usage](https://user-images.githubusercontent.com/55560694/65283449-11cb6900-db2f-11e9-8bb3-2ce39d48aa18.png)

## timeStamps2Vector.m
MATLAB function that accepts a vector of time stamps to be converted to a one-hot coded vector. Optionally resamples output vector to desired frequency.

## genPinkNoise.m
MATLAB function to create pink noise.

## onsetOffsetRamp.m
MATLAB function to add onset/offset volume ramping to beginning and ending of WAV.
