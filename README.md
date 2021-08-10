# MiscellaneousMatlab
Odds and ends that may be of use...

## breathTimes.m
MATLAB function for detecting the beginning and ending of inhalations during respiration, optimised for speech breathing.

![Example usage](https://user-images.githubusercontent.com/55560694/95772577-ddd06200-0cb4-11eb-8176-03d63cf7f005.png)


## calculateAsynchrony.m
MATLAB function for calculating asynchrony (euclidean distance) between vectors of datapoints (e.g., time stamps, drum hits, linguistic annotations). Plots your pairings and produces a series of absolute distances between paired points, as determined by serial closest two-way matches. Optionally returns indices of unpaired values from each input vector.

![Example usage](https://user-images.githubusercontent.com/55560694/95774839-1e31df00-0cb9-11eb-9233-a430d632ce5f.png)

## timeStamps2Vector.m
MATLAB function that accepts a vector of time stamps to be converted to a one-hot coded vector. Optionally resamples output vector to desired frequency.

## genPinkNoise.m
MATLAB function to create pink noise.

## onsetOffsetRamp.m
MATLAB function to add onset/offset volume ramping to beginning and ending of WAV.
