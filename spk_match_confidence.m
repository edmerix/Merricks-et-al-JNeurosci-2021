function [avg_probs,full_probs] = spk_match_confidence(orig,new)
% spk_match_confidence calculates the probability of new individual spikes 
% belonging to the same unit as a previously defined single unit. Usage:
%   [avg_probs,full_probs] = spk_match_confidence(orig,new);
% where orig and new are both [m x n] matrices of m spikes across n data
% points, with orig being your cleanly defined unit, and new being the
% spikes to calculate probabilities for.
%
% This works by fitting a Gaussian to the original unit's voltages at each
% data point (a clean unit should follow a normal distribution in clean
% recordings), and then for each value in a new waveform, it subtracts the
% probability of finding a value further from the mean in the original
% distribution than that from 1. Then these confidences are averaged across
% all data points, to provide a single match confidence for that waveform.
%
% Outputs:
%   1) The single value match confidence for each new waveform
%   2) The raw match confidences at each data point for each new waveform
%
% E. M. Merricks, Ph.D. 2021-02-03 (building on spk_gauss_probs,
% 2018-07-18)

% pre-allocate memory
full_probs = NaN(size(new));

% for each time point, calculate the probability the value for each new
% wave would come from the distribution of original spikes
for m = 1:size(orig,2)
    % fit a distribution to the original spikes' voltages at this datapoint
    pd = fitdist(orig(:,m),'Normal');
    % calculate distance from the mean in the fit:
    rngVals = abs(new(:,m) - pd.mu);
    % calculate the probability of finding a "true" value further than this
    % distance from the mean, given the fitted distribution's SD:
    full_probs(:,m) = 1 - (normcdf(rngVals,0,pd.sigma) - normcdf(-rngVals,0,pd.sigma));
end
% we now have the probability for each data point in each new spike, make
% the average probability across each spike:
avg_probs = mean(full_probs,2);
