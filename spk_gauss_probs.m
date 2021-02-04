function [avg_probs,full_probs] = spk_gauss_probs(orig,new)
% spk_gauss_probs calculates the probability of new individual spikes 
% belonging to the same unit as a previously defined single unit. Usage:
%   [avg_probs,full_probs] = spk_gauss_probs(orig,new);
% where orig and new are both [m x n] matrices of m spikes across n data
% points, with orig being your cleanly defined unit, and new being the
% spikes to calculate probabilities for.
%
% This takes advantage that the distribution of voltages at each data point
% in a clean unit should be Gaussian, but there is a fallacy in thinking
% that those probabilities relate to the probability that that value came
% from that distribution. That is to say, if we had a ground truth and knew
% that a given waveform came from the tail of that neuron's distributions,
% say with a probability of 0.1, that gives us the probability of selecting
% that waveform from the unit's true waveforms, but NOT the probability
% that that waveform came from that unit. We are therefore using the
% Gaussian primarily for its shape, and scale it such that its peak value
% is equal to 1, thereby saying a waveform that went through the most
% likely voltage at each data point has a match confidence of 1 for coming
% from that unit.
%
% For more advanced match confidence calculations we recommend the use of a
% Gaussian mixture model instead, as employed when the setting prob_method
% is set to 'gmm' (default) in the function "template_match_convhull" in 
% the same repository as this.
%
% E. M. Merricks, Ph.D. 2018-07-18

% make the voltage range we calculate the probabilities over (rounded to go
% slightly above and below the max/min voltage across all spikes)
yscale = [min([new(:); orig(:)]) max([new(:); orig(:)])];
yscale(1) = 10*floor(yscale(1)/10);
yscale(2) = 10*ceil(yscale(2)/10);

% pre-allocate memory
pdfs = zeros(length(yscale(1):yscale(2)),size(orig,2));
full_probs = NaN(size(new));
probs_raw = NaN(size(new));

% for each time point, calculate the probability the value for each new
% wave would come from the distribution of original spikes
for m = 1:size(orig,2)
    % fit a distribution to the original spikes' voltages at this datapoint
    pd = fitdist(orig(:,m),'Normal');
    % calculate the probabilities for each new spike's voltage at this time
    probs_raw(:,m) = pd.pdf(new(:,m));
    % calculate the actual probability density for a normal distribution
    % around the original spikes at this datapoint (because noise should be
    % from background in a cleanly sorted unit, which should be gaussian)
    pdfs(:,m) = pd.pdf(yscale(1):yscale(2));
    % divide new probabilities by max in prob density because at that
    % voltage at that time point, it would be an exact match:
    full_probs(:,m) = probs_raw(:,m)/max(pdfs(:,m));
    % (A discrete equivalent of multiplying by sqrt(2*SD*pi))
end
% we now have the probability for each data point in each new spike, make
% the average probability across each spike:
avg_probs = mean(full_probs,2);
