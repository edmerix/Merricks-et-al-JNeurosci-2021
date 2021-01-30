function [dodge,pow,frq,thr] = spike_denoise(spks,varargin)
% spike_denoise: find suspected non-physiological artifacts in a matrix of
% detected spikes. Uses the FFT of each detected spike and looks for
% outliers above a set threshold in a set frequency band.
%
% Inputs:   
%   1) spks (an [n x m] matrix of waveforms where n is the number of 
%   waveforms and m is the number of data points per waveform) [Required]
%   2+) name, value pairs of settings from:
%       Fs:             sampling frequency of the waveforms Default: 3e4
%       threshold:      multiple of SDs above which in the specified 
%                       frequencies are treated as outliers. Default: 5
%                       (if raw_threshold is set, that overrides this)
%       raw_threshold:  raw value in dB for the outlier cutoff in specified
%                       frequencies. Default: []
%       frequencies:    start and stop frequencies between which to look
%                       Default: [2500 5000] (Hz)
%       above:          look in any frequencies above this value. Default:
%                       2500 (Hz)
%       below:          look in any frequencies below this value. Default:
%                       500 (Hz)
%       nfft:           nfft for the FFT (8192)
%       combined:       whether or not to calculate the combination of
%                       frequency requests as one group or independently
%                       (i.e. with different SD calculations). Default:
%                       false
%
%   You can combine the different frequency input options, and the function
%   will treat them as a set if "combined" is true, or will run through
%   each in turn if not.
%   N.B. Be careful with combinations and the "above" option: frequencies
%   run all the way up to Fs/2, which can massively decrease the mean + SD
%   of the power in the full frequency set. Recommended to use the
%   "frequencies" version instead of "above" in order to set a high cutoff
%   rather than letting it run all the way up to Fs/2. "below" doesn't
%   suffer the same issue, as it just runs down to zero, so a combination
%   of "frequencies" and "below" with "combined" set to false is optimum in
%   most use-cases.
%
% Outputs:
%   1) dodge:   the indices of waveforms that were marked as outliers
%   2) pow:     the power from the FFT of all waveforms (in dB)
%   3) frq:     the frequencies that the above "pow" were calculated at
%   4) thr:     the actual threshold used across the specified frequencies
%
% E.M.Merricks, Ph.D. 2020

if nargin < 1 || isempty(spks)
    error('Need at least a matrix of spikes as input');
end
% Input options (these are the defaults, adjust them at runtime by setting
% them as name, value pairs):
settings.Fs = 3e4;
settings.threshold = 5;
settings.raw_threshold = [];
settings.frequencies = [2500 3500];
settings.above = 2500;
settings.below = 500;
settings.nfft = 8192;
settings.combined = false;

% Assign the user's requested settings in place of defaults:
allowable = fieldnames(settings);
if mod(length(varargin),2) ~= 0
    error('Inputs must be in name, value pairs');
end
for v = 1:2:length(varargin)
    if find(ismember(allowable,varargin{v}))
        settings.(varargin{v}) = varargin{v+1};
    else
        disp([9 'Not assigning ''' varargin{v} ''': not a setting of the spike_denoise function']);
    end
end

% Calculate the FFT across the whole matrix of waveforms:
z = fft(spks,settings.nfft,2);
halfn = floor(settings.nfft/2)+1;
deltaf = 1/(settings.nfft/settings.Fs);
frq = (0:(halfn-1)) * deltaf;
pow = NaN(size(spks,1),halfn);
pow(:,1) = abs(z(:,1))./settings.nfft;
pow(:,2:(halfn-1)) = abs(z(:,2:(halfn-1))) ./ (settings.nfft/2);
pow(:,halfn) = abs(z(:,halfn)) ./ settings.nfft;

% Find the requested frequency ranges, first with the bandpass version:
srchSize = [3 1]; % sneaky use of false+1/true+1 being cast to doubles...
srch = zeros(srchSize(settings.combined+1),length(frq));
if ~isempty(settings.frequencies)
    srch(1,:) = frq > settings.frequencies(1) & frq <= settings.frequencies(2);
end
% Then use logic "OR" gates to combine the "above" and "below" input 
% options with the bandpass version if using combined method, otherwise set
% them up as separate rows:
if ~isempty(settings.above)
    if settings.combined
        srch = srch | frq > settings.above;
    else
        srch(2,:) = frq > settings.above;
    end
end
if ~isempty(settings.below)
    if settings.combined
        srch = srch | frq < settings.below;
    else
        srch(3,:) = frq < settings.below;
    end
end
% If combined, they're all in one and srch will be 1 by n. Otherwise, run
% each in turn:
dodgeLogical = zeros(1,size(spks,1));
for s = 1:size(srch,1)
    % Calculate threshold, if not supplied as a raw value:
    if ~isempty(settings.raw_threshold)
        thr = settings.raw_threshold;
    else
        amp_subset = pow(:,srch(s,:) == 1);
        thr = mean(amp_subset(:)) + settings.threshold * (std(amp_subset(:)));
    end
    % find all waveforms that go above the threshold anywhere in the requested
    % ranges:
    [i,~] = ind2sub(size(pow(:,srch(s,:) == 1)),find(pow(:,srch(s,:) == 1) > thr));
    % no need to return the same waveform more than once:
    dodgeLogical(unique(i)) = 1;
end
dodge = find(dodgeLogical); % return the indices from the logical indices