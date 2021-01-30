function [matched, raw] = gmm_match(original,new,varargin)
% [matched, raw] = gmm_match(original, new, %name,value pairs of settings%)
% 
% N.B. 2021-01-30: this function has been superseded to overhaul how
% waveform match probabilities are calculated (using a combination of
% posterior probabilities and original Gaussian fits, rather than expanding
% the Gaussians when using 'joint' = false.
%
% Original help:
%
% Template matches new waveforms to previously defined single units, using
% PCA and fitting Gaussian mixture models.
%
% Designed for use with "NeuroClass" objects: a MultipleUnits object
% containing SingleUnit objects (see http://github.com/edmerix/NeuroClass)
% but the structure of those classes can be mimicked instead, along with
% some minor edits to this function.
%
% Inputs:
%   1:  the spike sorted MultipleUnits class that new waveforms should be 
%       matched to.
%
%   2:  new data to match, as an n-by-1 struct as follows, for n channels
%           new.waveforms:  the p-by-q matrix of p spikes to match (q
%                           datapoints per spike)
%           new.times:      the 1-by-p array of spike times
%           AND EITHER:
%           new.channel:    what channel these spikes were from (to match
%                           the original data's channel assignments).
%           OR:
%           new.electrode:  what electrode label these spikes were from
%
%   3:  name, value pairs of settings. Options are:
%           threshold:      probability of match above which a waveform
%                           will be matched to that unit [0.5]
%           total_pc:       number of PCs to use must explain this
%                           percentage of variance in original spikes [95]
%           GMM_expansion:  multiple to expand fitted Gaussian by to
%                           calculate match probabilities (only used if 
%                           settings.joint == false) [5]
%           joint:          boolean, whether or not to calculate GMM
%                           distribution across all channel units at once
%                           (true GMM) or one at a time (basically just
%                           fitting n-dimensional Gaussians to each cluster
%                           in turn) [false] (currently working on this)
% Outputs:
%   1:  a new MultipleUnits object containing the template 
%       matched data, where spikes above settings.threshold are assigned to
%       their most probable unit.
%   2:  The probability of matches for all spikes to all units on their
%       channel.
%
% E. M. Merricks, Ph.D. 2020-02-10

% check inputs are correct format:
if nargin < 2 || isempty(original) || isempty(new)
    error('Need at least 2 inputs: see help for more info');
end
if ~isfield(new,'waveforms') || ~isfield(new,'times') || ~isfield(new,'channel')
    error('New data input needs to have "waveforms", "times", and "channel" field (for each channel)')
end

% default settings:
settings.threshold = 0.5;   % probability in GMM above which a spike is considered "matched"
settings.total_pc = 95;     % number of PCs to use must explain this percentage of variance
settings.GMM_expansion = 5; % multiple to expand fitted Gaussian by to calculate probabilities
settings.joint = false;     % calculate the GMM on all channel units at once, or one-by-one?
% assign user input (name, value pair) settings:
allowable = fieldnames(settings);
if mod(length(varargin),2) ~= 0
    error('Extra inputs must be in name, value pairs');
end
for v = 1:2:length(varargin)
    if find(ismember(allowable,varargin{v}))
        settings.(varargin{v}) = varargin{v+1};
    else
        disp([9 'Not assigning ''' varargin{v} ''': not an input of gmm_match function']);
    end
end

% set up outputs:
raw = struct();
matched = MultipleUnits('patient',original.patient,'seizure',original.seizure);
matched.info = 'Automatically template-matched data using GMMs';
matched.epoch = original.epoch;
if min([new.times]) < matched.epoch(1)
    matched.epoch(1) = min([new.times]);
end
if max([new.times]) > matched.epoch(2)
    matched.epoch(2) = max([new.times]);
end
matched.extra = struct();
matched.extra.GM_settings = settings;

% work our way through all new matrices of waveforms:
for n = 1:length(new)
    spks = new(n).waveforms;
    spkt = new(n).times;
    if isfield(new(n),'channel')
        isElec = false;
        chan = new(n).channel;
        disp([9 'Working on input ' num2str(n) ' (channel ' num2str(chan) ')'])
        units = original.channel_units(chan);
    elseif isfield(new(n),'electrode')
        isElec = true;
        elec = new(n).electrode;
        disp([9 'Working on input ' num2str(n) ' (electrode ' elec ')'])
        units = original.units(strcmpi({original.units.electrodelabel}, elec));
    else
        error(['Missing electrode or channel info in item ' num2str(n)])
    end
    
    if isempty(units)
        disp([9 9 'No original units on channel, skipping'])
    else
        wvs = cell2mat({units.waveforms}');
        amts = cellfun(@length,{units.times});

        assigns = NaN(1,sum(amts));
        for a = 1:length(amts)
            reached = sum(amts(1:a-1));
            if isempty(reached)
                reached = 0;
            end
            assigns(reached+1:reached+amts(a)) = a;
        end

        % calculate the PCA on the original units:
        warning('off','stats:pca:ColRankDefX')
        [coef,pc,~,~,expl] = pca(wvs,'Centered',false);
        warning('on','stats:pca:ColRankDefX')
        amt_expl = cumsum(expl);
        inds = find(amt_expl >= settings.total_pc);
        nPC = inds(1);
        % calculate the PC scores for the new waveforms on the original unit
        % PCA:
        pc_post = spks * coef;

        unq = unique(assigns); % should always be 1:length(amts) but just in case
        match_probs = NaN(size(spks,1),length(unq));
        
        if settings.joint
            scale_probs = NaN(size(spks,1),length(unq));
            gmd = fitgmdist(pc(:,1:nPC),length(units));
            probs = posterior(gmd,pc_post(:,1:nPC));
            for u = 1:length(unq)
                tempGMD = gmdistribution(gmd.mu(u,:),gmd.Sigma(:,:,u));
                scale_probs(:,u) = tempGMD.pdf(pc_post(:,1:nPC))./tempGMD.pdf(tempGMD.mu);
            end
            % this won't work if only one scale_probs(1,:) is > 0:
            % (scaling will be 100% the original probs value in that case)
            match_probs = probs .* (scale_probs./sum(scale_probs,2));
            % OR, instead of the for loop and above line:
            % match_probs = probs .* gmd.ComponentProportion; % ...?
            % and then choose which compmonent on probs, not on match_probs
        else
            gmd = cell(1,length(unq));

            % for each unit, fit the GMM on requested number of PCs, and calculate
            % the probabilities for each new spike for that unit:
            for u = 1:length(unq)
                these_pc = pc(assigns == unq(u),1:nPC);
                % fitgmdist is overkill when only using 1 cluster, but it's speedy:
                gmd{u} = fitgmdist(these_pc,1);
                % expand the fit by the requested amount:
                big_gmd = gmdistribution(gmd{u}.mu, settings.GMM_expansion * gmd{u}.Sigma);
                % calculate the match probabilities, scaling the Gaussian such that
                % a waveform at its exact mean has a match probability of 1:
                match_probs(:,u) = big_gmd.pdf(pc_post(:,1:nPC))./big_gmd.pdf(big_gmd.mu);
            end
        end

        % store the raw GMM fits and probabilities (these should be consulted
        % rather than taking output at face value... probability thresholds
        % should probably be different for each unit based on the original
        % unit's density in PC space. Need to work out a mathematical basis for
        % defining this threshold):
        if isElec
            raw(n).electrode = elec;
        else
            raw(n).channel = chan;
        end
        raw(n).GMM_fits = gmd;
        if ~settings.joint
            raw(n).GMM_expansion = settings.GMM_expansion;
        end
        raw(n).UIDs = [units.UID];
        raw(n).probabilities = match_probs;

        % get maximal match probabilities for each spike, and which unit it was
        % assigned to:
        [mx,w] = max(match_probs,[],2);
        if settings.joint
            [~,w] = max(probs,[],2);
        end
        inds = find(mx > settings.threshold);
        new_unq = unique(w(inds));
        % for each spike that surpasses threshold in its max probability,
        % assign that spike to that unit:
        temp_unit = SingleUnit;
        for nu = 1:length(new_unq)
            temp_unit(nu) = SingleUnit('patient',original.patient,'seizure',original.seizure);
            matches = inds(w(inds) == new_unq(nu));
            temp_unit(nu).waveforms = spks(matches,:);
            temp_unit(nu).times = spkt(matches);
            if isElec
                temp_unit(nu).electrodelabel = elec;
            else
                temp_unit(nu).channel = chan;
            end

            temp_unit(nu).type = units(new_unq(nu)).type;
            temp_unit(nu).UID = units(new_unq(nu)).UID;
            if ~settings.joint
                % will need to match the results to the original units if
                % doing joint method, so don't know this yet:
                temp_unit(nu).extra.mean_original_waveform = mean(units(new_unq(nu)).waveforms);
            end
            temp_unit(nu).extra.probabilities = mx(matches);

            matched.add_unit(temp_unit(nu));
        end
    end
end