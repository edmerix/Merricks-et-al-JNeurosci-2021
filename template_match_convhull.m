function matched = template_match_convhull(original,new,varargin)
% matched = template_match_convhull(original, new, %name, value settings%)
%
% Template matches new waveforms to previously defined single units, using
% PCA and convex hulls. Calculates the probability that each waveform was a
% match to the original unit using either the Gaussian distribution of the
% original unit's waveforms at each time point, or building a
% Gaussian mixture model from the original unit's PCA and getting the
% posterior probabilities of the newly assigned data points, depending on
% supplied settings. (Defaults to Gaussian mixture model).
%
% Designed for use with "NeuroClass" objects: a MultipleUnits object
% containing SingleUnit objects (see http://github.com/edmerix/NeuroClass)
% but the structure of those classes can be mimicked instead, along with
% some minor edits to this function.
%
% Inputs:
%   1:  the spike sorted MultipleUnits class that new waveforms should be
%       matched to.
%   2:  new data to match, as an n-by-1 struct as follows, for n channels
%           new.waveforms:  the p-by-q matrix of p spikes to match (q
%                           datapoints per spike)
%           new.times:      the p-by-1 array of spike times
%           new.channel:    what channel these spikes were from (to match
%                           the original data's channel assignments).
%   3:  name, value pairs of settings. Options are:
%           prob_method:    Method by which to assess spike match
%                           probability to the original unit: ['gmm']
%                               'gmm': assess probabilities based on
%                                      Gaussian mixture model in PC space
%                               'voltage': assess probabilities based on
%                                      Gaussian distributions of voltages
%                                      at each data point in original unit.
%           threshold:      probability of match above which a waveform
%                           will be matched to that unit (set to zero to
%                           include all convex hull matches but store the
%                           probabilities in the output) [0]
%           total_pc:       number of PCs to use must explain this
%                           percentage of variance in original spikes [95]
%           max_pc:         maximum total PCs to use (help speed up
%                           processing at expense of variance explained)
%                           [5]
%           centered:       whether or not to center the data for PCA. When
%                           set to false, the new waveforms' scores can be
%                           calculated simply by multiplying them by the
%                           old waveform's PC coefficients. If true, then
%                           the new waveforms should be centered manually
%                           before calculating their scores (done
%                           automatically by the function)
%                           [true]
%           joint:          boolean, whether or not to calculate GMM
%                           distribution across all channel units at once
%                           (true GMM) or one at a time (basically just
%                           fitting n-dimensional Gaussians to each cluster
%                           in turn) [true]
% Outputs:
%   1:  a new MultipleUnits object containing the template
%       matched data, where spikes above settings.threshold are assigned to
%       their most probable unit.
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
settings.prob_method = 'gmm';
settings.threshold = 0;     % probability above which a spike is considered "matched"
settings.total_pc = 95;     % number of PCs to use must explain this percentage of variance
settings.max_pc = 5;        % maximum total PCs to use (help speed up processing at expense of variance explained)
settings.centered = true;   % automatically center the PCA or not
settings.joint = true;      % calculate the GMM on all channel units at once, or one-by-one?

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
matched = MultipleUnits('patient',original.patient,'seizure',original.seizure);
matched.info = 'Automatically template-matched data using convex hulls';
matched.epoch = original.epoch;
if min([new.times]) < matched.epoch(1)
    matched.epoch(1) = min([new.times]);
end
if max([new.times]) > matched.epoch(2)
    matched.epoch(2) = max([new.times]);
end
matched.extra = struct();
matched.extra.convhull_settings = settings;


% work our way through all new matrices of waveforms:
for n = 1:length(new)
    spks = new(n).waveforms;
    spkt = new(n).times;
    chan = new(n).channel;
    
    if isnumeric(chan)
        chanstr = num2str(chan);
    else
        chanstr = chan;
    end
    disp([9 'Working on input ' num2str(n) ' (channel ' chanstr ')'])
    
    units = original.channel_units(chan);
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
        [coef,pc,~,~,expl] = pca(wvs,'Centered',settings.centered);
        warning('on','stats:pca:ColRankDefX')
        amt_expl = cumsum(expl);
        inds = find(amt_expl >= settings.total_pc);
        nPC = inds(1);
        nPC = max(nPC,2); % at least 2 PCs should be used
        nPC = min(nPC,settings.max_pc); % no more than max_pc
        disp([9 9 'Using ' num2str(nPC) ' dimensions'])
        % calculate the PC scores for the new waveforms on the original unit
        if ~settings.centered
            pc_post = spks * coef;
        else
            post_waves = spks - mean(spks,1);
            pc_post = post_waves * coef;
        end
        
        unq = unique(assigns); % should always be 1:length(amts) but just in case
        
        clusterGM = cell(1,length(unq));
        mnClusterPC = NaN(length(unq),nPC);
        
        if settings.joint && strcmp(settings.prob_method,'gmm')
            gmd = fitgmdist(pc(:,1:nPC),length(units));
            id_probs = gmd.posterior(pc_post(:,1:nPC));
            for u = 1:length(unq)
                clusterGM{u} = gmdistribution(gmd.mu(u,:),gmd.Sigma(:,:,u));
                mnClusterPC(u,:) = clusterGM{u}.mu;
            end
            UIDtranslation = zeros(length(units),1);
        end
        
        for u = 1:length(unq)
            if exist('inhull','file')
                clusPts = pc(assigns == unq(u),1:nPC);
                tol = 1.e-13*mean(abs(clusPts(:)));
                matchedIDs = inhull(pc_post(:,1:nPC),clusPts,[],tol);
            else
                %TODO: replace this error with a warning that we're limited
                %      to 3 dimensions, and use convhull instead. (Could
                %      use convhulln for higher dimensions, but
                %      inpolyhedron doesn't work beyond 3 dimensions.)
                error('No inhull function on path');
            end
            switch settings.prob_method
                case 'gmm'
                    if settings.joint
                        mnMatchPC = mean(pc_post(matchedIDs,1:nPC));
                        % need to find which gaussian was this one (I don't
                        % like this, it's a bit hacky. Keep thinking...)
                        dsts = pdist2(mnMatchPC,mnClusterPC);
                        [~,g] = min(dsts);
                        [~,wh] = max(id_probs(matchedIDs,:),[],2);
                        if mode(wh) ~= g
                            warning(['Mean PC scores of matched waveforms did not align with most likely cluster in channel ' num2str(chan) ' unit ' num2str(unq(u))]);
                            disp(['Minimum distance thought it was component ' num2str(g) ', modal match probability thought it was ' num2str(mode(wh)) '...'])
                            disp([9 '...going with the modal max match probabilities'])
                            g = mode(wh);
                        end
                        inner_probs = clusterGM{g}.pdf(pc_post(matchedIDs,1:nPC))./clusterGM{g}.pdf(clusterGM{g}.mu);
                        UIDtranslation(u) = units(unq(g)).UID;
                    else
                        % for each unit, fit the GMM on requested number of PCs, and calculate
                        % the probabilities for each new spike for that unit:
                        these_pc = pc(assigns == unq(u),1:nPC);
                        % fitgmdist is overkill when only using 1 cluster, but it's speedy:
                        gmd = fitgmdist(these_pc,1);
                        %{
                        % expand the fit by the requested amount:
                        big_gmd = gmdistribution(gmd.mu, settings.GMM_expansion * gmd.Sigma);
                        %}
                        % calculate the match probabilities, scaling the Gaussian such that
                        % a waveform at its exact mean has a match probability of 1:
                        inner_probs = gmd.pdf(pc_post(matchedIDs,1:nPC))./gmd.pdf(gmd.mu);
                        % unit identities are lost in the joint-GMM method,
                        % so we have to use g from here on out. This method
                        % doesn't lose identity, so g == u:
                        g = u;
                    end
                    
                case 'voltage'
                    matchW = spks(matchedIDs,:);
                    %TODO: move spk_gauss_probs to a subfunction, or inline:
                    inner_probs = spk_gauss_probs(units(u).waveforms,matchW);
                    % unit identities are lost in the joint-GMM method, so
                    % we have to use g from here on out. This method
                    % doesn't lose identity, so g == u:
                    g = u;
                otherwise
                    error(['Unknown probability calculation method: ' settings.prob_method])
            end
            
            matchW = spks(matchedIDs,:);
            matchT = spkt(matchedIDs);
            matchW(inner_probs < settings.threshold,:) = [];
            matchT(inner_probs < settings.threshold) = [];
            
            match_probs = inner_probs(inner_probs > settings.threshold);
            
            
            temp_unit = SingleUnit('patient',original.patient,'seizure',original.seizure);
            temp_unit.UID = units(unq(g)).UID;
            temp_unit.waveforms = matchW;
            temp_unit.times = matchT;
            temp_unit.channel = units(unq(g)).channel;
            temp_unit.type = units(unq(g)).type;
            temp_unit.extra.probabilities = match_probs;
            temp_unit.extra.mean_original_waveform = mean(units(unq(g)).waveforms);
            temp_unit.extra.sd_original_waveform = std(units(unq(g)).waveforms);
            
            if settings.joint && strcmp(settings.prob_method,'gmm')
                temp_unit.extra.raw_ID_prob_UIDs = UIDtranslation;
                temp_unit.extra.raw_ID_probs = id_probs(matchedIDs,:);
                % i.e. probability that spike came from raw_ID_prob_UIDs(1)
                % is raw_ID_probs(n,1)
            end
            
            matched.add_unit(temp_unit);
        end
        
    end
end