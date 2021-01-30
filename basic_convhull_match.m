function [matched, pc_pre, pc_post] = basic_convhull_match(original,new,uid,centered)
% 3D convex hull matching. Use convhull_match_ND to work in > 3 dimensions.
%   Takes two necessary inputs: original and new (described below), an 
%   optional third for the UID (assigned value) in the original dataset to 
%   match to (allowing passing in a full population of cells for 
%   calculating the PCA, then only using one cell within that population to
%   match to), and an optional fourth to determine if PCA should auto-
%   center the data.
%
%   Returns a struct containing the waveforms and times of the matched data
%   
%   Input struct design:
%       original.waveforms  = all waveforms to calculate PCA on
%       original.times      = spiketimes of those waveforms
%       original.assigns    = the assignment IDs for each waveform
%                           
%       new.waveforms       = all waveforms that you want to find matches
%                             within
%       new.times           = the spike times of those waveforms
%
% For a more full-featured version that works in higher than 3-dimensional
% space, see the function "template_match_convhull" in the same repository
% as this function, which also calculates waveform match probabilities, but
% requires extra dependencies.
%
% E. M. Merricks, Ph.D. 2019-05-09

if nargin < 3 || isempty(uid)
    unit_inds = 1:length(original.times);
else
    unit_inds = find(original.assigns == uid);
end

if nargin < 4 || isempty(centered)
    centered = true;
end

% calculate the PCA on the original units:
[coef, pc_pre] = pca(original.waveforms,'Centered',centered);
% calculate the PC scores for the new waveforms on the original unit's PCA
% space:
if centered
    pc_post = new.waveforms * coef;
else
    post_waves = new.waveforms - mean(new.waveforms,1);
    pc_post = post_waves * coef;
end

fv = []; 
fv.faces = convhull(pc_pre(unit_inds,1),pc_pre(unit_inds,2),pc_pre(unit_inds,3),'simplify',true);
fv.vertices = pc_pre(unit_inds,1:3);

matchedIDs = inpolyhedron(fv,pc_post(:,1:3));

matched.waveforms = new.waveforms(matchedIDs == 1,:);
matched.times = new.times(matchedIDs == 1);
matched.indices = find(matchedIDs == 1);