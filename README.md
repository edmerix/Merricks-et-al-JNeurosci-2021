# Code from Merricks et al., JNeurosci (2021)

A collection of MATLAB code used in, or derived from, our paper "[Neuronal Firing and Waveform Alterations through Ictal Recruitment in Humans](https://doi.org/10.1523/JNEUROSCI.0417-20.2020)", published in the Journal of Neuroscience.

## Overview

|             Function            |                                                                                                                Description                                                                                                                |
|:--------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ```template_match_convhull.m``` | The core convex hull template matching method, with multiple features. _Has dependencies (see below)_ |
| ```basic_convhull_match.m```    | A simplified version of convex hull template matching that only works in ≤ 3 dimensions, but has no dependencies.                                                                                                                         |
| ```spike_denoise.m```           | The novel method of automatically detecting electrical noise and transients in detected spikes using FFTs.                                                                                                                                |
| ```spk_gauss_probs.m```         | The method of building spike match confidences using scaled Gaussian fits.                                                                                                                                                                |
| ```gmm_match.m```               | The original Gaussian mixture model method of template matching _(outdated but functional!)_                                                                                                                                              |

### Notes:

```template_match_convhull``` is more full-featured, including options to calculate match confidence by wave shape or a Gaussian mixture model, but has the following dependencies:

- [NeuroClass](https://github.com/edmerix/NeuroClass)  (My object-oriented toolbox for analyses using populations of single units)
- [InHull.m](https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull) by John D'Errico, to allow for higher-dimensional convex hulls

To avoid these dependencies, if only using ≤ 3 dimensions in PC space, the simplified version, ```basic_convhull_match``` can be used instead, in combination with ```spk_gauss_probs``` or a Gaussian mixture model.

```gmm_match``` uses Gaussian mixture models to do template matching as well as the match confidences, bypassing the convex hull. This is useful in stable recordings, to track neurons in a probabilistic manner in the absence of noise that approaches cluster boundaries. __N.B.__ This function is in the process of being updated to allow for match confidences derived from a combination of the posterior probabilities and the original Gaussian fits themselves.


## Usage

In general, units need to be well-isolated in the peri-ictal period without any outliers to enable accurate template matching. Depending on the firing rates of the neurons, we find that typically requires between 10 and 30 minutes of continuous data. Due to waveform changes because of micro-motion of the cell bodies relative to the electrode tips, this should be done as close as possible to the period to be template matched.

For analyses using template matched data and the original epoch, the template matching should be run on both the _new and the original epoch_, to avoid biasing analyses by spike sorting method.

When denoising detected spikes with ```spike_denoise```, this should be done after spike detection but before clustering or template matching. The default parameters are set up to fairly comprehensively detect slow transients (< 500 Hz) or electrical noise (> 2.5 kHz) in our data, but the discarded waveforms should always be plotted and checked to make sure good waveforms are not being caught, depending on the quality of your recordings. The method uses standard deviations in your specified frequency ranges, so particularly clean or noisy recordings may trip it up.

See the help sections of each function for more in-depth notes.


## Bonus animation!

![Animation of convex hull matching in progress](superfluous/convexhullbuild.gif?raw=true "Convex hulls in action on PCA clusters")




