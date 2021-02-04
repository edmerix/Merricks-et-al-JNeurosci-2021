# Code from Merricks et al., JNeurosci (2021)

A collection of MATLAB code used in, or derived from, our paper "[Neuronal Firing and Waveform Alterations through Ictal Recruitment in Humans](https://doi.org/10.1523/JNEUROSCI.0417-20.2020)", published in the Journal of Neuroscience.

## Overview

|             Function            |                                                                                                                Description                                                                                                                |
|:--------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ```template_match_convhull.m``` | The core convex hull template matching method, with multiple features. _Has dependencies (see below)_ 		                                                                                                                              |
| ```basic_convhull_match.m```    | A simplified version of convex hull template matching that only works in ≤ 3 dimensions, but has no dependencies.                                                                                                                         |
| ```spike_denoise.m```           | The novel method of automatically detecting electrical noise and transients in detected spikes using FFTs.                                                                                                                                |
| ```spk_gauss_probs.m```         | The original method of building spike match confidences, using scaled Gaussian fits.                                                                                                                                                      |
| ```spk_match_confidence.m```      | The new method of building spike match confidences, using the probability of finding a value further from the mean in the original unit's Gaussian fit at each data point                                                                 |
| ```gmm_match.m```               | The original Gaussian mixture model method of template matching _(outdated but functional!)_                                                                                                                                              |

### Notes:

```template_match_convhull()``` is more full-featured, including options to calculate match confidence by wave shape or a Gaussian mixture model, but has the following dependencies:

- [NeuroClass](https://github.com/edmerix/NeuroClass)  (My object-oriented toolbox for analyses using populations of single units)
- [InHull.m](https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull) by John D'Errico, to allow for higher-dimensional convex hulls

To avoid these dependencies, if only using ≤ 3 dimensions in PC space, the simplified version, ```basic_convhull_match()``` can be used instead, in combination with ```spk_gauss_probs()```, ```spk_match_confidence()``` or a Gaussian mixture model.

```gmm_match()``` uses Gaussian mixture models to do template matching as well as the match confidences, bypassing the convex hull. This is useful in stable recordings, to track neurons in a probabilistic manner in the absence of noise that approaches cluster boundaries. __N.B.__ This function is in the process of being updated to allow for match confidences derived from a combination of the posterior probabilities and the original Gaussian fits themselves.

The original ```spk_gaus_probs()``` is what was used in the original paper, but calculates match confidences based on scaled Gaussians to attempt to approximate probabilities of matches. This is necessary because the a value found in the fitted Gaussian corresponds to the probability of picking that waveform from the _true_ original waveforms, __not__ the probability that a new waveform came from that neuron. 

```spk_match_confidence()``` updates this approach to make better use of the Gaussian distributions of the original neuron's waveforms without scaling them. It instead calculates what the probability of selecting a true match from the original unit _farther_ from the mean at each distribution (in either direction) is. See [Fig. 1](#fig1) immediately below:

<a name="fig1">![Explanation of spk_match_confidence.m method](superfluous/gauss_probs.png?raw=true "Explanation behind spk_match_confidence method")</a>

<small>__Figure 1.__ _Match confidences using ```spk_match_confidence()```_</small>

<small>__A.__ Gaussian distributions are fitted to the range of voltages in the original neuron's waveforms at each datapoint. The original waveforms are seen below the 3D plot in shades of blue, with the PDF of each Gaussian overlaid in the 3D heat map. An example new waveform that is a "good" match to this neuron is shown overlaid on the 3D map in black, and a "bad" one in red. __B. i.__ Overview of how the probabilities are calculated for each distribution at each datapoint. The true voltage value for the black and red waveforms are shown on each distribution (large dots) with their distance from the mean value (stems). The area under the Gaussian beyond both this point and the equivalent distance on the opposite side of the mean (small dots; both areas shown with shading), gives the probability of finding a waveform from the original distribution with a voltage beyond this distance from the mean. __ii.__ A zoom of the datapoint highlighted with the triangle in i. In this example we can see the proximity of the black waveform's value to the original neuron's mean at this datapoint (blue dashed line), and the area under the PDF beyond this distance from the mean is thus larger (black shading, which continues under the red shading). Conversely, the red waveform's value is farther from the mean, producing a lower probability of selecting a waveform from the original dataset with a value farther from the mean. These give values of 0.56 and 0.07 for the black and red waveforms respectively at this datapoint. The overall match confidence can be calculated as the mean of these areas over all datapoints. </small>

## Usage

In general, units need to be well-isolated in the peri-ictal period without any outliers to enable accurate template matching. Depending on the firing rates of the neurons, we find that typically requires between 10 and 30 minutes of continuous data. Due to waveform changes because of micro-motion of the cell bodies relative to the electrode tips, this should be done as close as possible to the period to be template matched.

For analyses using template matched data and the original epoch, the template matching should be run on ___both__ the new and the original epoch_, to avoid biasing analyses by spike sorting method.

When denoising detected spikes with ```spike_denoise()```, this should be done after spike detection but before clustering or template matching. The default parameters are set up to fairly comprehensively detect slow transients (< 500 Hz) or electrical noise (> 2.5 kHz) in our data, but the discarded waveforms should always be plotted and checked to make sure good waveforms are not being caught, depending on the quality of your recordings. The method uses standard deviations in your specified frequency ranges, so particularly clean or noisy recordings may trip it up.

See the help sections of each function for more in-depth notes.


## Bonus animation!

![Animation of convex hull matching in progress](superfluous/convexhullbuild.gif?raw=true "Convex hulls in action on PCA clusters")

