# Microphone spotformer
This code is the matlab implementation of a microphone spotformer. A microphone spotformer is basically a beamformer that attempts to attenuate or amplify sounds from certain regions-of-interest.

As input, the spotformer takes the target locations, the interfer locations and the microphone locations. Some additional user parameters are required as well, such as the size of the region. 
The inputs are used to construct a spatial correlation or covariance matrix. These matrices are used to minimise the energy from the interferer regions with respect to that from the target regions. 

A few examples are provided. 
- Example 1: A reverberant scenario with a single interferer and a single target source. The weights are precomputed. This is handy when needing to reuse the same weights for multiple sound excerots.
- Example 2: A reverberant scenario with a single interferer and a single target source. The weights are computed when using the spotformer. 

The covariance matrices are computed through spatial integration. For this, a Clenshaw-Curtis quadrature method is needed. I used the [Fast Clenshaw-Curtis quadrature function](https://nl.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature) by G. von Winckel, published under a permissive license. 

The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). 

The spotformer implemented here is my interpretation of the following paper:
J. Martinez, N. Gaubitch and W. B. Kleijn, "A robust region-based near-field beamformer," 2015 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), South Brisbane, QLD, Australia, 2015, pp. 2494-2498, doi: 10.1109/ICASSP.2015.7178420. keywords: {Robustness;Microphones;Speech;Reverberation;Arrays;Array signal processing;near-field beamformer;generalized eigenvalue problem;robust beamformer;microphone arrays;reverberation},


 
