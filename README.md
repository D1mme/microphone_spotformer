# Microphone spotformer
This code is the matlab implementation of a microphone spotformer. A microphone spotformer is basically a beamformer that attempts to attenuate or amplify sounds from certain regions-of-interest.

As input, the spotformer takes the target locations, the interfer locations and the microphone locations. Some additional user parameters are required as well, such as the size of the region. 
The inputs are used to construct a spatial correlation or covariance matrix. These matrices are used to minimise the energy from the interferer regions with respect to that from the target regions. 

A few examples are provided. 
- Example 1: An anechoic scenario with a single interferer and a single target source. The weights are precomputed. This is handy when needing to reuse the same weights for multiple sound excerots.
- Example 2: A slightly reverberant (T60 $\approx$ 120 ms)  scenario with a single interferer and a single target source. The weights are computed when using the spotformer. 
- Example 3: A more reverberant scenario (T60 $\approx$ 300 ms) with one interfer and one target source.
- Example 4: A slightly reverberant scenario (T60 $\approx$ 120 ms) with four interferers and two target sources.
All examples use eight microphones and have white Gaussian noise at the interferers. Note that other types of playback signal and different numbers of microphones also work. 

The covariance matrices are computed through spatial integration. For this, a numerical integration method is needed. I used the [Fast Clenshaw-Curtis quadrature function](https://nl.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature) by G. von Winckel, published under a permissive license. Note that this is a different integration method than that used in the original paper ([1]). I did some small experiments (nothinh on wbhich you can draw a definitive conclusion!) and it appeard that the Clenshaw-Curtis quadrature gave more accurate results than the Gauss-Hermite quadrature method used in [1] (though it should be noted that they integrate over a different range). 

The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). 

The spotformer implemented in this repisotory is my interpretation of the following paper:
[1] J. Martinez, N. Gaubitch and W. B. Kleijn, "A robust region-based near-field beamformer," <em> 2015 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) <em>, South Brisbane, QLD, Australia, 2015, pp. 2494-2498, doi: 10.1109/ICASSP.2015.7178420.


 
