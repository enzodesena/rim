# (Randomized) Image Method

This project provides a Matlab implementation of the **image method**, a widely used acoustic model to calculate the room impulse response (RIR) of a rectangular room. The software allows multiple observation points and fractional delays.

If needed, the software generates can randomise the position of the image sources to avoid the occurrence of **sweeping echoes**, resulting in the so-called **randomized image method**. See ``De Sena et al. "On the modeling of rectangular geometries in  room acoustic simulations." IEEE/ACMT TASLP (2015).``


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This version has been tested on Matlab2020a, but I have no reason to believe it wouldn't work on previous versions too. 
No additional toolbox is needed. 

### Installing

No installation needed. 

## Running the tests

To run the tests simply run

```
rim_tests
```

## Using the software

The main script is `rim.m`:

``RIM(mic_pos, source_pos, room_dim, beta, rir_length, Fs)``
- `mic_pos` is the 3xM matrix with the position of M omnidirectional microphones. The dimenions are in [meters]. 
- `source_pos` is the 3x1 vector with the position of the omni sound source [meters]. 
- `room_dim` is the 3x1 vector with the dimensions of the room [meters].
- `beta` is the 2x3 vector with the reflection coefficient of the walls, in this order: [x1, y1, z1; x2, y2, z2], with e.g. x1 and x2 are the reflection coefficient of the surface orthogonal to the x-axis and with the subscript 1 referring to walls adjacent to the coordinate origin, and the subscript 2 referring to the opposite wall. In the anechoic case beta=zeros(2,3). 
- `rir_length` is the length of the RIR [seconds].
- `Fs` is the sampling frequency [Hz].

The following are optional inputs:

``RIM(mic_pos, source_pos, room_dim, beta, rir_length, Fs, rand_dist, Tw, Fc, c)``
- `rand_dist` is the random distance added to the position of the image sources in [meters]; in the Transaction paper we found that rand_dist=0.08 m was sufficient to remove sweeping echoes in most cases studied (rand_dist=0 for the standard image method) (default is 0 cm).
- `Tw` is the length of the fractional delay filter in [seconds] (default is Tw=0.004 s)
- `Fc` is the cut-off frequency of the fractional delay filter in [Hz] (default is Fc=0.9*(Fs/2))
- `c` is the speed of sound in m/s (default is c=343)

## Generating the plots in the 2015 Transaction paper

Run ``rim_taslp_plots``.

Note with regard to Figure 11a (Highly regular setup with Randomized Image Method): this repository uses a different randomisation strategy compared to the paper. While in the paper the source was moved on a line between the original image sound and the observer by +-8 cm (i.e. 16 cm spread), here, the source is moved within a cube with edge 16 cm and centered at the original image source. Visual inspection suggests the effect is similar. If you'd like to replicate this figure exactly as in the Transaction paper please refer to [desena.org/sweep/](https://desena.org/sweep/).

## Contributing

Contributions are welcome. 

## Author

* **Enzo De Sena** - [desena.org](https://desena.org)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details


If you use this code, please cite the following paper:

``De Sena, E., Antonello, N., Moonen, M. and Van Waterschoot, T., 2015. On the modeling of rectangular geometries in room acoustic simulations. IEEE/ACM Transactions on Audio, Speech, and Language Processing, 23(4), pp.774-786.`` 

[IEEE Xplore](https://ieeexplore.ieee.org/document/7045580)
