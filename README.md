# Animation for Lattice QCD

2023/10/30 akio@yukawa.kyoto-u.ac.jp A. Tomiya

This is a package for making an animation of a lattice configuration.

- Input format: ILDG, JLD2 (small modification is needed. See [here](https://github.com/akio-tomiya/Gaugefields.jl) )

- Output: GIF.

Plaquette for each site $(x,y)$ is plotted and interpolated by B-spline for presentation.

$\mu$ - $\nu$ and $z$ axis are averaged.

The $t$ direction is imaginary-time but treated as a real-time.



A sample file is for L16160432_beta6.1 with the quenched Wilson plaquette action, namely $a\times N_x \approx 1$ fm.

The scale $a(\beta=6.1)$ is taken from arXiv:hep-lat/9806005.



A file ``misc/my_parameters_large.toml`` is a parameter file and a configuration for LatticeQCD.jl and this will generate the sample used in above. The configuration is in ILDG format.

Ref: https://arxiv.org/abs/hep-lat/0004025 
