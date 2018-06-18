# volumelearning
Matlab scripts to demonstrate learning spatial and temporal signals in spiking networks. Code is structured around the main file [volumeNetwork](volumeNetwork.m) which can run when invoked with a matlab interpreter.

More specificallly this is a two layer network of LIF neurons with gaussian current based, time delayed synapses. Ideas and inspiration are taken from [1, 2].
General information about quirks and issues during development can be found on [this google doc](https://docs.google.com/document/d/1mJY0HtCmyt_8qT5UmYbZxuVqA2B_u3qlE3r0deTQxLY/edit?usp=sharing).

## Interesting milestones
- Demo of Multiple gaussian output synapses [working](https://github.com/jotia1/volumelearning/commit/ba17a45eb05bd605f451f3f4510fb89df9917564)
- STDP and Synaptic scaling working [demo](https://github.com/jotia1/volumelearning/commit/a4d03ad126d3b05776fe476b21816eaa8e9c9952)

## Current features
- [x] Gaussian synapses
- [x] Multiple outputs
- [x] STDP
- [ ] Delay Learning 
- [ ] Lateral inhibition in output
- [ ] Parameter sweep scripts
- [ ] Event-based data
- [ ] Party ... wait no, publish

## Useful spiking resources
Additionally I'd like to acknowledge other bits of code that inspired/influenced this implementation as there is a need for more accessible resources around spiking networks (especially those with delays). 
In particular the [Neural Data Modeling Group's tutorial](http://www.mjrlab.org/2014/05/08/tutorial-how-to-write-a-spiking-neural-network-simulation-in-matlab-from-scratch/) and Izhikevich's [Polychronous code](https://www.izhikevich.org/publications/spnet.htm). [3]

## References
[1] P. W. Wright and J. Wiles, “[Learning transmission delays in spiking neural networks: A novel approach to sequence learning based on spike delay variance](https://ieeexplore.ieee.org/document/6252371/)” in The 2012 International Joint Conference on Neural Networks (IJCNN), 2012, pp. 1–8.\
[2] T. Matsubara, “[Conduction Delay Learning Model for Unsupervised and Supervised Classification of Spatio-Temporal Spike Patterns](https://ieeexplore.ieee.org/abstract/document/7966073/)” Front. Comput. Neurosci., vol. 11, 2017.\
[3] E. M. Izhikevich, “[Polychronization: Computation with Spikes](https://www.izhikevich.org/publications/spnet.htm)” Neural Computation, vol. 18, no. 2, pp. 245–282, Feb. 2006.
