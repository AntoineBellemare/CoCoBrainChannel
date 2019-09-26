# CoCoBrainChannel
Pure Data and Python tools for creative and aesthetic neurofeedback

This repository is a work-in-progress.
It builds on EEGsynth architecture to easily compute features of brain signals and communicate to PureData (or any OSC compatible software) 
in order to facilitate the fluency of complex computational processes of EEG signal (Python) with real-time audio processing (PureData).

The BioTuning Library makes uses of brain signal spectral information (spectral peaks) to easily construct musical harmonic structures.  

PureData Installation.
Works with Pd Vanilla.
Required libraries:
- Cyclone
- List-abs
- Mrpeach
- iemlib

CCBC patch is the main sonification patch.
CCBC_optimize is an optimized version of the pathc to ensure CPU efficiency.
