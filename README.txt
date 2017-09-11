The folder contains the MATLAB codes to extract features using unsupervised data-
driven modulation filtering approach (for noise robust speech recognition system).

- The modulation filters (rate and scale separately) are learned from mel spectrograms 
  using Convolutional Restricted Boltzmann Machine (CRBM).
- Multiple filters are learned using residual approach.
- The filter selection criteria uses average hidden activation probability values.
- The mel spectrograms are filtered using the selected modulation filters 
  which are then fed as features to train DNN for building ASR.

Sequence to run :

0. Zero_example_mel_spectrogram_speech.m
1. One_example_filtLearning_Rate_crbm.m
2. Two_example_filtLearning_Scale_crbm.m
3. Three_example_filtSelection_validation.m
4. Four_hidden_prob_avg.m
5. Five_example_feature_extraction_speech_forASR.m