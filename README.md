# Threat-Brain-Activity

These MATLAB scripts concern the preprocessing, feature extraction and functional connectivity of EEG data acquired during the visualization of an unexpected visual threat implemented in VR. Below is a brief description of each of the scripts.



pre_processing.m
  - Removes excess data (time instants of recordings before and after the VR video presentation)
  - Applies filters to the EEG data (bandpass and notch filters)
  - Removes artifacts with ICA 


feature_extraction_functional_connectivity.m
  - Extracts relative PSD values
  - Applies the Kolmogorov-Smirnov test to find data distribution
  - Applies the Kruscal-Wallis test to find the significant EEG channels and frequency bands
  - Calculates the functional connectivity (COH, iCOH, wPLI and MPC) between the significant EEG channels (NOTE: DTF is calculated through the eConnectome MATLAB toolbox)


findDataDistribution.m
  - Applies the Kolmogorov-Smirnov test to find data distribution


findSignificantConnections.m
  - Applies the Kruscal-Wallis test to find the significant EEG channels connections
  - If wanted, performs a multiple comparation correction to the p-values obtained from the Kruscal-Wallis test


findHigherConnectivityConnections.m
  - Finds the channel pairs whose connectivity values is higher than a defined threshold
  - Plots the connectivity matrices and a brain connectivity schematic representation
    
