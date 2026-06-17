# ECG-Dynamic-Model-Parameter-Extraction-Using-PSO
ECG Dynamic Model Parameter Extraction Using PSO


This repository provides a MATLAB implementation for modeling ECG morphology using a Gaussian Mixture Model (GMM) and Particle Swarm Optimization (PSO).

The method detects cardiac cycles, converts the ECG signal to a phase representation, estimates the average ECG morphology, and approximates it using a set of Gaussian kernels optimized with Particle Swarm Optimization (PSO).

The resulting Gaussian parameters can be used to reconstruct or synthesize ECG signals and analyze ECG morphology in a compact parametric form.
Algorithm Pipeline

The workflow of the proposed approach is summarized below:

    Load ECG signal
    R‑peak detection using the Pan–Tompkins algorithm
    Phase estimation based on RR intervals
    Phase-domain ECG averaging
    Gaussian kernel extraction using PSO
    Selection of strongest Gaussian components
    Synthetic ECG reconstruction

The method models the ECG morphology as a sum of Gaussian functions defined in the cardiac phase domain.
Mathematical Model

The ECG morphology in the phase domain is approximated by a sum of Gaussian kernels:

Z(θ)=∑i=1Naiexp⁡(−(θ−θi)22bi2)
Z(θ)=i=1∑N​ai​exp(−2bi2​(θ−θi​)2​)

Where:

    aiai​ : amplitude of the Gaussian component
    bibi​ : width (standard deviation)
    θiθi​ : center phase of the component
    NN : number of Gaussian components

These parameters describe the morphology of the ECG waveform (P wave, QRS complex, T wave, etc.).
Input Data

The program expects a MATLAB .mat file containing:

    x → ECG signal matrix
    fs → sampling frequency (Hz)

Example:

                                                                    text
x : [channels × samples]
fs : sampling frequency

The first channel of x is used as the ECG signal.
Main Processing Steps
1. R‑Peak Detection

R‑peaks are detected using the Pan–Tompkins algorithm, a well‑known method for QRS detection.

This step determines cardiac cycle boundaries.
2. Phase Calculation

The ECG signal is transformed into a phase representation using RR intervals.

Each heartbeat is mapped into the interval:

[−π,π]
[−π,π]

This representation makes it easier to analyze ECG morphology independently of heart rate.
3. Phase‑Domain ECG Averaging

The ECG signal is divided into phase bins and the following statistics are computed:

    Mean ECG waveform
    Standard deviation per phase bin

This produces a phase‑aligned average ECG morphology.
4. Gaussian Parameter Estimation

The ECG mean waveform is approximated by multiple Gaussian kernels.

Optimization is performed using Particle Swarm Optimization (PSO) to estimate:

    Gaussian amplitudes ai
    Gaussian widths bi
    Gaussian phase centers θi

Initially many Gaussian kernels are extracted, and then the strongest components are selected.

MaxNumGaussian = 25

5. Synthetic ECG Generation

Using the estimated Gaussian parameters, a synthetic ECG waveform is generated.

This synthetic signal approximates the morphology of the original ECG.
Output

The script produces several figures

The optimization framework used in this repository is related to our research on ECG modeling and enhancement using Particle Swarm Optimization within a Bayesian framework.

If you use this code in your research, please cite:

H. Danandeh Hesar and A. Danandeh Hesar

“ECG enhancement using a modified Bayesian framework and particle swarm optimization”

Biomedical Signal Processing and Control, Vol. 80, 104280, 2023.

https://doi.org/10.1016/j.bspc.2022.104280
