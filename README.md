:::writing
# ECG Modeling Using Gaussian Mixture in Phase Domain

This repository provides a MATLAB implementation for **modeling ECG morphology in the phase domain using a Gaussian Mixture Model (GMM)**.  
The method detects cardiac cycles, converts the ECG signal to a **phase representation**, estimates the **average ECG morphology**, and approximates it using a set of **Gaussian kernels optimized with Particle Swarm Optimization (PSO)**.

The resulting Gaussian parameters can be used to **reconstruct or synthesize ECG signals** and analyze ECG morphology in a compact parametric form.

---

# Algorithm Pipeline

The workflow of the proposed approach is summarized below:

1. **Load ECG signal**
2. **R‑peak detection** using the Pan–Tompkins algorithm  
3. **Phase estimation** based on RR intervals  
4. **Phase-domain ECG averaging**  
5. **Gaussian kernel extraction using PSO**  
6. **Selection of strongest Gaussian components**  
7. **Synthetic ECG reconstruction**

The method models the ECG morphology as a **sum of Gaussian functions defined in the cardiac phase domain**.

---

# Mathematical Model

The ECG morphology in the phase domain is approximated by a sum of Gaussian kernels:

``` 
\[
Z(\theta) = \sum_{i=1}^{N} a_i \exp\left(-\frac{(\theta-\theta_i)^2}{2b_i^2}\right)
\]
`
Where:

- \(a_i\) : amplitude of the Gaussian component  
- \(b_i\) : width (standard deviation)  
- \(\theta_i\) : center phase of the component  
- \(N\) : number of Gaussian components  

These parameters describe the morphology of the ECG waveform (P wave, QRS complex, T wave, etc.).

---

# Input Data

The program expects a MATLAB `.mat` file containing:

- **x** → ECG signal matrix  
- **fs** → sampling frequency (Hz)

Example:

```
x : [channels × samples]
fs : sampling frequency
```

The first channel of `x` is used as the ECG signal.

---

# Main Processing Steps

## 1. R‑Peak Detection

R‑peaks are detected using the **Pan–Tompkins algorithm**, a well‑known method for QRS detection.

This step determines cardiac cycle boundaries.

---

## 2. Phase Calculation

The ECG signal is transformed into a **phase representation** using RR intervals.

Each heartbeat is mapped into the interval:

\[
[-\pi, \pi]
\]

This representation makes it easier to analyze ECG morphology independently of heart rate.

---

## 3. Phase‑Domain ECG Averaging

The ECG signal is divided into **phase bins** and the following statistics are computed:

- Mean ECG waveform
- Standard deviation per phase bin

This produces a **phase‑aligned average ECG morphology**.

---

## 4. Gaussian Parameter Estimation

The ECG mean waveform is approximated by multiple Gaussian kernels.

Optimization is performed using **Particle Swarm Optimization (PSO)** to estimate:

- Gaussian amplitudes \(a_i\)
- Gaussian widths \(b_i\)
- Gaussian phase centers \(θ_i\)

Initially many Gaussian kernels are extracted, and then the **strongest components** are selected.

```
MaxNumGaussian = 25
```

---

## 5. Synthetic ECG Generation

Using the estimated Gaussian parameters, a **synthetic ECG waveform** is generated.

This synthetic signal approximates the morphology of the original ECG.

---

# Output

The script produces several figures:

**Figure 1**

ECG signal with detected R‑peaks.

**Figure 2**

Comparison between:

- Original ECG mean waveform
- Gaussian mixture approximation.

**Figure 3**

Comparison between:

- Original ECG signal
- Synthetic ECG signal generated from Gaussian parameters.

---

# Saved Parameters

The script stores estimated parameters in a MATLAB file:

```
*_params_linear_Phase.mat
```

Saved variables include:

- `OptimumParams` → Gaussian parameters ``` math[a_i , b_i , θ_i] `
- `ECGmean` → ECG mean waveform in phase domain
- `fs` → sampling frequency
- `x` → original ECG signal

---

# Requirements

- MATLAB
- Global Optimization Toolbox (for `particleswarm`)
- Pan–Tompkins QRS detection implementation

---

# Applications

This framework can be used in several ECG processing tasks:

- Synthetic ECG generation  
- ECG morphology modeling  
- ECG compression  
- Biomedical signal simulation  
- ECG denoising frameworks  
- Cardiac dynamics analysis  

---

# Related Publication

The optimization framework used in this repository is related to our research on ECG modeling and enhancement using **Particle Swarm Optimization within a Bayesian framework**.

If you use this code in your research, please cite:

H. Danandeh Hesar and A. Danandeh Hesar  
**“ECG enhancement using a modified Bayesian framework and particle swarm optimization”**  
Biomedical Signal Processing and Control, Vol. 80, 104280, 2023.

https://doi.org/10.1016/j.bspc.2022.104280

---

# BibTeX

```
@article{danandeh2023ecg,
  title={ECG enhancement using a modified Bayesian framework and particle swarm optimization},
  author={Danandeh Hesar, Hamed and Danandeh Hesar, Amin},
  journal={Biomedical Signal Processing and Control},
  volume={80},
  pages={104280},
  year={2023},
  publisher={Elsevier}
}
```

---

# License

This project is intended for **research and educational purposes**.
:::
