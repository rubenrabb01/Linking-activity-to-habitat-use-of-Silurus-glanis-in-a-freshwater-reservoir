
# Catfish Movement in Římov Reservoir

> **Graphical Abstract (top)**  
> *Graphic depiction of the expected turning angle (radians) and step length (m) distributions for a three‑state Hidden Markov Model.*  

<p align="center">
  <img src="figures/graphical_abstract.png" alt="Graphical Abstract" width="80%"/>
</p>

---

## Table of Contents

- [Overview](#overview)
- [Study Area](#study-area)
- [Data & Telemetry](#data--telemetry)
- [Environmental Data](#environmental-data)
- [Analysis Workflow (End-to-End)](#analysis-workflow-end-to-end)
  - [4. Build & Select HMMs](#4-build--select-hmms)
- [References](#references)

---

## Overview

We analyzed the summer movement behavior of European catfish (*Silurus glanis*) in the Římov Reservoir (Czech Republic) using high-resolution acoustic telemetry and a three-state Hidden Markov Model (HMM).

---

## Analysis Workflow (End-to-End)

### 4. Build & Select HMMs

We fit a 3‑state HMM using gamma step lengths and von Mises turning angles. Covariates: cyclic time (cos, sin) and standardized epilimnion temperature. Models compared via AIC/log‑likelihood.

#### Model components and equations

### Model components and equations

**Emission distributions (per state \(s\)):**

$$
L_t \mid S_t=s \sim \text{Gamma}(k_s,\ \theta_s)
$$

$$
\Theta_t \mid S_t=s \sim \text{von Mises}(\mu_s,\ \kappa_s)
$$

**Time-of-day cyclic covariates (seconds since midnight \(t\)):**

$$
\mathrm{cos\_time}_t = \cos\!\left(\frac{2\pi t}{86400}\right), \quad
\mathrm{sin\_time}_t = \sin\!\left(\frac{2\pi t}{86400}\right)
$$

**Temperature standardization (for numerical stability):**

$$
\texttt{stand\_mean\_temp} = \frac{\texttt{temp\_mean} - \mu}{\sigma}
$$

where \(\mu\) and \(\sigma\) are the mean and SD of epilimnion temperatures.

**State process with covariates (multinomial logit):**

$$
\log\frac{p_{ij}(t)}{p_{i3}(t)} =
\beta_{0,ij}
+ \beta_{1,ij}\,\mathrm{cos\_time}_t
+ \beta_{2,ij}\,\mathrm{sin\_time}_t
+ \beta_{3,ij}\,\texttt{stand\_mean\_temp}_t
+ \beta_{4,ij}\,\mathrm{cos\_time}_t \cdot \texttt{stand\_mean\_temp}_t
+ \beta_{5,ij}\,\mathrm{sin\_time}_t \cdot \texttt{stand\_mean\_temp}_t
$$

for \(j=1,2\) (baseline \(j=3\)). Reduced models drop terms accordingly.

---

## References

- Albeke, S. E. (2025). *rKIN: (Kernel) Isotope Niche Estimation* (Version 1.0.4). [CRAN](https://cran.r-project.org/web/packages/rKIN/index.html)
- Bacheler, N. M., Michelot, T., Cheshire, R. T., & Shertzer, K. W. (2019). Fine-scale movement patterns of gray triggerfish *Balistes capriscus*. *Fisheries Research, 215*, 76–89. [https://doi.org/10.1016/j.fishres.2019.02.014](https://doi.org/10.1016/j.fishres.2019.02.014)
- Copp, G. H., Britton, J. R., Cucherousset, J., García-Berthou, E., Kirk, R., Peeler, E., & Stakėnas, S. (2009). A review of European catfish *Silurus glanis*. *Fish and Fisheries, 10*(3), 252–282.
