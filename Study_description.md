### Model components and equations

**Emission distributions (per state \(s\)):**

$$
L_t \mid S_t=s \sim \mathrm{Gamma}(k_s,\ \theta_s)
$$

$$
\Theta_t \mid S_t=s \sim \mathrm{von\ Mises}(\mu_s,\ \kappa_s)
$$

**Time-of-day cyclic covariates (seconds since midnight \(t\)):**

$$
\cos\_time_t = \cos\!\left(\frac{2\pi t}{86400}\right), \quad
\sin\_time_t = \sin\!\left(\frac{2\pi t}{86400}\right)
$$

**Temperature standardization (for numerical stability):**

$$
stand\_mean\_temp = \frac{temp\_mean - \mu}{\sigma}
$$

where \(\mu\) and \(\sigma\) are the mean and standard deviation of epilimnion temperatures.

**State process with covariates (multinomial logit):**

$$
\begin{aligned}
\log\frac{p_{ij}(t)}{p_{i3}(t)} = {} & \beta_{0,ij} \\
 &+ \beta_{1,ij}\cos\_time_t \\
 &+ \beta_{2,ij}\sin\_time_t \\
 &+ \beta_{3,ij}stand\_mean\_temp_t \\
 &+ \beta_{4,ij}\cos\_time_t \cdot stand\_mean\_temp_t \\
 &+ \beta_{5,ij}\sin\_time_t \cdot stand\_mean\_temp_t
\end{aligned}
$$



### 4. Build & Select HMMs

We fit a 3-state HMM using gamma step lengths and von Mises turning angles. To determine the influence of temperature on fish activity, epilimnion temperature was included in the state process of the HMM after standardization (variable `stand_mean_temp`):

$$
stand\_mean\_temp = \frac{temp\_mean - \mu}{\sigma}
$$

where \(\mu\) and \(\sigma\) are the mean and standard deviation of epilimnion temperatures.

In addition, time of day was included using trigonometric transformations of seconds since midnight (\(t\)):

$$
tod\_trigo = \cos\!\left(\frac{2\pi t}{86400}\right) \;+\; \sin\!\left(\frac{2\pi t}{86400}\right)
$$

Thus, the following five HMMs were fit using the **moveHMM** package (Michelot et al. 2016) and compared based on AIC and log-likelihood:

1. No covariates  
2. `stand_mean_temp`  
3. `tod_trigo`  
4. `tod_trigo + stand_mean_temp`  
5. `tod_trigo * stand_mean_temp`

The model which best predicted the activity states included an interaction between cyclic time and epilimnion temperature (see Appendix I, Table A1.2). The most likely sequence of states and associated probabilities were extracted from the selected model using:

- `viterbi()` — implements the Viterbi algorithm (Zucchini et al., 2016)  
- `stateProbs()` — computes state probabilities  

Weeks and individuals where only two states were identified were removed from further analysis.
