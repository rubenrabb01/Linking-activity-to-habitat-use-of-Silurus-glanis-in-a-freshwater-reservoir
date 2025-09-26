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
\log\frac{p_{ij}(t)}{p_{i3}(t)} =
\beta_{0,ij} + \beta_{1,ij}\cos\_time_t
+ \beta_{2,ij}\sin\_time_t
+ \beta_{3,ij}stand\_mean\_temp_t
+ \beta_{4,ij}\cos\_time_t \cdot stand\_mean\_temp_t
+ \beta_{5,ij}\sin\_time_t \cdot stand\_mean\_temp_t
$$

