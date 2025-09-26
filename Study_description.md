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
\beta_{0,ij}
+ \beta_{1,ij}\,\text{cos\_time}_t
+ \beta_{2,ij}\,\text{sin\_time}_t
+ \beta_{3,ij}\,\text{stand\_mean\_temp}_t
+ \beta_{4,ij}\,\text{cos\_time}_t \cdot \text{stand\_mean\_temp}_t
+ \beta_{5,ij}\,\text{sin\_time}_t \cdot \text{stand\_mean\_temp}_t
$$

for \(j=1,2\) (baseline \(j=3\)). Reduced models drop terms accordingly.
