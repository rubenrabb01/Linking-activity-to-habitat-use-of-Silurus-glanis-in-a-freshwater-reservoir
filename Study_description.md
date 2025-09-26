# Linking activity to habitat use of *Silurus glanis* in a freshwater reservoir

![Graphical Abstract](figures/Graphical_Abstract.png)  
**Graphical Abstract.** Conceptual depiction of step length (m) and turning angle (radians) distributions for a three-state Hidden Markov Model (inactive, low, high activity). Adapted from Rabaneda-Bueno, unpublished.

---

## Methods

### Study area
Římov Reservoir (South Bohemia, Czech Republic; 48°51.00978′ N, 14°29.47462′ E) is a canyon-shaped reservoir built in 1978 for drinking water storage and flood control. It extends ~10 km, covers 210 ha, and reaches 45 m maximum depth (Figure 1). The system shows mesotrophic–eutrophic gradients (Šeda & Devetter 2000), strong summer stratification (Říha et al. 2022), and fluctuating water levels that prevent littoral macrophyte growth (Říha et al. 2015).  
The estimated catfish population (~200 individuals) and restricted human access make this an ideal system for studying natural predator behaviour.  

During the study (July–August 2017), epilimnion temperature ranged between 19.7–22.7 °C (mean: 21.0 ± 0.6 °C). Thermocline depth averaged 5.5 m.

![Figure 1](figures/Figure%201.png)  
**Figure 1.** Placement of telemetry receivers in Římov Reservoir, Czech Republic. Blue shading represents bathymetry. Adapted from Říha et al. 2025a.

---

### Data collection
- **Telemetry array:** 91 Lotek WHS3250 receivers covering reservoir, tributary, and bay.  
- **Temperature:** HOBO loggers deployed every 1 m (0–13 m, plus 20 m) at four sites.  
- **Fish capture & tagging:** 15 *S. glanis* captured by long-lining, anaesthetized (2-phenoxy-ethanol, 0.7 ml l−1), and tagged with Lotek MM-M-11-28-TP transmitters (13 g, burst rate 15 s).  
- **Final dataset:** 13 individuals retained after filtering.

---

### Data processing
Positions estimated with Lotek UMAP v1.4.3, filtered, and averaged at 5-min intervals.  
Habitat covariates calculated:  

```r
# Example processing code
library(rgeos)
# Distance to shore
positions$dist_shore <- gDistance(positions, reservoir_shoreline, byid = TRUE)

# Distance to reservoir bottom
positions$dist_bottom <- positions$res_depth - positions$fish_depth
```

- Missing values linearly interpolated for ≤1.5 h gaps.  
- Analysis period: July–August (peak activity, excluding spawning).  

---

### Hidden Markov Models (HMMs)
To classify activity states from step lengths and turning angles, we fitted HMMs using `moveHMM`. Missing data were interpolated using CTCRW with `momentuHMM::crawlWrap()`.

```r
library(moveHMM)
library(momentuHMM)

# Prepare track with CTCRW interpolation
interp_track <- crawlWrap(data, theta = c(log(10), log(1)))

# Define HMM distributions
step_dist <- "gamma"
angle_dist <- "vm"

# Fit 3-state HMM with covariates
hmm_model <- fitHMM(interp_track,
                    nbStates = 3,
                    dist = list(step = step_dist, angle = angle_dist),
                    formula = ~ cos_time + sin_time * stand_mean_temp)
```

- Candidate models:  
  1. Null  
  2. Temperature only  
  3. Time of day only  
  4. Additive (temp + time)  
  5. Interaction (temp × time)  

- Best model selected by AIC = interaction.  
- Most likely state sequences estimated with `viterbi()`.  

---

### Cumulative Link Mixed Models (CLMMs)
To link decoded HMM states (inactive < low < high) to habitat covariates, we fitted CLMMs with `ordinal`. Fish ID included as a random intercept.

```r
library(ordinal)

clmm1 <- clmm(state ~ scale(dist_shore) * diel + (1|fish_id), data = df)
clmm2 <- clmm(state ~ scale(botdepth) * diel + (1|fish_id), data = df)
clmm3 <- clmm(state ~ scale(dist_bottom) * diel + (1|fish_id), data = df)
clmm4 <- clmm(state ~ scale(fish_depth) * diel + (1|fish_id), data = df)
```

- Separate CLMMs for each habitat covariate.  
- Diel period (day/night) as effect modifier.  
- Model selection by AIC.  

---

### Three-dimensional habitat use
We estimated 50% (core) and 95% (home range) kernel utilization distributions (KUDs) per state × diel period using the `rKIN` package.

```r
library(rKIN)

kud50 <- kernelUD(track, h = "href", grid = 100)
kud95 <- kernelUD(track, h = "href", grid = 100, percent = 95)
```

- Overlaps calculated between states (day vs night).  
- Results plotted as depth–bottom cross-sections.

---

## Results

### Activity states
Three distinct activity states:  
- **Inactive:** short steps (2.4 ± 2.2 m), wide turning angles.  
- **Low activity:** medium steps (21.3 ± 15.3 m).  
- **High activity:** long steps (72.1 ± 29.5 m), directed movement.  

![Figure 5](figures/Figure%205.png)  
**Figure 5.** Final HMM model: (a) step length and (b) turning angle distributions per state.

---

### Habitat use and overlap
State home range areas increased with activity (inactive = 8.0 ha; low = 37.8 ha; high = 74.5 ha).  

![Figure 7](figures/Figure%207.png)  
**Figure 7.** Horizontal area used per state (N = 13). * indicates significant differences (ANOVA, p < 0.001).

---

### Time allocation
Most individuals spent 30–48% of time inactive, concentrated in few areas.  

![Figure 9](figures/Figure%209.png)  
**Figure 9.** Proportion of detections per state for each fish. Totals per individual shown above bars.

---

### Temperature and diel effects
At mean temp (21 °C): inactivity dominated day, activity peaked at night. Warmer conditions extended nocturnal activity; cooler advanced transitions.  

![Figure 11](figures/Figure%2011.png)  
**Figure 11.** State probabilities across diel cycle at low (19.7 °C), mean (21 °C), and high (22.7 °C) epilimnion temperatures. Grey = twilight; light grey = sunrise/sunset.

---

### Habitat CLMM results
Activity increased with distance from shore, bottom depth, and distance to bottom, but effects weakened at night. Absolute fish depth had weaker effects.

![Figure 12](figures/Figure%2012.png)  
**Figure 12.** Predicted state probabilities by habitat covariates: (a) distance to shore, (b) bottom depth, (c) distance to bottom, (d) fish depth. Solid = day, dashed = night.

**Table 3. CLMM results (distance to shore example).**

| Response | Predictors                        | Odds Ratio | 95% CI       | p-value |
|----------|-----------------------------------|------------|--------------|---------|
| Activity (Inactive \| Low Active) | dist_shore_scaled              | 0.86       | 0.63–1.18 | 0.348   |
| Activity (Low Active \| High Active) | diel period [night]            | 7.89       | 5.77–10.77| <0.001  |
| …        | …                                 | …          | …            | …       |

*(full model outputs in Appendix)*

---

### Three-dimensional habitat use
- Day: inactivity nearshore/shallow, high activity spanned deep zones.  
- Night: states shifted shallower, overlaps increased.  

![Figure 13](figures/Figure%2013.png)  
**Figure 13.** Spatial overlay of 50% (core) and 95% (home range) kernel distributions by diel period. Dashed line = thermocline; shaded area = reservoir bottom.

**Table 4. Maximum spatial extent per state × diel period.**

| State        | CI  | Diel | Max Fish Depth (m) | Max Reservoir Depth (m) |
|--------------|-----|------|---------------------|--------------------------|
| Inactivity   | 50  | Day  | 4.05 | 5.92 |
| Inactivity   | 50  | Night| 3.64 | 4.20 |
| Low Activity | 95  | Day  | 15.08| 36.30|
| High Activity| 95  | Night| 7.51 | 36.42|

**Table 5. Percent overlap of 3D areas (day vs night).**

| Ref. State | Overlap State | Day Core (%) | Night Core (%) |
|------------|---------------|--------------|----------------|
| Inactivity | High Activity | 50.9 | 58.3 |
| Inactivity | Low Activity  | 100  | 83.2 |
| Low Activity| High Activity| 46.0 | 40.8 |

---

## Discussion
- Catfish displayed three consistent activity states, with inactivity dominating daytime and high activity at night.  
- Temperature modified diel cycles, extending or advancing activity transitions.  
- Habitat strongly structured behaviour: nearshore/littoral zones linked with inactivity, deeper pelagic zones with activity.  
- Spatial overlap between states suggests flexible use of reservoir habitats, with management implications for predator–prey interactions in lentic systems.  

---

## Appendix

![Figure A1.1](figures/Figure%20A1.1.png)  
**Figure A1.1.** Hourly mean epilimnion temperature. Red dashed line = global average.

![Figure A1.2](figures/Figure%20A1.2.png)  
**Figure A1.2.** Diurnal change in epilimnion temperature. Dark grey = twilight; light grey = sunrise/sunset.

---

**Table A1.1. Initial HMM parameters**

| Parameter | State 1 | State 2 | State 3 |
|-----------|---------|---------|---------|
| Step Mean | 3       | 25      | 75      |
| Step SD   | 3       | 25      | 75      |
| Turning Angle Mean | π | 0 | π/2 |

**Table A1.2. HMM model selection (AIC).**

| Model                        | AIC     | LogLik   |
|-------------------------------|---------|----------|
| No covariates                 | 916760  | -458459  |
| Cos+Sin*time × temp (best)    | 915517  | -457709  |

---

- **Table 2**: CLMM coefficients for habitat × diel period  
- **Table 3**: Summary of overlap metrics  
- **Table 4**: Depth characteristics of 3D habitat use  

---
