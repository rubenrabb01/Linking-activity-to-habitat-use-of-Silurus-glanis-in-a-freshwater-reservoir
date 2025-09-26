# Linking activity to habitat use of *Silurus glanis* in a freshwater reservoir

## Table of Contents
1. [Methods](#methods)
   - [Study area](#study-area)
   - [Data collection](#data-collection)
   - [Statistical analysis](#statistical-analysis)
2. [Results](#results)
   - [Activity states](#activity-states)
   - [Time allocation and spatial distribution](#time-allocation-and-spatial-distribution)
   - [Influence of temperature and time of day](#influence-of-temperature-and-time-of-day)
   - [Habitat use and diel period](#habitat-use-and-diel-period)
   - [Three-dimensional habitat use](#three-dimensional-habitat-use)
3. [Figures](#figures)
4. [Tables](#tables)

---

## Methods

### Study area
Rimov Reservoir is a canyon-shaped reservoir in South Bohemia, Czech Republic (N 48°51.00978′, E 14°29.47462′), built in 1978 for drinking water storage and flood control. The main reservoir body extends 10 km, covers 210 ha, and reaches 45 m maximum depth (Figure 1). It exhibits mesotrophic–eutrophic gradients (Šeda & Devetter, 2000), strong summer stratification (Říha et al. 2022), and fluctuating water levels that prevent macrophyte growth in the littoral zone (Říha et al. 2015).  
The reservoir hosts ~200 catfish with limited human disturbance, making it an ideal lentic system to study predator activity.  
During the study, epilimnion temperature ranged 19.7–22.7 °C (mean 21.0 ± 0.6 °C; Figure S1).

### Data collection
An array of 91 Lotek WHS3250 receivers was deployed on 18 April 2017 across the reservoir, tributary, and bay (Figure 1). Temperature was monitored with HOBO loggers at 1-m intervals (0–13 m plus 20 m).  
Fifteen *S. glanis* were captured, anaesthetized (2-phenoxy-ethanol, 0.7 ml l−1), tagged with Lotek MM-M-11-28-TP transmitters, and released between 18–25 April 2017. After filtering, 13 individuals were retained for analysis.

### Statistical analysis

#### Data processing
- Positions calculated with Lotek UMAP v1.4.3, averaged to 5-min intervals.  
- Period selected: 1 July–31 August 2017 (summer peak).  
- Habitat variables: fish depth, distance to shore (via `rgeos`), reservoir depth (1 × 1 m bathymetry), distance to bottom.  
- Missing values interpolated if <1.5 h gaps; otherwise excluded.  

#### Identifying activity states
- HMM fitted with `moveHMM` (Michelot et al. 2016).  
- Step length = distance between consecutive locations.  
- Turning angle = directional change between steps.  
- Missing data interpolated with CTCRW (`momentuHMM::crawlWrap`) for gaps ≤1.5 h.  
- Three states (inactive, low, high activity) chosen as ecologically appropriate.  
- Covariates: epilimnion temperature (standardized) + time of day (sin/cos).  
- Five candidate models compared (no covariates, temp, time, additive, interaction).  
- Best model: interaction between cyclic time and temperature (AIC).  
- Most likely sequence estimated with `viterbi()`.  

#### Mapping activity states to home range
- Kernel UD (95% contours) via `adehabitat` (Calenge 2006).  
- Overlap and areas calculated, compared with ANOVA.  

#### Predicting habitat use
- CLMM (`ordinal`) with activity state as ordered response (inactive < low < high).  
- Fish ID random intercept.  
- Separate models for each habitat variable (distance to shore, depth, distance to bottom, absolute fish depth).  
- Diel period treated as effect modifier.  
- Model selection by AIC.  

#### Diel differences in habitat use
- Kernel isopleths (50% = core, 95% = home range) via `rKin`.  
- Polygons compared by diel period and plotted with depth profiles.  

---

## Results

### Activity states
Three states identified:  
- **High activity**: long steps (72.1 ± 29.5 m), small angles (κ = 1.7).  
- **Low activity**: medium steps (21.3 ± 15.3 m), small angles (κ = 0.6).  
- **Inactive**: short steps (2.4 ± 2.2 m), large angles (κ = 0.2).  

Areas differed significantly (p < 0.001): inactive = 8.0 ha, low = 37.8 ha, high = 74.5 ha. High activity overlapped lower states (Table 1).

### Time allocation and spatial distribution
Most individuals spent 30–48% of time inactive (Figure 5), often clustered in concentrated reservoir areas, but locations varied among fish.

### Influence of temperature and time of day
At mean temp (21 °C):  
- High activity peaked sunset–night (16:30–01:30).  
- Inactivity dominated daytime (02:00–15:30).  
- Low activity dominated transitions.  
At high temp (22.7 °C): activity extended later into night.  
At low temp (19.7 °C): transitions earlier, more low activity during day (Figure 6).  

### Habitat use and diel period
- All CLMMs: diel period strong predictor (p < 0.001).  
- Distance to shore, bottom depth, and distance from bottom strongly increased activity (OR ~5–6, all p < 0.001).  
- Effects weaker at night.  
- Absolute fish depth weakest predictor (lowest AIC).  

### Three-dimensional habitat use
Day: inactivity nearshore/shallow, high activity spanning deep zones, low activity intermediate (Figure 8, Table 4).  
Night: all states shifted shallower; overlaps increased.  

---

## Figures
- **Figure 1**: Map of Rimov reservoir and receiver array  
- **Figure 2**: Example HMM state classification  
- **Figure 3**: Step length and turning angle distributions per state  
- **Figure 4**: State home range areas (boxplot)  
- **Figure 5**: Time allocation per individual  
- **Figure 6**: Predicted state probabilities vs diel cycle × temperature  
- **Figure 7**: Habitat CLMM effects  
- **Figure 8**: Day/night 3D core areas and home ranges  

---

## Tables
- **Table 1**: Area and overlap of activity states  
- **Table 2**: CLMM coefficients for habitat × diel period  
- **Table 3**: Summary of overlap metrics  
- **Table 4**: Depth characteristics of 3D habitat use  

---
