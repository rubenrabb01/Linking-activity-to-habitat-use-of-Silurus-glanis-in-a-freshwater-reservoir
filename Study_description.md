# Catfish Movement in Římov Reservoir

> **Graphical Abstract (top)**  
> *Graphic depiction of the expected turning angle (radians) and step length (m) distributions for a three-state Hidden Markov Model.*  

<p align="center">
  <img src="figures/graphical_abstract.png" alt="Graphical Abstract" width="80%"/>
</p>

---

## Table of Contents

- [Overview](#overview)
- [Study Area](#study-area)
- [Data \& Telemetry](#data--telemetry)
- [Environmental Data](#environmental-data)
- [Analysis Workflow (End-to-End)](#analysis-workflow-end-to-end)
  - [0. Environment Setup](#0-environment-setup)
  - [1. Load Libraries \& Helper Functions](#1-load-libraries--helper-functions)
  - [2. Load \& Filter Telemetry Data](#2-load--filter-telemetry-data)
  - [3. Clean \& Prepare Tracks for HMM](#3-clean--prepare-tracks-for-hmm)
  - [4. Build \& Select HMMs](#4-build--select-hmms)
  - [5. Decode States \& Predict Stationary Probabilities](#5-decode-states--predict-stationary-probabilities)
  - [6. Habitat Modeling (CLMMs)](#6-habitat-modeling-clmms)
  - [7. 3D Habitat Use (rKIN)](#7-3d-habitat-use-rkin)
- [Results](#results)
  - [Activity States](#activity-states)
  - [Time Allocation \& Spatial Distribution](#time-allocation--spatial-distribution)
  - [Influence of Temperature \& Time of Day](#influence-of-temperature--time-of-day)
  - [Habitat Use and Diel Period (CLMMs)](#habitat-use-and-diel-period-clmms)
  - [Three-Dimensional Habitat Use](#three-dimensional-habitat-use)
- [Figures](#figures)
- [Tables](#tables)
- [Discussion (brief)](#discussion-brief)
- [Appendix I](#appendix-i)
- [Appendix II](#appendix-ii)

---

## Overview

We analyzed the summer movement behavior of European catfish (*Silurus glanis*) in the Římov Reservoir (Czech Republic) using high-resolution acoustic telemetry and a three-state Hidden Markov Model (HMM). We integrated temperature (epilimnion) and diel cycling, mapped state-specific spatial use, modeled habitat effects with cumulative link mixed models (CLMMs), and visualized 3D habitat overlaps (core vs home-range). The full computational pipeline and figures are embedded below.

[Back to top](#table-of-contents)

---

## Study Area

Římov Reservoir is a canyon-shaped reservoir in South Bohemia, Czech Republic (N 48°51.00978′, E 14°29.47462′), built in 1978 for drinking water and flood control. The main reservoir body extends ~10 km; surface area 210 ha; mean depth ~16 m (max 45 m). A longitudinal gradient exists from mesotrophic conditions at the dam to eutrophic conditions toward the main inflow (Šeda & Devetter, 2000), with parallel gradients in algae, zooplankton, and fish density. The system typically stratifies from April–September (Říha et al., 2022). Due to fluctuating water levels and steep banks, submerged aquatic macrophytes are typically absent in the littoral zone (Říha et al., 2015).

Public access is restricted, with an estimated ~200 catfish and minimal angling disturbance, making the system well-suited for studying natural movement in a lentic environment. During the study period (July–August 2017), epilimnion temperatures ranged 19.7–22.7 °C (mean 21 ± 0.6 °C). Thermocline depth averaged 5.47 m.

[Back to top](#table-of-contents)

---

## Data & Telemetry

An array of 91 acoustic receivers (Lotek WHS3250) was deployed on 2017-04-18 to cover the entire reservoir (87), the inflowing tributary (3), and one small bay near the dam (1). This design allowed fine-scale high-resolution positions (see Figure 1). Temperature stratification was monitored using HOBO Pendant temp/light 64K loggers at four longitudinal locations (loggers placed every 1 m from surface to 13 m; plus 20 m at dam and mid sections). For analysis, we used the mean 0–6 m temperature (epilimnion).

Fifteen *S. glanis* were captured by long-lining (per Vejřík et al., 2017) at four reservoir locations (2017-04-18 to 2017-04-25), anesthetized with 2-phenoxy-ethanol (0.7 ml l⁻¹), surgically implanted with Lotek MM-M-11-28-TP transmitters (65×12 mm, 13 g; burst 15 s; plus pressure/temperature sensors), and released at capture sites. Locations were post-processed with UMAP v1.4.3 (Lotek); data cleaning follows Říha et al. (2021, 2022). For this study, we used 2017-07-01 to 2017-08-31 (peak activity/growth; Copp et al., 2009). After filtering, 13 individuals and 85,650 positions remained.

[Back to top](#table-of-contents)

---

## Environmental Data

- Epilimnion temperatures: mean of 0–6 m across four stations.  
- Diurnal/cyclic time derived from timestamps (UTC) via trigonometric transforms.  
- Thermocline depth for 3D visualization context.

[Back to top](#table-of-contents)

---

## Analysis Workflow (End-to-End)

> **Tip:** Each section below includes the exact code used in the analysis. 

### 0. Environment Setup

- **Project structure**
