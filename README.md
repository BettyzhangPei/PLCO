#  APDI: Assessing the Performance of Dietary Assessment Instruments
Food Frequency Questionnaires (FFQs) are extensively utilized in nutritional epidemiology to gauge dietary intake over specified periods. However, FFQs often encounter measurement errors. 

Given the availability of a reference instrument (F), APDI evaluates the performance of dietary assessment  instrument (Q) like FFQs using two common metrics:
1. The correlation (ρ) between true intake (T) and the assessment instrument (Q).
2. The attenuation factor (λ).

Building on the measurement error model by Thompson et al. (2008), our R package APDI employs the method of moments for estimation. This package estimates:
- The correlation between nutrient intake (Q_N) and true nutrient intake (T_N) for each nutrient (N).
- The correlation between energy intake (Q_E) and true energy intake (T_E).   
- The correlation between nutrient residuals (Q_R) and true residuals (T_R).
- The attenuation factors for Q_N, Q_E, and Q_R.

# Installation
The current version can be installed from source using the package `devtools`
```r
devtools::install_github("BettyzhangPei/AARP")
```

# Usage Examples
- `AI-ReML` function

The `NEM` function provides estimates based on various types of input nutrients. For instance:
```r
NME(Q.N.1= FFQ.Protein.1, Q.N.2= FFQ.Protein.2,
 Q.E.1= FFQ.Energy.1, Q.E.2= FFQ.Energy.2,
 F.N.1= Ref.Protein.1, F.N.2= Ref.Protein.2,
 F.E.1= Ref.Energy.1, F.E.2= Ref.Energy.2,
 data=input_data)
```
where
- Q.N.1: An assessment instrument for nutrient (N) at time 1.
- Q.N.2: An assessment instrument for nutrient (N) at time 2.
- Q.E.1: An assessment instrument for total energy (E) at time 1.
- Q.E.2: An assessment instrument for total energy (E) at time 2.
- F.N.1: A reference instrument for nutrient (N) at time 1.
- F.N.2: A reference instrument for nutrient (N) at time 2.
- F.E.1: A reference instrument for total energy (E) at time 1.
- F.E.2: A reference instrument for total energy (E) at time 2.

# Usage Notes
1. We recommend transforming each variable (e.g., using a log transformation) to approximate a normal distribution before applying our functions.
2. Our calculations assume there is no missing data in the dataset.
3. Following Thompson et al. (2008), nutrients and energy are measured twice for each Q and F.
   
# References
Thompson, F.E., Kipnis, V., Midthune D., Freedman, L.S., Carroll, R.J., Subar, A.F., Brown, C.C., Butcher, M.S., Mouw, T., Leitzmann, M., Schatzkin, A.(2008) Performance of a food-frequency questionnaire in the US NIH-AARP (National Institutes of Health-American Association of Retired Persons) Diet and Health Study. Public Health Nutrition, 11(2): 183-95.

Carroll, R.J., Midthune, D., Freedman, L. S., and Kipnis, V. (2006) Seemingly unrelated measurement error models,
with application to nutritional epidemiology. Biometrics, 62(1):75–84
   
