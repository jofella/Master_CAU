/*---------------------------------------------------------------------------

Group Member:
- Josef Fella / Quantiative Finance (enrolement?) / stu245231
- Robert Hennings / Quantiative Finance (enrolement?) / ...
- ...

---------------------------------------------------------------------------*/
clear


* ============== Excercise 1) ===============================================

**** a) Import and briefly describe the dataset:

cd "C:\Users\josef\Documents\GitHub\Master_CAU\Semester_1\Eco_1"
use HA_smoking.dta

* Get sense of data
describe


**** b)  Regress with OLS with heteroscedasticity robust SD. ^Beta^ significant?

* New variable age^2
gen agesq = age * age

* OLS regression
regress smoker smkban age agesq female hsdrop hsgrad colsome colgrad, vce(robust)


/* Interpretation:

- If we look purely on the p-value of 0.000 we would argue that its statistically significant, since its below the significance level of 0.05. So we have statistical evidence to reject our H0 (less smoking with smoking ban).
- On the other hand we the coefficient is -.048052, so 1 increased unit of smoking ban would decrease the predicted probability of being a smoker by .048052. In absolute terms this is quiet small.

--> Statistically relevant BUT practical importance will be dependent of the study and scale of variable involved.

 */ 

 

**** c) Wald test

test hsdrop = hsgrad
* Solution of sheet 4 d)



**** d)

// Assuming these are your coefficients
local alpha   -0.0268721
local beta1   -0.048052
local gamma_age    0.0097597
local gamma_agesq  -0.0001291
local gamma_female -0.032133
local gamma_hsgrad 0.2245546
local gamma_colgrad 0.0433244

// Given values
local smkban   1
local age      70
local agesq    age * age
local female   1
local hsgrad   1
local colgrad  1

// Calculate the predicted linear probability
local lin_pred $alpha + $beta1 * `smkban' + $gamma_age * `age' + $gamma_agesq * `agesq' + ///
                      $gamma_female * `female' + $gamma_hsgrad * `hsgrad' + $gamma_colgrad * `colgrad'

// Calculate the predicted probability using the logistic function
gen predicted_prob = exp(`lin_pred') / (1 + exp(`lin_pred'))

// Display the predicted probability
display "Predicted Probability of Smoking: " predicted_prob







* ============== Excercise 2) ===============================================

* a)

* b)
