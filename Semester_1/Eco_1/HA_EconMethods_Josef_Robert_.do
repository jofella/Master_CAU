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

 
**** c) Interpretation and Wald test:

/* Does probability of smoking increase or decrease? Given the coefficients:
- hsdrop (dropout): Coefficient = 0.2821909
- hsgrad (graduate): Coefficient = 0.2245546

--> Model implies higher probability for hsfropouts than for graduates
*/


* Wald test
test hsdrop = hsgrad

/* Interpretation:
- H0: hsdrop = hsgrad (they are equal)

- Result: We reject the H0 --> F = 10.25 and p-value 0.0014 (less than 0.05)
- So both variables are statistically significantly different
*/


**** d) Predicted probability for the women: (logic is to store the values & coefficients and calc the probability)



regress smoker smkban age agesq female hsdrop hsgrad colsome colgrad, vce(robust)

* Set the values for the woman's characteristics
scalar smkban_value = 1
scalar age_value = 70
scalar agesq_value = age_value^2
scalar female_value = 1
scalar hsdrop_value = 0
scalar hsgrad_value = 1
scalar colsome_value = 0
scalar colgrad_value = 1

* Set the coefficients from your regression output
scalar smkban_coeff = -0.048052
scalar age_coeff = 0.0097597
scalar agesq_coeff = -0.0001291
scalar female_coeff = -0.032133
scalar hsdrop_coeff = 0.2821909
scalar hsgrad_coeff = 0.2245546
scalar colsome_coeff = 0.1556154
scalar colgrad_coeff = 0.0433244

* Calculate the linear prediction
scalar linear_prediction = smkban_coeff * smkban_value + ///
                            age_coeff * age_value + ///
                            agesq_coeff * agesq_value + ///
                            female_coeff * female_value + ///
                            hsdrop_coeff * hsdrop_value + ///
                            hsgrad_coeff * hsgrad_value + ///
                            colsome_coeff * colsome_value + ///
                            colgrad_coeff * colgrad_value

* Calculate the predicted probability (logit)
scalar predicted_probability = exp(linear_prediction) / (1 + exp(linear_prediction))

* According to documentaion of stata
di "Predicted Probability:", %6.4f predicted_probability



/* Interpretation:
- Predicted probability of 0.5593 --> over 50%

Problem:
- We only asume linear relationship of predictor and regressor
- Especially the positvie coefficient of education may suggest a more complex relationship --> not considered here
*/


* ============== Excercise 2) ===============================================

**** a)
ssc install estout


*1) Linear Probability Model:
regress smoker smkban age agesq female hsdrop hsgrad colsome colgrad, vce(robust)

eststo linear_model

*2) Probit Model:
probit smoker smkban age agesq female hsdrop hsgrad colsome colgrad, vce(robust)
margins, dydx(smkban)

eststo logit_model

*3) Logit Model:
logit smoker smkban age agesq female hsdrop hsgrad colsome colgrad, vce(robust)
margins, dydx(smkban)

eststo probit_model

* Display a table with all three estimation results
esttab linear_model logit_model probit_model, se




* b)

* Fit the probit model
probit smoker smkban age agesq female hsdrop hsgrad colsome colgrad, vce(oim)

* Store the probit results
eststo probit_model

* (i) Calculate the effect for a male, 40 years old, college graduate
* Replace the values with those of the specific group
predict smoker_male40_coll, xb
local smkban_value_male40_coll 1
local age_value_male40_coll 40
local agesq_value_male40_coll = `age_value_male40_coll' * `age_value_male40_coll'
local female_value_male40_coll 0
local hsdrop_value_male40_coll 0
local hsgrad_value_male40_coll 1
local colsome_value_male40_coll 1
local colgrad_value_male40_coll 1

* Calculate the average partial effects for the male group
margins, dydx(smkban) at(`smkban_value_male40_coll' `age_value_male40_coll' `agesq_value_male40_coll' ///
                            `female_value_male40_coll' `hsdrop_value_male40_coll' `hsgrad_value_male40_coll' ///
                            `colsome_value_male40_coll' `colgrad_value_male40_coll')


* Display the results
esttab probit_model, margins

* (ii) Calculate the effect for a female, 20 years old, high school dropout
* Replace the values with those of the specific group
predict smoker_female20_drop, xb
local smkban_value_female20_drop 1
local age_value_female20_drop 20
local agesq_value_female20_drop = `age_value_female20_drop' * `age_value_female20_drop'
local female_value_female20_drop 1
local hsdrop_value_female20_drop 1
local hsgrad_value_female20_drop 0
local colsome_value_female20_drop 0
local colgrad_value_female20_drop 0

* Calculate the average partial effects for the female group
margins, dydx(smkban) at(`smkban_value_female20_drop' `age_value_female20_drop' `agesq_value_female20_drop' ///
                            `female_value_female20_drop' `hsdrop_value_female20_drop' `hsgrad_value_female20_drop' ///
                            `colsome_value_female20_drop' `colgrad_value_female20_drop')

* Display the results
esttab probit_model, margins
