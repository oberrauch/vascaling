# Linearity of the volume/area scaling model

The equilibrium run under constant and random climate scenarios produced highly symmetrical results for a positive/negative mass balance (or better for a negative/positive temperature bias of 0.5°C). The suggested explanation is the linear nature of the volume/area scaling (VAS) model or rather of the underlying equations... Hereafter I take a closer look at the individual equations to confirm the stated hypothesis. This includes experiment on the mass balance model as well as on the evolution model, using the same equilibrium and biased climate scenarios as for the equilibrium runs.

## Governing equations



## Mass balance model

### Constant climate scenario

The `ConstantMassBalance` model simulates an average constant climate based on the observations over a 31-year period centred around a given year `y0`​. Consequentially, a `ConstantMassBalance` model with `y0` = $t^*$ keeps the glacier in equilibrium as a result from the mass balance calibration. Hence, adding a negative/positive temperature bias will result in a positive/negative glacier mass balance.

![](../plots/linearity/constant_mb.pdf)

As during the equilibrium experiments a temperature bias of ±0.5°C relative to the equilibrium climate scenario is used. The specific mass balance changes by +196 mm w.e. year$^{-1}$ and -204 mm w.e. year$^{-1}$ for the negative and positive temperature bias, respectively. That is a difference in 4%, so even though the result looks/seems symmetrical it is not.

- [ ] TODO: compare changes to the flowline mass balance model.

### Sensitivity to temperature & precipitation bias

#### Mass balance as function of temperature bias

Let's state the obvious: decreasing/increasing the air temperature will increase/decrease the mass balance. Hence, a negative/positive temperature bias applied to the equilibrium climate scenario will result in a positive/negative mass balance. Hereafter I'll try to quantify said effect.

![](/Users/oberrauch/work/master/plots/linearity/lin_reg_detail.pdf)

For a sensible range of temperature biases from -5°C to +5°C, the constant mass balance function reacts as predicted. The result is almost symmetrical even though slightly positively skewed and can be approximated by a linear function to an acceptable degree of accuracy. For temperature biases $\beta_{T}$ between -1°C and +1°C the function $b(\beta_T) \simeq -400 \beta_T - 5.9$ approximates the actual mass balance function $b$ with a coefficient of determination of 0.9995.

Hence, the symmetry hypothesis is sustainable for small temperature biases. But what about the bigger picture?!

![](../plots/linearity/lin_reg_temp.pdf)

Let's exploit the two extreme cases:

1. Once the terminus temperature is below the threshold of solid precipitation as well as below the temperature threshold for ice melt, the specific mass balance will be constant. Even lower temperatures will not increase the amount of solid precipitation, neither will the melting rate decrease below zero. The mass balance equation reduces to
   $$
   B = \sum_{i=1}^{12} P_{i}^{\text{solid}} =\sum_{i=1}^{12} P_{i}^{\text{}}.
   $$
   All precipitation falls solid, hence $P_{i}^{\text{solid}} = P_{i}^{\text{}}$.

2. Once the air temperature at the maximum glacier surface elevation is above the threshold for solid precipitation, the glacier wont get any mass input in form of solid precipitation. Increasing the temperature further will increase the melt rate, explaining the linear decrease of specific mass balance. The mass balance equation reduces to
   $$
   B = \sum_{i=1}^{12}\left[
   			- \mu^* \cdot
   					\left(T_{i}^{\text{terminus}}
   					- T_{\text{melt}}\right)
   	\right].
   $$
   
   The slope is given by the temperature sensitivity $\mu^*$.

Since the air temperature affects the melting rate and the amount of solid precipitation, the transition between those two extreme cases is not linear.

#### Mass balance as function of precipitation bias

While the input air temperature is used to calculate different mass balance variables (terminus temperature, solid precipitation amount), the precipitation bias is just a multiplication factor for the precipitation input (analog to the precipitation scaling factor). Hence, the mass balance as function of the precipitation bias is linear.

![](../plots/linearity/lin_reg_prcp.pdf)

## Evolution model

![](../plots/linearity/volume_norm.pdf)

|                            | Positive MB       | Negative MB       |
| :------------------------- | :---------------- | :---------------- |
| Specific mass balance      | +196 mm w.e. year | -204 mm w.e. year |
| Initial volume             | 0.596 km          | 0.596 km          |
| Min/Max volume (overshoot) | 0.729 km (122%)   | 0.480 km (81%)    |
| Final/Equilibrium volume   | 0.697 km (117%)   | 0.502 km (84%)    |

## Discussion



## Conclusion

