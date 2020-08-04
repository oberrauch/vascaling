# Equilibrium experiments

**Introduction, or 'What am I doing?! And why?!'**

Equilibrium runs are a useful tool to asses the behavior of glacier models. The OGGM provides two convenient mass balance models (or better climate scenarios) for equilibrium experiment, the `ConstantMassBalance` model and the `RandomMassBalance` model. So first of all, those two mass balance model are implemented to work with the volume area scaling model.

All experiments are performed on the Hintereisferner (RGI60-11.00897) in the Ötztal Alp in Austria. The first qualitative conclusions are drawn from the evolution of the different glacier geometries (length, area, volume) under different climatic forcings. Scaling methods applied to a single glacier give only an order of magnitude estimation (cf. [Bahr et al. (2015)][], section 8.5), which is accounted in the following analysis.

More quantitative results are drawn from an auto correlation analysis and a power spectral density analysis, inspired by [Roe, Baker (2014)][Roe, Baker (2014)].

**Research question, or What do I expect?!:**

How do the models differ and why?

- [ ] Best case scenario: the two models produces comparable results, corresponding with expectations/observations/knowledge... Therefore I would not have to argue/find explanations. Very unlikely, though...
- [ ] Most likely scenario: The model results differ in a coherent and explainable fashion. Only thing to do is to find said explanations.
- [ ] Worst case scenario: There is no discernible pattern in the model behavior and no explanations can be found.

## Mass balance models

**Temperature index mass balance model:**

In a nutshell, a glaciers (annual specific surface) mass balance $B$ is the difference between accumulation (i.e., mass gain by snowfall, avalanches, drift, ...) and ablation (i.e., mass loss via ice melt, sublimation, ...) over the course of a year. The used temperature index mass balance model relies solely on the monthly solid precipitation onto the glacier surface $P_i^{\text{solid}}$ and the monthly mean air temperature $T_i^{\text{}}$ as input. The details/differences between VAS model and flowline model will be explained later. 
$$
B = \left[\sum_{i=1}^{12}\left[
			P_{i}^{\text{solid}} 
			- \mu^* \cdot
					\max\left(T_{i}^{\text{}}
					- T_{\text{melt}},\ 0\right)
	\right]\right] - \beta^*
$$
Hereby, $\mu^*$ is the glacier specific temperature sensitivity, $\beta^*$ is the mass balance residual (compared to observations) and $T_{\text{melt}}$ the mean air temperature above which ice melt occurs. So per definition, $\mu^*$ is the temperature sensitivity to keep the glacier in equilibrium over the 31-year climate period around the year $t^*$ while neglecting a potential mass balance residual $\beta^*$.

**Differences between the flowline model and the volume/area scaling model:**

The flowline model requires a (point) mass balance value for every grid point of the flowline (i.e., for each elevation band). Therefore, the mass balance is a function of elevation and the heights of the grid points must be supplied. Solid precipitation and air temperature are computed on the given heights.

The volume/area scaling model on the other hand computes an average value for the entire glacier. The mass balance model requires only the minimal and maximal glacier elevation as input parameters, to compute the terminus temperature and the area averaged amount of solid precipitation.

As mentioned before, there are two suitable climate scenarios to perform equilibrium experiments. Both mass balance models simulating said scenarios are implemented to work with both evolution models.

### Constant climate scenario

The `ConstantMassBalance` model simulates an average constant climate based on observations over a 31-year period centered on a given year `y0`​. Hence, the specific mass balance does not change from year to year. The task `run_constant_climate` takes an additional temperature bias (and a possible precipitation bias, for that matter) as parameters, to alter the observed climatic conditions.

The same idea of a constant climate is used during the mass balance calibration, solving the mass balance equation for the temperature sensitivity $\mu^*$. The temperature sensitivity $\mu^*$ is calibrated, in order to keep the model glacier in equilibrium under a constant climate centered on $t^*$.
$$
\mu^* = 
	\frac{P_{\text{clim. avg}}^{\text{solid}}}
		{\max(T_{\text{clim. avg}}^{\text{}} - T_{\text{melt},}\ 0)},
$$
whereby $P_{\text{clim. avg}}^{\text{solid}}$ and $T_{\text{clim. avg}}^{\text{}}$ are the average yearly solid precipitation amount and average yearly air temperature during the climatological period centered on $t^*$, respectively. Consequentially, a `ConstantMassBalance` model with `y0` = $t^*$ keeps the glacier in equilibrium.

### Random climate scenario

The `RandomMassBalance` model takes the climate information from a randomly chosen year within a 31-year period centered on a given year `y0` to compute the specific mass balance. Hence, the model runs on a synthetic (random) climate scenario based on actual observations.

Using the climatological period around the "equilibrium" year $t^*$, the model glacier should stay in an equilibrium state (underlying minor fluctuations). In analogy to the `ConstantMassBalance` model, the `run_random_climate` task allows to add a temperature bias. Increasing/decreasing the temperature of the equilibrium period should result in a retreating/advancing model glacier, reaching a new equilibrium after some years.

## Methods, experimental setup

All experiments are performed on the Hintereisferner (RGI60-11.00897, Ötztal, Austria) using the HistAlp climate data [Auer et al. (2007)][] with the corresponding hyper parameters (see [Mass-balance model calibration for the Alps](https://oggm.org/2018/08/10/histalp-parameters/) on the OGGM blog for more information).

The needed preprocessing includes GIS tasks (computing a local grid using the SRTM DEM and the RGI outline, computing centerlines), climate tasks (preparing the HistAlp data), mass balance calibration (computing the temperature sensitivity $\mu^*$) as well as the inversion tasks (estimating a bed topography) for the flowline model.

As explained above, the mass balance model calibration depends on the chosen "equilibrium" year $t^*$. Hence, $t^*$ must be equal for both evolution models since we want to use the same climatic period. Rather than relying on the reference tables, which are different for the different evolution models, the temperature sensitivity $\mu^*$ is computed using $t^*$ =  1927​ as equilibrium year and no mass balance residual ($\beta^*$ = 0). This corresponds to the reference year for the flowline model and is not too far off from the reference year for the VAS model (which is 1905). 

Both evolution models run for 3'000 years with the `ConstantMassBalance` model and for 10'000 years with the `RandomMassBalance` model, both initialized around $t^*$. Furthermore, each climate scenario runs with three different temperature biases of 0°C, +0.5°C and -0.5°C resulting in an equilibrium run, a negative and a positive mass balance run, respectively.

The yearly geometric properties (length, area and volume) of the model glacier are stored to allow further investigations. In addition to the absolute values, a dataset with normalized values (with respect to the initial value) is produced, allowing better comparability.

## Results

### Time series

The following section tries to explain the model behavior. The plots show a comparison between VAS and flowline model time series (length and volume), both for the constant and random climate scenario. The results are normalized with respect to the initial value for better comparability.

**Overall findings:**

1. Both model behave as expected and produce the same qualitative results. The glacier stays in an equilibrium state using the climate around $t^*$ and decreases/increases in size for a positive/negative temperature bias. Plots with absolute values can be found in the appendix. (@TODO)

2. The glacier size changes (dramatically) less under the VAS model than under the flowline model (true for length, area, and volume).

   *Note*: However, volume estimations from volume/area scaling of a single glaciers must be considered as order of magnitude result. The scaling constant $c$ is a random variable which varies (drastically) from glacier to glacier. Apparently, the global mean value of $c=0.034\ \mathrm{km^{3-2\gamma}}$ is a bad fit for the characteristics of Hintereisferner.

   *Second Note*: Changing the scaling constants changes the absolute values of ice volume (as well as surface area and glacier length). A higher volume/area scaling constant results in a larger initial ice volume. Subjected to the same climate perturbation (temperature step change), an initially larger glacier will gain/loose more ice and reach a higher equilibrium ice volume than a smaller one. However, when normalized with initial ice volume there are no more discernible differences in the magnitude of ice volume change. The temporal evolution, i.e., the oscillation behavior, is comparable, even if smaller glaciers react faster than larger ones (which is to be expected).

   *TL;DR; Turns out, the scaling constant does not change the magnitude of the normalized volume change.* 

3. The glacier length of the VAS model has to be seen more as a model parameter, rather than as an actual glacier property. The VAS glacier length decreases/increases only by about 8 percent compared to its initial value, for a positive/negative temperature bias of 0.5 °C. This correspond to an absolute length change of less than 400 m, which is very little compared to the 3 to 4 km in length change (~40% of the initial value) produced by the flowline model. (*Note*: May change with different $c$ parameter. *Second note*: Turns out, it does not when normalized with initial length.)

4. The result of VAS model under a constant climate scenario with a non-zero temperature bias reminds of a damped oscillating signal. The modeled length reaches its maximum after ~200 years, overshooting the equilibrium result by more than 1%. Followed by two minor but still discernible peaks until the new equilibrium is reached. Both, glacier surface area and glacier volume reach their maximum earlier and overshoot by more.

#### Constant climate scenario

**Equilibrium climate:** While the VAS model glacier stays absolutely the same under equilibrium climate, the flowline model seems not be in equilibrium from the get go. The volume of the flowline glacier increases by about 8% over the first 750 years, before it reaches a new equilibrium. The volume change is accompanied by a corresponding change in glacier surface area and glacier length, however both are showing an initial negative swerve.

**Positive and negative mass balance:** ~~The VAS model produces almost perfectly symmetric results for the positive and negative mass balance scenario, compared to the equilibrium run. Given the linear nature of the underlying equations, this behavior is expected.~~ Both VAS runs show a damped oscillating behavior, with an initial positive/negative overshoot. Compared to the flowline model the changes in length, area and volume are much smaller, ranging between 10% and 30% with respect to the initial value. For the flowline model, a negative temperature bias of 0.5°C results in a volume increase of about 74% (from 0.76 km$^3$ to 1.33 km$^3$), while a positive temperature bias of 0.5°C results in a volume decrease of about 36% (from 0.76 km$^3$ to 0.49 km$^3$). Both models react stronger to the positive mass balance (negative temperature bias) than to the negative mass balance (positive temperature bias), while it is more distinct for the flowline model. Given that the OGGM considers the actual glacier geometry, a stepwise mass balance change triggers a positive feedback loop. A positive mass balance perturbation results in a greater accumulation area, which itself increases the mass balance furthermore, explaining the asymmetrical behavior.

![](../plots/eq_runs/volume_norm_comparison_constant.pdf)

![](../plots/eq_runs/length_norm_comparison_constant.pdf)

#### Random climate scenario

The model behavior under the random climate scenario is essentially the same as under the constant climate scenario, just with noisier results. The glacier geometries of the VAS model and the flowline model under the same temperature bias are highly correlated. The glacier lengths shows the highest correlation, whereas the areas shows the lowest. Overall the correlation coefficient ranges between 0.57 and 0.81, for the different glacier geometries under different climate scenarios.

![](../plots/eq_runs/length_norm_comparison_random.pdf)

![](../plots/eq_runs/volume_norm_comparison_random.pdf)

### Time series analysis

All analysis in the following paragraphs are based on the modeled length changes under the random climate scenarios with different temperature biases.

#### Auto correlation function

Overall, the  VAS model shows less auto correlation at shorter lag times (< 200 years) than the flowline model, indicating more high frequency variability. The auto correlation functions of the different models converge after about 500 years

The auto correlation functions of the VAS model length are quite similar for the different temperature biases (correlation between 0.99 and 0.96), and virtually equal for the first 150 years of lag time (RMSD below 0.03). This is again explained by the underlying linear equations and the resulting symmetry of the different runs. Additionally, the auto correlation function extracts the same oscillating behaviour seen under the constant climate scenario from the noisier signal produced by the random climate scenario.

The flowline model run with the negative mass balance has a similar auto correlation function to the VAS model for the first 65 years. This can probably be explained by a faster response time to the positive temperature bias, compared to the other two flowline runs. Thereafter, the auto correlation it is again higher than for the VAS model runs, but still somewhat different from the other two flowline model runs.

![](../plots/eq_runs/acf_length_random.pdf)

#### Power spectral density

As a general observation, the power of all length change signals decreases with increasing frequency. This means, that the change in glacier length (and therefore also change in surface area and change in glacier volume) is driven mainly by long term trends rather than by the inter-annual variability of the climate forcing. As seen before, the flowline model glacier changes (much) more than the VAS model glacier, reflected in the PSD analysis with an overall higher power.

The PSD of the VAS model is almost identical between runs with different temperature biases, attributed to the linear (and symmetric) behaviour seen before. The flowline model produces a less coherent PSD between the runs with different temperature bias, whereby the equilibrium run and the positive mass balance run are more alike.

**Flowline model:** The flowline model run with a positive mass balance is a textbook example for a low pass filter. A low pass filter passes lower frequencies while attenuation (or rejecting) higher frequencies. And so does the glacier (model). After an initial decrease in power density with increasing frequencies, the power density stays fairly constant at around 4.5$\cdot10^4$ (@TODO: unit) for frequencies above 0.1 year$^{-1}$. The longterm changes of the climate system are driving the glacial evolution, while the short term changes in monthly or yearly climatic parameters have little to no effect. The PSD suggests that trends  shorter than 10 years (corresponding to a frequency of 0.1 year$^{-1}$) are filtered by the Hintereisferner.

The spectral analysis of the equilibrium run and the negative mass balance run are quite similar. While both share the overall behaviour of a low pass filter, there is a strong increase in power density for frequencies higher than 0.4 year$^{-1}$ (trends of 2.5 years or less). This is explained by... what exactly?! @TODO.

**VAS model:** The PSD of the different VAS model runs show no discernible differences between the different temperature biases. The power density $P$ decreases constantly on the logarithmic scale. The energy is transferred from larger scales to smaller scales (direct energy cascade), in contrast to the low pass behavior of the flowline model. A linear regression analysis yields the power relationship $P \propto f^{-3.5}$ with a coefficient of determination r$^2$ = 0.98, for frequencies $f$ between 0.01 and 0.5 year$^{-1}$.



![](../plots/eq_runs/psd_length_random.pdf)

## Discussion

@TODO: explain/discuss the difference in behaviour between the models.

### Open questions:

**Time series:**

- [ ] Underestimation of the overall/longterm length/area/volume change by the VAS model, but more high frequency variability?!
- [ ] Explain damped oscillating behavior...!

**Auto correlation function:**

- [ ] Does the time scale makes sense?!
- [ ] Explain behaviour of flowline model run with positive mass balance bias.

**Power spectral density**

- [ ] Unit of power density (y-axis)?!
- [ ] VAS model: what does the power relation tell me?!
- [ ] Flowline model: positive mass balance vs. equilibrium and negative mass balance.



## References

[Roe, Baker (2014)]: https://doi.org/10.3189/2014jog14j016	"Glacier response to climate perturbations: an accurate linear geometric model"
[Bahr et al. (2015)]: http://doi.org/10.1002/2014RG000470	"A review of volume–area scaling scaling of glaciers"
[Auer et al. (2007)]: http://www.zamg.ac.at/histalp/index.php	"Auer I, Böhm R, Jurkovic A, Lipa W, Orlik A, Potzmann R, Schöner W, Ungersböck M, Matulla C, Briffa K, Jones PD, Efthymiadis D, Brunetti M, Nanni T, Maugeri M, Mercalli L, Mestre O, Moisselin J-M, Begert M, Müller-Westermeier G, Kveton V, Bochnicek O, Stastny P, Lapin M, Szalai S, Szentimrey T, Cegnar T, Dolinar M, Gajic-Capka M, Zaninovic K, Majstorovic Z, Nieplova E, 2007. HISTALP – Historical instrumental climatological surface time series of the greater Alpine region 1760-2003. International Journal of Climatology 27: 17-46"

