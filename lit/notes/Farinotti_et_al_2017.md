### How accurate are estimates of glacier ice thickness? Results from ITMIX, the Ice Thickness Model Intercomparison eXperiment

**Farinotti, Daniel et al (2013) - The Cryosphere - https://doi.org/10.5194/tc-11-949-2017**

#### **Abstract**:

The ITMIX is the first systematic approach in comparing 17 different model and assessing their uncertainties. As was to be expected, individual ice thickness estimations can vary substantially. Averaging the results of different models, however, significantly improves the thickness estimation, with deviations of 10±24% compared to direct measurements. Additionally, models using sets of multiple input data are more sensitive to initial conditions.

#### Conclusions:

> The goal was to assess model performance for cases in which no a priori information about ice thickness is available. The experiment included 15 glaciers and 3 ice caps spread across a range of different climatic regions, as well as 3 synthetically generated test cases.

Experimental setup: estimate ice thickness without prior knowledge; 15 glaciers, 3 ice caps and 3 synthetic test cases representative for a variety of climatic conditions; 17 models using different approaches (minimization, mass conservation, shear-stress, ice flow velocity and others)

> [...] results highlighted the large deviations between individual solutions and even between solutions of the same model category. The local spread often exceeded the local ice thickness. [...] Substantial improvements in terms of accuracy, however, could be achieved when combining the results of different models. Locally, the mean deviation between an average composite solution and the measured ice thickness was on the order of 10 ± 24% of the mean ice thickness (1σ estimate).

Results of individual model (even of the same type) vary significantly, with local spreads larger than the local ice thickness. The average of all results, however, significantly improves the thickness estimation, with deviations of 10±24% compared to direct measurements.

> [...] models that include SMB [surface mass balance], ∂h/∂t, or surface flow velocity fields in addition to the glacier outline and DEM did not perform better when compared to approaches requiring less data, in particular for real-world cases. [...] importance for mutually consistent data sets and suggests that improved observational capabilities could help to improve the performance of the next generation of ice thickness estimation methods.

Additional input data has a detrimental effect on the model accuracy, most likely due to uncertainties and inconsistencies in observational data.

Summary recommendations for future ice thickness estimations:

- ensemble methods produce more robust and reliable results, averaging out random errors in individual models
- observational uncertainties of input data has to be considered (Bayesian approach, see [Brinkerhoff et al. (2016)][])
- new datasets of surface ice flow velocity (e.g., [Scambos et al. (2016)][]) should be considered as input data (with attention to observational uncertainties)
- ice divides are a major source of uncertainties when modeling ice caps and must be handled with care
- a centralized and mutually consistent set of ice thickness measurements (e.g., GlaThiDa) is paramount for model calibration and validation

#### Participation models:

A total of 17 model (from 13 research groups) is used, whereby two modeling approaches are used twice with different/independent implementations. All model can be classified into one of the following categories:

**Minimization approach:** Ice thickness inversion as a minimization problem. A forward model for glacier ice flow (of any type) predicts the surface elevation (or another observable quantity) based on an initial guess for the bedrock topography. Differences in the modeled and observed surface elevation (represented by a cost function) inform the updated model parameters for the next run. This is repeated until the cost function reaches a minimum. Models in this category: Brinkerhoff-v2, [Van Pelt et al. (2013)][], [Fürst et al. (2017)][]

**Mass conserving approaches:** The continuum equation under the assumption of incompressibility states that the flux divergence must be equal to the sum of thickness change and mass balance rate as
$$
\nabla \cdot q = \dpar{h}{t} - \dot b.
$$
Quantifying the ice flux from estimated distributions of thickness change and mass balance rate and using Glen's flow law as constitutive equation, it is possible to solve for the ice thickness $h$ as
$$
h = \sqrt[n+2]{\frac{q}{2A}\frac{n+2}{(f\rho g \sin\alpha)^n}},\label{eq:thickness_mass_conservation}
$$
with ice density $\rho$, gravitational acceleration $g$, surface slope $\alpha$, shape factor $f$ and constant $A$ and exponent $n$ from Glen's flow law. In general, the ice thickness is estimated along flowlines and then interpolated over the entire glacier. Models in this category: [Farinotti et al. (2009)][], [Maussion et al. (2019)][], [Huss and Farinotti (2012)][], [Clarke et al. (2013)][], [Morlighem et al., 2011][]

**Shear-stress-based approaches:** The shallow ice approximation related the basal shear stress $\tau$ to ice thickness $h$ and glacier slope $\alpha$ as $\tau = \rho g h \sin\alpha$. The shear stress is estimated using an empirical relation by [Haeberli and Hoelzle (1995)][], relating the driving stress $\tau$ to the elevation range of the glacier. This allows to solve for the ice thickness $h$, the main difference to prior category is the missing account for mass conservation. Models in this category: [Linsbauer et al., 2009][], [Frey et al., 2014][] and an independent reimplementation of it.

**Velocity-based approaches:** Assuming parallel flow and Glen's flow law as constitutive equation, the vertical velocity profile can be derived assuming a linear increase in shear stress with depth. Integrating from bed up to the surface yields a relation for the surface velocity as
$$
u_s = u_b + \frac{2A}{n+1}\tau^nh.
$$
An alternative form is obtained by replacing the ice flux $q = \bar u h$ in Eq. ($\ref{eq:thickness_mass_conservation}$) with the product of depth average velocity $\bar u$ and ice thickness $h$. Given an assumption that relates the surface velocity $u_s$ to the basal velocity $u_b$ or  the average velocity $\bar u$ allows to solve for the ice thickness. Models in this category: three different implementations based on [Gantayat et al. (2014)][] and [Rabatel et al. (2018)][]

**Other approaches:** There are two additional models/approaches that fall outside all above mentioned categories. [Clarke et al. (2009)][] use an artificial neural network, disregarding any kind of ice physics. [Brinkerhoff et al. (2016)][] approach the ice thickness inversion from the framework of Bayesian interference, modeling bed elevations and ice flux divergence as Gaussian random fields with assumed covariance and unknown mean.

#### Results and Discussion:

**Between-model comparison**

Even models of the same category showed a strong spread amongst them, which hints that the models are independent even if based on similar/the same conceptual principles.

**Comparison to ice thickness measurements:** 

The ensemble average estimation matches direct measurements quite well, with an average deviation of less than 10% for 17 out of 21 cases. This can be explained by the law of large numbers, i.e., the random unbiased errors of different models cancel each other out.

For the following test cases the ensemble averaged ice thickness is significantly smaller than observed: Unterarr, Tasman, Washmawapta and Urumqi. Possible explanations include debris cover, the highly branched nature of the glacier, accumulation from surrounding steep and ice free slopes (avalanches, snow drift, ...) and the cold nature of the glacier (wrong flow rate factor).

The ice thickness estimation along the flowline is generally closer to the observations, then the distribution across. This can explained by the use of surface slope as predictor, which controls mostly the along-flow thickness.

In general, the synthetic cases yield better result. Given that the synthetic glaciers "grew" from similar models that are used for the inversion, and that the input data is not subjected to any observational uncertainties this was to be expected. Two additional points are worth pointing out: (1) the no slip assumption is well justified for the synthetic cases, less so for the real glaciers and (2) the synthetic glaciers are all in steady state, while the real glaciers hardly are.

As a comparison, the well know and highly used volume/area scaling model with $\bar h = c A^{\gamma-1}$ (using standard values of $c=0.034\ (0.054)$ and $\gamma = 1.36\ (1.25)$ for glaciers (ice caps), see e.g., [Bahr et al. (2015)][]) deviates from the measured ice thickness by –42±59% which is more than double the spread produced by the ensemble model average of 10±24%. The highly negative mean could be changed by using a different scaling constant $c$, however the spread will stay the same

**Individual model performance** ...



#### Introduction:

...

#### Experimental setup:

...

#### Considered test cases and data:

...

#### References:

[Brinkerhoff et al. (2016)]: https://doi.org/10.3389/feart.2016.00008	"Brinkerhoff DJ, Aschwanden A and Truffer M (2016) Bayesian Inference of Subglacial Topography Using Mass Conservation. Front. Earth Sci. 4:8. doi: 10.3389/feart.2016.00008"
[Van Pelt et al. (2013)]: https://doi.org/10.5194/tc-7-987-2013	"van Pelt, W. J. J., Oerlemans, J., Reijmer, C. H., Pettersson, R., Pohjola, V. A., Isaksson, E., and Divine, D.: An iterative inverse method to estimate basal topography and initialize ice flow models, The Cryosphere, 7, 987–1006, https://doi.org/10.5194/tc-7-987-2013, 2013."
[Fürst et al. (2017)]: https://doi.org/10.5194/tc-11-2003-2017	"Fürst, J. J., Gillet-Chaulet, F., Benham, T. J., Dowdeswell, J. A., Grabiec, M., Navarro, F., Pettersson, R., Moholdt, G., Nuth, C., Sass, B., Aas, K., Fettweis, X., Lang, C., Seehaus, T., and Braun, M.: Application of a two-step approach for mapping ice thickness to various glacier types on Svalbard, The Cryosphere, 11, 2003–2032, https://doi.org/10.5194/tc-11-2003-2017, 2017."
[Farinotti et al. (2009)]: https://doi.org/10.3189/002214309788816759	"Farinotti, D., Huss, M., Bauder, A., Funk, M., &amp; Truffer, M. (2009). A method to estimate the ice volume and ice-thickness distribution of alpine glaciers. Journal of Glaciology, 55(191), 422-430. doi:10.3189/002214309788816759"
[Maussion et al. (2019)]: https://doi.org/10.5194/gmd-12-909-2019	"Maussion, F., Butenko, A., Champollion, N., Dusch, M., Eis, J., Fourteau, K., Gregor, P., Jarosch, A. H., Landmann, J., Oesterle, F., Recinos, B., Rothenpieler, T., Vlug, A., Wild, C. T., and Marzeion, B.: The Open Global Glacier Model (OGGM) v1.1, Geosci. Model Dev., 12, 909–931, https://doi.org/10.5194/gmd-12-909-2019, 2019."
[Huss and Farinotti (2012)]: https://doi.org/10.1029/2012JF002523	"Huss, M., and Farinotti, D. (2012), Distributed ice thickness and volume of all glaciers around the globe, J. Geophys. Res., 117, F04010, doi:10.1029/2012JF002523."
[Clarke et al. (2013)]: 
[Morlighem et al., 2011]: https://doi.org/10.1029/2011GL048659	"Morlighem, M., Rignot, E., Seroussi, H., Larour, E., Dhia, H. B., and Aubry, D.: A mass conservation approach for map- ping glacier ice thickness, Geophys. Res. Lett., 38, L19503, doi:10.1029/2011GL048659, 2011."
[Haeberli and Hoelzl (1995)]: 
[Linsbauer et al. (2009)]: 
[Frey et al. (2014)]: 
[Gantayat et al. (2014)]: 
[Rabatel et al. (2018)]: 
[Clarke et al. (2009)]: https://doi.org/10.1175/2008JCLI2572.1	"Clarke, G. K. C., Berthier, E., Schoof, C. G., and Jarosch, A. H.: Neural networks applied to estimating subglacial topography and glacier volume, J. Climate, 22, 2146–2160, doi:10.1175/2008JCLI2572.1, 2009."