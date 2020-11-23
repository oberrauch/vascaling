# Glacier volume-area scaling papers

- **Adhikari and Marsahall (2012):** *Glacier volume-area relation for high order mechanism and transient glacier states*

  - [ ] Read

  **Strength**: experiments compared to a full 3-D model; restriction to ensembles of glaciers (not individuals) recognized. Correctly identifies glacier shape and slope as drivers of volume (closure conditions in [Bahr et al. (1997)][] and [Bahr et al (2015)][])

  **Errors**: the scaling exponent $\gamma$ is treated as variable and some derived values fall outside theoretical meaningful limits (e.g., $\gamma > 2$); additionally claims that volume-area scaling theory is limited to SIA and valid only for steady state conditions.

  **Abstract**: Glacier volume is known for less than 0.1% of the world's glaciers, but this information is needed to quantify the impacts of glacier changes on global sea level and regional water resources. Observations indicate a power-law relation between glacier area and volume, with an exponent 1.36. Through numerical simulations of 3D, high-order glacier mechanics, we demonstrate how different topographic and climatic settings, glacier flow dynamics, and the degree of disequilibrium with climate systematically affect the volume-area relation. We recommend more accurate scaling relations through characterization of individual glacier shape, slope and size. An ensemble of 280 randomly generated valley glaciers spanning a spectrum of plausible glaciological conditions yields a steady-state exponent = 1.46. This declines to 1.38 for glaciers that are 100 years into a sustained retreat, which corresponds exceptionally well with the observed value for present-day glaciers.

  **Comment**: It seems to be, that the sole purpose of this paper is to calibrate the volume-area scaling exponent. Since the exponent is a theoretical constant $\gamma$, the study is of little interest... A possible benefit could be to use high and low values $\gamma$ to identify glacier characteristics for which the volume-area scaling is less accurate.

- **Farinotti et al. (2009):** *A method to estimate the ice volume and ice thickness distribution of alpine glaciers*

  ==New method to determine alpine glacier ice thickness/volume based on mass turnover and principles of ice-flow mechanics.== 

  - [x] Read (yes, but years ago while working on my Bachelors thesis)

  **Strength**: the developed method is applicable to single glaciers, and compared to the scaling method.

  **Errors**: volume-area scaling theory wrongly assumed to be limited to steady-state conditions and applied to single glaciers. (8.1, 8.5)

  **Abstract**:  Sound knowledge of the ice volume and ice-thickness distribution of a glacier is essential for many glaciological applications. However, direct measurements of ice thickness are laborious, not feasible everywhere and necessarily restricted to a small number of glaciers. In this paper, we present a method to estimate the ice-thickness distribution and the total ice volume of alpine glaciers. This method is based on glacier mass turnover and principles of ice-flow mechanics. The required input data are the glacier surface topography, the glacier outline and a set of borders delineating different "ice-flow catchments". Three parameters describe the distribution of the "apparent mass balance", which is defined as the difference between the glacier surface mass balance and the rate of ice-thickness change, and two parameters define the ice-flow dynamics. The method was developed and validated on four alpine glaciers located in Switzerland, for which the bedrock topography is partially known from radio-echosoundings. The ice thickness along 82 cross-profiles can be reproduced with an average deviation of about 25% between the calculated and the measured ice thickness. The cross-sectional areas differ by less than 20% on average. This shows the potential of the method for estimating the ice-thickness distribution of alpine glaciers without the use of direct measurements.

  **Comments**: The OGGM ice thickness inversion is based on this study. The estimated ice volume for Rhone-, Silvretta- and Unteraar Glacier is larger (by 37% on average) than the volume determined by the volume-area scaling method. The estimated ice volume of the highly branched Glacier de Zinal is lower than the scaling counterpart. However, volume-area scaling should not be applied to single glaciers and a "case study" with $n=4$ has little to no meaning.

- **Farinotti and Huss (2013):** *An upper bound estimate for the accuracy of glacier volume-area scaling*

  ==Estimate of the statistical accuracy of total volume estimations derived from volume-area scaling==

  - [x] Read

  **Strength**: Performs a rigorous statistical analysis of the accuracy of volume-area scaling. Correctly concludes that scaling should be applied to large populations of glaciers and that allowing the scaling exponent $\gamma$ to vary with time is less accurate than treating $\gamma$ as a constant.

  **Errors**: Assumes scaling exponent $\gamma$ is a spatial and temporal variable or a constant other than the theoretical value of 1.375. Estimates $\gamma$ from data. The variability of $\gamma$ is built into the estimate of volume-area scaling accuracy. A balance gradient is specified for the 3D Stokes model, but its value may be inconsistent with the volume-area scaling closure condition. (8.2, 8.3, 8.7, 8.9, 8.11)

  **Abstract**: Volume–area scaling is the most popular method for estimating the ice volume of large glacier samples. Here, a series of resampling experiments based on different sets of synthetic data is presented in order to derive an upper-bound estimate (i.e., a level achieved only within ideal conditions) for its accuracy. For real-world applications, a lower accuracy has to be expected. We also quantify the maximum accuracy expected when scaling is used for determining the glacier volume change, and area change of a given glacier population. A comprehensive set of measured glacier areas, volumes, area and volume changes is evaluated to investigate the impact of real-world data quality on the so-assessed accuracies. For populations larger than a few thousand glaciers, the total ice volume can be recovered within 30% if all data currently available worldwide are used for estimating the scaling parameters. Assuming no systematic bias in ice volume measurements, their uncertainty is of secondary importance. Knowing the individual areas of a glacier sample for two points in time allows recovering the corresponding ice volume change within 40% for populations larger than a few hundred glaciers, both for steady-state and transient geometries. If ice volume changes can be estimated without bias, glacier area changes derived from volume–area scaling show similar uncertainties to those of the volume changes. This paper does not aim at making a final judgement on the suitability of volume–area scaling as such, but provides the means for assessing the accuracy expected from its application.

  **Comments**: 

- **Grinsted (2013)**: *Global glacier volume estimate*

  - [ ] Read

  **Strength**: Good compendium of previously published analyses; theoretical insights on role of surface slope; acknowledges that treating glacier complexes will degrade accuracy

  **Errors**: Attempts to assign multiple values to fixed scaling exponent $\gamma$ by empirical means and introduces additional parameters with no theoretical basis; treats glacier complexes as single entities.
  
  **Abstract**: I assess the feasibility of using multivariate scaling relationships to estimate glacier volume from glacier inventory data. Scaling laws are calibrated against volume observations optimized for the specific purpose of estimating total global glacier ice volume. I find that adjustments for continentality and elevation range improve skill of area–volume scaling. These scaling relationships are applied to each record in the Randolph Glacier Inventory, which is the first globally complete inventory of glaciers and ice caps. I estimate that the total volume of all glaciers in the world is 0.35 ± 0.07 m sea level equivalent, including ice sheet peripheral glaciers. This is substantially less than a recent state-of-the-art estimate. Area–volume scaling bias issues for large ice masses, and incomplete inventory data are offered as explanations for the difference.
  
  **Comments**:
  
- **Harrison (2013)**: *How do glaciers respond to climate? Perspectives from the simplest models*

  - [ ] Read

  **Strength**: A macroscopic scaling-type derivation with a non-dimensional parameter consistent with [Bahr et al. (2015)][]. Correctly notes that volume-area scaling exponent is a constant, and combines volume scaling with response time scaling.

  **No errors or weaknesses**

  **Abstract**: We study the approximations underlying the macroscopic theory of glacier response to climate, and illustrate some basic properties of glacier response using two simple examples, a block on an inclined plane and a block with a more realistic terminal region. The properties include nonlinearity and the usefulness of linear approximations, sensitivity to bed slope, timescale, stability, characteristic elevations of the equilibrium line at which the properties of the response change, the limit of fast response in a changing climate, the minimum sustainable size, increased sensitivity to climate change as a glacier retreats, and finally the similarity in the responses of simple glaciers with the same product of length and square of bed slope.

  **Comments**: Difficult to read and understand...

- **Huss and Farinotti (2012)**: *Distributed ice thickness and volume of all glaciers around the globe*

  ==First physically based global ice volume estimate 170 × 10^3^ ± 21 × 10^3^ km^3^ (0.43 ± 0.06 m SLE)==

  - [ ] Read

  **Strength**: Method is applicable to single glaciers. Volume-area scaling correctly noted to be inapplicable to multiple-glacier complexes treated as single entity.

  **Errors**: Volume-area scaling wrongly assumed to be limited to steady state conditions and to be unsuitable for complex geometries such as dendritic glaciers. Volume-area scaling exponent assumed to be spatially variable and treated as tunable parameter

  **Abstract**: A new physically based approach for calculating glacier ice thickness distribution and volume is presented and applied to all glaciers and ice caps worldwide. Combining glacier outlines of the globally complete Randolph Glacier Inventory with terrain elevation models (Shuttle Radar Topography Mission/Advanced Spaceborne Thermal Emission and Reflection Radiometer), we use a simple dynamic model to obtain spatially distributed thickness of individual glaciers by inverting their surface topography. Results are validated against a comprehensive set of thickness observations for 300 glaciers from most glacierized regions of the world. For all mountain glaciers and ice caps outside of the Antarctic and Greenland ice sheets we find a total ice volume of 170 × 10^3^ ± 21 × 10^3^ km^3^, or 0.43 ± 0.06 m of potential sea level rise.

  **Comments**: Compares findings to ice volume estimates from scaling.

- **Lüthi (2009)**: *Transient response of idealized glaciers to climate change*

  - [ ] Read

  **Strength**: An alternative theoretical derivation of volume-area, volume-length and response time scaling for a limited geometry. Results are largely consistent with [Bahr et al. (2015)][] and with a 3D Stokes model. Contains a theoretical derication of the scaling constant $c$ for the specified geometry. Correctly applies response time and volume-area scaling simultaneously.

  **Errors**: Theory is limited to an idealized geometry.

  **Abstract**:  The transient response of glaciers to climate variations is investigated with a novel low-order model of an idealized glacier resting on uniformly inclined bedrock. The model consists of two volumes representing the accumulation and ablation areas, which are joined by a flux gate controlling mass flow according to the shallow-ice approximation. Under the assumption of a constant vertical mass balance gradient, a volume-length scaling relation is derived which depends explicitly on mass balance gradient and bedrock slope. Analytic expressions for the volume and area timescales are given, which are inversely proportional to the mass-balance gradient and a geometric factor which is the ratio between the vertical extent of the ablation area and the ice thickness at the equilibrium line. From the low-order model, a dynamical system in length and volume is obtained. Results from this system are in good agreement with solutions obtained from a transient finite-element model solving the full forcebalance and mass-conservation equations. Under periodic forcing there are significant deviations of the length response, which show that the usual relaxation-type parameterization of length change is not well suited for short-term reactions. Switching on and off periodic climate forcing, the model glaciers show surprisingly large initial and final transient responses that have not been investigated before. These results are of significance for the interpretation of length and thickness changes observed on glaciers.

  **Comments**:

- **Marzeion et al. (2012)**: *Past and future sea-level change from the surface mass balance of glaciers*

  - [ ] Read

  **Strength**: Robust treatment of time scale of transients, proper scaling parameters adopted.

  **No errors or weaknesses**

  **Comment**: The main paper which sparked the idea for my thesis...

- **Radic and Hock (2010)**: *Regional and global volumes of glaciers derived from statisticall upscaling of glacier inventory data*

  - [ ] Read

  **Strength**: Selects appropriate values for the scaling exponent $\gamma$ and for ice caps.

  **Errors**: Assumes scaling exponent $\gamma$ has associated uncertainty.

  **Abstract**: Very few global‐scale ice volume estimates are available for mountain glaciers and ice caps, although such estimates are crucial for any attempts to project their contribution to sea level rise in the future. We present a statistical method for deriving regional and global ice volumes from regional glacier area distributions and volume area scaling using glacier area data from ∼123,000 glaciers from a recently extended World Glacier Inventory. We compute glacier volumes and their sea level equivalent (SLE) for 19 glacierized regions containing all mountain glaciers and ice caps on Earth. On the basis of total glacierized area of 741 × 10^3^ ± 68 × 10^3^ km^2^, we estimate a total ice volume of 241 × 10^3^ ± 29 × 10^3^ km3, corresponding to 0.60 ± 0.07 m SLE, of which 32% is due to glaciers in Greenland and Antarctica apart from the ice sheets. However, our estimate is sensitive to assumptions on volume area scaling coefficients and glacier area distributions in the regions that are poorly inventoried, i.e., Antarctica, North America, Greenland, and Patagonia. This emphasizes the need for more volume observations, especially of large glaciers and a more complete World Glacier Inventory in order to reduce uncertainties and to arrive at firmer volume estimates for all mountain glaciers and ice caps.

  **Comments**:

- **Radic and Hock (2011)**: *Regionally differentiated contribution of mountain glaciers and ice caps to future sea-level rise*

  - [ ] Read

  **Strength**: Selects appropriate values for the scaling exponent $\gamma$ and for ice caps; notes the potential time variation in the scaling constant $c$.

  **Errors**: Assumes scaling exponent $\gamma$ has associated uncertainty.

  **Abstract**: The contribution to sea-level rise from mountain glaciers and ice caps has grown over the past decades. They are expected to remain an important component of eustatic sea-level rise for at least another century, despite indications of accelerated wastage of the ice sheets. However, it is difficult to project the future contribution of these small-scale glaciers to sea-level rise on a global scale. Here, we project their volume changes due to melt in response to transient, spatially differentiated twenty-first century projections of temperature and precipitation from ten global climate models. We conduct the simulations directly on the more than 120,000 glaciers now available in the World Glacier Inventory, and upscale the changes to 19 regions that contain all mountain glaciers and ice caps in the world (excluding the Greenland and Antarctic ice sheets). According to our multi-model mean, sea-level rise from glacier wastage by 2100 will amount to 0.124±0.037 m, with the largest contribution from glaciers in Arctic Canada, Alaska and Antarctica. Total glacier volume will be reduced by 21±6%, but some regions are projected to lose up to 75% of their present ice volume. Ice losses on such a scale may have substantial impacts on regional hydrology and water availability.

- **Radic et al. (2007)**: *Volume–area scaling vs flowline modelling in glacier volume projections*

  - [ ] Read

  **Strength**: Acknowledges potential for errors from use of 1D model.

  **Errors**: A 1D glacier model is used (volume-area scaling parameters from this mode will not necessarily apply to real glaciers); unrealistic wide range of values for scaling exponent $\gamma$ found. Assumes scaling exponent $\gamma$ is time dependent for non-equilibrium conditions. Assumes volume-area scaling applies only to steady state conditions. (8.1, 8.2, 8.3, 8.9)

  **Abstract**: Volume–area scaling provides a practical alternative to ice-flow modelling to account for glacier size changes when modelling the future evolution of glaciers; however, uncertainties remain as to the validity of this approach under non-steady conditions. We address these uncertainties by deriving scaling exponents in the volume–area relationship from one-dimensional ice-flow modelling. We generate a set of 37 synthetic steady-state glaciers of different sizes, and then model their volume evolution due to climate warming and cooling as prescribed by negative and positive mass-balance perturbations, respectively, on a century timescale. The scaling exponent derived for the steady-state glaciers ( $\gamma$ = 1.56) differs from the exponents derived for the glaciers in transient (non-steady) state by up to 86%. Nevertheless, volume projections employing volume–area scaling are relatively insensitive to these differences in scaling exponents. Volume–area scaling agrees well with the results from ice-flow modelling. In addition, the scaling method is able to simulate the approach of a glacier to a new steady state, if mass-balance elevation feedback is approximated by removing or adding elevation bands at the lowest part of the glacier as the glacier retreats or advances. If area changes are approximated in the mass-balance computations in this way, our results indicate that volume–area scaling is a powerful tool for glacier volume projections on multi-century timescales.

  **Comments**:

- **Radic et al. (2008)**: *Analysis of scaling methods in deriving future volume evolutions of valley glaciers*

  - [ ] Read

  **Strength**: Volume-area scaling adapted to non-steady conditions; application of multiple scaling relations to single glacier.

  **Errors**: Scaling exponent $\gamma$ treated as adjustable parameter and assumes that $\gamma$ has to vary to account for non-equilibirum conditions.

  **Abstract**: Volume–area scaling is a common tool for deriving future volume evolutions of valley glaciers and their contribution to sea-level rise. We analyze the performance of scaling relationships for deriving volume projections in comparison to projections from a one-dimensional ice-flow model. The model is calibrated for six glaciers (Nigardsbreen, Rhonegletscher, South Cascade Glacier, Sofiyskiy glacier, midre Love´nbreen and Abramov glacier). Volume evolutions forced by different hypothetical mass-balance perturbations are compared with those obtained from volume–area (V-A), volume–length (V-L) and volume–area–length (V-A-L) scaling. Results show that the scaling methods mostly underestimate the volume losses predicted by the ice-flow model, up to 47% for V-A scaling and up to 18% for V-L scaling by the end of the 100 year simulation period. In general, V-L scaling produces closer simulations of volume evolutions derived from the ice-flow model, suggesting that V-L scaling may be a better approach for deriving volume projections than V-A scaling. Sensitivity experiments show that the initial volumes and volume evolutions are highly sensitive to the choice of the scaling constants, yielding both over- and underestimates. However, when normalized by initial volume, volume evolutions are relatively insensitive to the choice of scaling constants, especially in the V-L scaling. The 100 year volume projections differ within 10% of initial volume when the V-A scaling exponent commonly assumed, ? ¼ 1.375, is varied by –30% to +45% (? ¼ [0.95, 2.00]) and the V-L scaling exponent, q ¼ 2.2, is varied by –30% to +45% (q ¼ [1.52, 3.20]). This is encouraging for the use of scaling methods in glacier volume projections, particularly since scaling exponents may vary between glaciers and the scaling constants are generally unknown.

  **Comments**:

- **Slangen and van de Wal (2011)**: *An assessment of uncertainties in using volume-area modelling for computing the twenty-first century glacier contribution to sea-level change*

  - [ ] Read

  **Strength**: Fixes scaling exponent $\gamma$ at appropriate values and allows the scaling constant $c$ to vary.

  **Errors**: Incorrectly performs a sensitivity analysis on the scaling exponent $\gamma$ which should remain constant.

  **Abstract**: A large part of present-day sea-level change is formed by the melt of glaciers and ice caps (GIC). This study focuses on the uncertainties in the calculation of the GIC contribution on a century timescale. The model used is based on volume-area scaling, combined with the mass balance sensitivity of the GIC. We assess different aspects that contribute to the uncertainty in the prediction of the contribution of GIC to future sea-level rise, such as (1) the volume-area scaling method (scaling factor), (2) the glacier data, (3) the climate models, and (4) the emission scenario. Additionally, a comparison of the model results to the 20th century GIC contribution is presented.

  We find that small variations in the scaling factor cause significant variations in the initial volume of the glaciers, but only limited variations in the glacier volume change. If two existing glacier inventories are tuned such that the initial volume is the same, the GIC sea-level contribution over 100 yr differs by 0.027 m or 18 %. It appears that the mass balance sensitivity is also important: variations of 20 % in the mass balance sensitivity have an impact of 17 % on the resulting sea-level projections. Another important factor is the choice of the climate model, as the GIC contribution to sea-level change largely depends on the temperature and precipitation taken from climate models. Connected to this is the choice of emission scenario, used to drive the climate models. Combining all the uncertainties examined in this study leads to a total uncertainty of 0.052 m or 35 % in the GIC contribution to global mean sea level. Reducing the variance in the climate models and improving the glacier inventories will significantly reduce the uncertainty in calculating the GIC contributions, and are therefore crucial actions to improve future sea-level projections.

- **van de Wal and Wild (2001)**: *Modelling the response of glaciers to climate change by applying volume-area scaling in combination with a high resolution GCM*

  - [ ] Read

  **Strength**: One of the first papers to pioneer the use of volume-area scaling to estimate sea level (change). Fixes exponent $\gamma$ at appropriate values for glaciers and ice caps. Recognizes that volume-area scaling applies in non-equilibrium conditions.

  **Errors**: Incorrectly performs a sensitivity analysis on the scaling exponent $\gamma$ which should remain constant.

  **Abstract**: A seasonally and regionally differentiated glacier model is used to estimate the contribution that glaciers are likely to make to global sea level rise over a period of 70 years. A high resolution general circulation model (ECHAM4 T106) is used to estimate temperature and precipitation changes for a doubled CO~2~ climate and serves as input for the glacier model. Volume-area relations are used to take into account the reduction of glacier area resulting from greenhouse warming. Each glaciated region has a specified glacier size distribution, defined by the number of glaciers in a size class and a mean area. Changes in glacier volume are calculated by a precipitation dependent mass balance sensitivity. The model predicts a global sea level rise of 57 mm over a period of 70 years. This corresponds to a sensitivity of 0.86 mm yr^-1^K^-1^. Assuming a constant glacier area as done in earlier work leads to an overestimation of 19% for the contribution to sea level rise.

  

#### References:

[Bahr et al. (1997)]: http://doi.org/10.1029/97JB01696	"Bahr, D. B., M. F. Meier, and S. D. Peckham (1997), The physical basis of glacier volume-area scaling, J. Geophys. Res., 102(B9), 20,355–20,362"
[Bahr et al. (2015)]: http://doi.org/10.1002/2014RG000470	"A review of volume-area scaling of glaciers"