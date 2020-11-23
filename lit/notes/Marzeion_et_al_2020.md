# Partitioning the uncertainty of ensemble projections of global glacier mass change

**Ben Marzeion et al. (2020) - Earth's future - https://doi.org/10.1029/2019EF001470**

## Plain language summary:

> Melting glaciers (outside the large ice sheets in Greenland and Antarctica) contribute strongly to rising sea level and are expected to continue to do so throughout this century. However, the amount of future sea level rise from glaciers is not well known. One of the causes for uncertainty is the lack of knowledge of future greenhouse gas emissions. This uncertainty is growing steadily during the 21st century and constitutes the most important uncertainty by 2100. Another cause of uncertaint[ies] are the glacier models themselves, since they rely on approximations and simplifications of complex glacier processes. This uncertainty is very important until the middle of the 21st century, but less important in the second half of the 21st century. Overall, glaciers will lose around 18 % of their ice mass in a low-emission scenario, or around 36 % in a high-emission scenario, contributing roughly 79 or 159 mm to sea level rise by 2100.
>
> Marzeion et al. (2020)

## Abstract:

> Glacier mass loss is recognized as a major contributor to current sea level rise. However, large uncertainties remain in projections of glacier mass loss on global and regional scales. We present an ensemble of 288 glacier mass and area change projections for the 21st century based on 11 glacier models using up to 10 general circulation models and four Representative Concentration Pathways (RCPs) as boundary conditions. We partition the total uncertainty into the individual contributions caused by glacier models, general circulation models, RCPs, and natural variability. We find that emission scenario uncertainty is growing throughout the 21st century and is the largest source of uncertainty by 2100. The relative importance of glacier model uncertainty decreases over time, but it is the greatest source of uncertainty until the middle of this century. The projection uncertainty associated with natural variability is small on the global scale but can be large on regional scales. The projected global mass loss by 2100 relative to 2015 (79 ± 56 mm sea level equivalent for RCP2.6, 159 ± 86 mm sea level equivalent for RCP8.5) is lower than, but well within, the uncertainty range of previous projections.
>
> Marzeion et al. (2020)

An ensemble of modeled glacier mass change is used to asses the uncertainty in future projections. The ensemble consists of 11 glacier models, up to 10 GCMs (general circulation models) and 4 RPCs (representative concentration pathways). The following key points are found:

- The contribution of glacier melt to sea level rise from 2015 until 2100 ranges from 79 ± 56 mm under RCP2.6 to 159 ± 86 mm under RCP8.6
- The difference between glacier models accounts for the main part of uncertainty in the first half of the century
- The uncertainty in global emission scenarios is the major contributor to the total uncertainty for the second half of the century

## 6. Conclusions

A major part of the total uncertainty in projected glacier mass loss for the second half of the 21st century can be attributed to the emission scenarios. Uncertainties in GCMs and glacier model play a secondary role, while natural variability may be neglected.

For the first half of the 21st century, the total uncertainty is dominated by differences in glacier models. These differences diminish over time. Hence, improving glacier models should have low priority for long term projections but a higher priority for short term projections. (Different wording: This indicates that long term projections won't significantly benefit from improved glacier models, while short term projections will.) The main problem for short term projections seems to be the imbalance between glacier geometry and the concurrent climate (cf. [Mernild et al. (2013)][], [Marzeion et al. (2018)][]).

The projected glacier mass loss is below past studies but well inside the uncertainty range. For detailed values see key points in abstract.

## 2. Data and methods

### 2.2 Glacier models

A total of 11 glacier models provide projections of glacier area and ice mass until 2100 for all primary RGI regions. The following list shows the model characteristics. If not further specified, the models resolve each individual glacier, do not account for frontal ablation and use monthly time steps for mass balance and annual time steps for geometry changes.

#### Global models (with or without Antarctica)

The following models cover either all RGI regions or all regions except Antartica and Sub-Antarctica (MAR2012, JULES, OGGM).

- [x] **van de Wal and Wild (2001) - WAL2001**: volume/area scaling, mass balance parametrization based on summer and winter temperatures, uniform precipitation, annual time steps for geometry and mass balance
- [x] **Marzeion et al. (2012) - MAR2012**: volume/area, volume/length and response time scaling, temperature-index model, linear precipitation gradient and threshold temperatures, accounts for elevation range of each glacier
- [x] **Radic et al. (2014) - RAD2014**: volume/area and volume/length scaling, temperature-index model, linear precipitation gradient and threshold temperatures, resolves 20 m elevation bands for each glacier
- [x] **Huss and Hock (2015) - GloGEM**: mass-conserving parameterized surface elevation change, temperature-index model, linear precipitation gradient and threshold temperatures, resolves 10 m elevation bands for each glacier
- [x] **Sakai and Fujita (2017) - GLIMB**: volume-area scaling, energy balance model, uniform precipitation and liner probability function of temperature, resolving individual geometry for each glacier but mass balance only on a 0.5° grid with 50 m elevation bands, daily time steps for the mass balance model
- [x] **Shannon et al. (2019) - JULES**: thickness change equates with surface mass balance (in ice equivalent), energy balance model, linear precipitation gradient and threshold temperatures, runs on a 0.5° grid with 250 m elevation bands and hourly time steps
- [x] **Maussion et al. (2019) - OGGM**: flowline model (SIA), temperature-index model, uniform precipitation and threshold temperatures, resolves each glacier with variable elevation bands (distance of 20 - 400 m on flowline), monthly time steps for mass balance and variable time steps for geometry model (to ensure numerical stability) but annual geometry update

The models WAL2001 (or SLA2012), MAR2012, RAD2014 and GloGEM already took part in the first glacier model intercomparison (GlacierMIP, Hock et al., 2019). Only the OGGM and the GloGEMflow update the glacier geometries using a flowline model. The $\Delta h$-method used by the GloGEM model is based on an empirically based mass/geometry change relation. However, it is the only method to account for frontal ablation of water-terminating glaciers. All the other models rely on different forms of volume/area scaling to update glacier geometries.

Most mass balance models are a combination of temperature-index models for melt and temperature thresholds to distinguish between solid and liquid precipitation. Different bias corrections are applied to the given GCM boundary conditions, model parameters are calibrated to maximize agreement with observations, independent uncertainty estimations are performed.

The mass loss is converted directly into sea-level equivalent, not accounting for any interception of terrestrial runoff (e.g., ground water replenishment, ...) or ice grounded below sea level (only GloGEM and OGGM are capable of it anyway).

#### Regional models

...

- **Anderson and Mackintosh (2012) - AND2012**: local model for New Zealand, flowline model using SIA, temperature-index model, uniform precipitation for each glacier and threshold temperatures, accounts for frontal ablation using the height above floatation criterion, 100 m grid, daily time steps for mass balance
- **Kraaijenbrink et al. (2017) - KRA2017**: local model for High Mountain Asia, volume/area relation per elevation band based on bed topography, mass balance gradient based on temperature index, uniform precipitation, resolves each glacier with individual elevation bands, annual time steps
- **Zekollari et al. (2019) - GloGEMflow**: local model for Central Europe, temperature-index model, linear precipitation gradient and threshold temperatures, resolves each glacier with variable elevation bands (distance of 10 - 202 m on flowline),
- **Rounce et al. (2020) - PyGEM**: local model for High Mountain Asia, mass-conserving parametrized surface elevation change, temperature-index model, linear precipitation gradient with threshold temperatures, resolves 10 m elevation bands for each glacier

## 3. Glacier mass and area projections

All glacier models are provided with the same initial glacier volume. Due to the different starting times and individual ice evolution of the different model, the global ice mass estimate at the start of the study period in 2015 varies between 391 mm SLE (GLIMB) and 551 mm SLE (RAD2014). Additionally, the mean of all global models 439 mm SLE is higher than the latest consensus estimate of 382 ± 99 mm SLE (Farinotti et al., 2019) but within the uncertainty range.

As was to be expected, the specific mass balance depends strongly on the emission scenario. For all regions, the specific mass balance is constantly negative or slightly increasing towards zero under RCP2.6 and strongly decreasing under RCP8.5. Projected glacier area decreases accordingly over the 21st century for all regions and emission scenarios, with stronger relative area loss for regions with smaller glaciers (e.g., Alps). This may result in a new equilibrium under RCP2.6 or an almost complete glacier loss under RCP8.5.



## References

[Mernild et al. (2013)]: https://doi.org/10.5194/tc-7-1565-2013	"Mernild, S. H., Lipscomb, W. H., Bahr, D. B., Radi´c, V., &amp; Zemp, M. (2013). Global glacier changes: A revised assessment of committed mass losses and sampling uncertainties. The Cryosphere, 7, 1565–1577."
[Marzeion et al. (2018)]: https://doi.org/10.1038/s41558-018-0093-1	"Marzeion, B., Kaser, G., Maussion, F., &amp; Champollion, N. (2018). Limited influence of climate change mitigation on short-term glacier mass loss. Nature Climate Change, 8, 305–308."

