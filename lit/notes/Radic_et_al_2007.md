## Volume–area scaling vs. flowline modeling in glacier volume projections

**Valentina Radić, Regine Hock and Johannes Oerlemans (2007) - Annals of Glaciology - http://doi.org/10.3189/172756407782871288**

#### **Abstract:**

The scaling exponent for transient glaciers is derived from a set of 37 synthetic glacier modeled by one-dimensional flowline model. The glacier evolution is forced by positive and negative mass balance perturbations on a century timescale. The resulting *transient* scaling exponent differs from the *steady state* exponent by up to 86%, however the ice volume projections from volume/area scaling are rather insensitive to said differences. Approximating mass-balance elevation feedback by adding or removing elevation band enables the simulation of a new steady state.

**Note:** Volume/area scaling does not assume anything, especially not a steady state condition. Volume/area scaling must be used in conjuncture with response time scaling if used for volume projections. Additionally, the scaling exponent $\gamma$ is a physically based constant (cf. Buckingham Pi theorem) and cannot be changed without changing the underlying closure conditions. If changes in the scaling exponent can be justified, there is still a physical sensible range $1 \leq \gamma \leq 1.5$ which cannot be exceeded. The scaling constant $c$ is a random variable and changes drastically from glacier to glacier. Therefore, a set of 37 glaciers is most likely to small for a sound determination of the scaling parameters. For more details see [Bahr et al. (2015)][] section 8.1, 8.2, 8.3 and 8.9.

**Private note:** On of the most nonsensical paper I've read... either I'm too stupid or ignorant to understand it (which is more than likely) or it really isn't worth the paper it is printed on (given that it's a PDF that means something). What bothers me the most, is the lack of explanations (or even just speculation) on why different scaling methods (calling them model is a bit of a stretch) behave the way they do... Could be interesting in performing similar experiments but with proper usage of the volume/area scaling regime (including response time scaling).

#### **Conclusions:**

> Scaling exponent $\gamma$ = 1.56 in the volume–area relationship obtained from 37 synthetic steady-state glaciers of different sizes differed from $\gamma$  = 1.375 derived theoretically by Bahr and others (1997) and from the exponents ($\gamma$ = [1.80, 2.90]) derived for each of 24 investigated glaciers under non-steady-state conditions, i.e. responding to hypothetical mass-balance perturbations.

> However, the range of differences in scaling exponent by up to 86% is shown to make negligible differences, <6%, in 100 year volume changes derived from the scaling approach.

> [...] simulating the glacier approaching a new steady state by simulating the feedback between area-averaged mass-balance and glacier geometry/elevation changes resulting from retreat or advance of the glacier. This feedback is captured by excluding area from or adding area to the lowest part of the glacier.

#### **Methods:**

**Flowline model:**

The one dimensional flowline model by [Oerlemans (1997)][] uses the vertically integrated continuity equation under the assumption of incompressibility and Glen's flow law as constitutive equation.
$$
\renewcommand{dpar}[2]{\frac{\partial#1}{\partial#2}}

\dpar{H}{t} = -\frac{1}{w}\dpar{}{x}\left[D\dpar{h}{x}\right] + b
$$
Hereby, $H$ is the ice thickness, $w$ the glacier width, $h$ the surface elevation, $b$ the mass balance rate and $D$ the diffusivity as
$$
D = w(\rho g)^3H^3\left(\dpar{h}{x}\right)^2\ \left(f_dH^2 + f_s\right),
$$
with $\rho$ as ice density and $g$ as gravitational acceleration. The vertical mean ice velocity is completely determined by internal deformation and sliding, as suggested by the deformation parameter $f_d$ and sliding parameter $f_s$.

**Set of synthetic glaciers:**

...

**Model-derived volume/area relationship:**

...

**Volume projections using volume/area scaling:**

The following "model" is used to estimate the ability of future volume projections using the volume/area scaling relation. The volume change $\Delta V(t) = \bar b(t) A(t)$ is given by the product of annual area-average net mass balance $\bar b(t)$ and surface area $A(t)$. The mass balance is calculated as sum of the mass balance of each elevation band $b_i$, weighted by its area $a_i$ and averaged over the total glacier area $A$. If total area is kept constant ($A(t) = A(t=0) = \text{const.}$) it is called the *reference surface mass balance*, and if the surface area is updated by inverting the scaling relation ($A(t) = [V(t-1) + \Delta V(t-1)]^{1/\gamma}$) it is called the *conventional mass balance*. Thereby, the area change is assumed to occur at the tongue area and implemented by adding or removing elevation bands of 100 m length (equal to the horizontal grid). The same distinction between constant surface area and updated surface area is made for the calculation of the volume change. The following three combinations are tested: (a) reference surface mass balance and no scaling (i.e constant area for volume change), (b) reference surface mass balance with scaling (i.e. updated area for volume change), and (c) conventional mass balance with scaling.

#### **Results and Discussion:**

**Volume/area relationship in steady state:**

The scaling exponent $\gamma_\text{synt}$ derived from the set of 37 synthetic glaciers equals $\gamma_\text{synt} = 1.56$ and thereby differs from the exponent derived by Bahr et al. (1997) of $\gamma_\text{Bahr}$ = 1.375 by 14%. Given the simplification in the one-dimensional flowline model and the uniform geometry, a certain difference from the empirical values was expected. **Note:** However, the value for the exponent falls outside the physical sensible range of $1\leq \gamma \leq 1.5$.

**Volume/area relationships in non-steady state:**

For each of the 24 glaciers subjected to mass balance perturbations, the yearly pairs of volume and area after a 50 year spin up period are used to compute volume/area scaling exponents for transient glaciers, yielding values $\gamma \in [1.80, 2.90]$. Thereby, larger values of $\gamma$ can be associated with warming scenarios and larger initial glacier size. **Note:** IMHO, this is completely nonsensical, given that volume/area scaling should not be used on single glaciers and must be used in combination with proper response time scaling.

**Volume evolutions: sensitivity experiments**

While the volume evolution (normalized to initial values) using model (a) and (b) follows the flowline model better for the initial 50 years, the overall difference at the end of the simulation is smallest for model (c) and around 12%. This is to be expected, since it is the most complex model that updates area and allows for a new equilibrium to form. In general, larger mass balance perturbations and smaller initial glaciers are more sensitive to the chosen scaling method. All calculations are performed with the scaling exponent $\gamma = 1.56$. The same experiment performed on glaciers that are initially in a non-steady state (i.e. the respective mass balance is nonzero for some decades prior to the perturbation) yields stronger deviations between different scaling methods and from the flowline model.

The difference in final volume projection using different scaling exponents ($\gamma \in {1.56, 1.89, 2.27, 1.375}$) is about 6% of the initial value. As before, the value increases for larger mass balance perturbations and smaller initial glaciers.

An additional experiment extends the observation period by 200 years, while keeping the climate constant (i.e. at the perturbation level of year 100). As expected, only the glacier using scaling method (c) reaches a new steady state (like the flowline model). 

#### **References:**

[Bahr et al. (2015)]: http://doi.org/10.1002/2014RG000470	"Bahr, D. B., Pfeffer, W. T., and Kaser, G. ( 2015), A review of volume‐area scaling of glaciers. Rev. Geophys., 53, 95– 140. doi:10.1002/2014RG000470."
[Oerlemans (1997)]: http://doi.org/10.3189/S0260305500012489	"Oerlemans, J. (1997). A flowline model for Nigardsbreen, Norway: Projection of future glacier length based on dynamic calibration with the historic record. Annals of Glaciology, 24, 382-389. doi:10.3189/S0260305500012489"

