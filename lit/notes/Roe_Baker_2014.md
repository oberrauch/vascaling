## Glacier response to climate perturbations: an accurate linear geometric model

**Gerard H. Roe and Marcia B. Baker (2014), Journal of Glaciology - https://doi.org/10.3189/2014jog14j016**

### Abstract

Using a numerical **flowline model** solving the **shallow-ice equations** with sliding, the paper investigates the fundamental parameters governing glacier advance and retreat as well as the spectral properties of fluctuations in glacier length. Thereby they find, that

> [...] the time evolution and spectral shape of glacier excursions depend on a single parameter, a time constant determined by the geometrical properties of the glacier. (Roe and Baker, 2014)

The glacier flow as a response to mass balance perturbation can be solved by a three-stage linear model (third order linear differential equation), which captures the following (overlapping) stages:

1. changes in interior thickness
2. changes in terminus flux
3. changes in glacier length

### Summary and Discussion

> At low frequencies and long timescales, a glacier responds in quasi-equilibrium with the climate forcing. The glacier’s length is dictated by the geometry it must attain to achieve a quasi-balanced mass budget, and the ice dynamics primarily affects glacier thickness. (Roe and Baker, 2014)

What does that mean?! @TODO: Try to understand this...

Ice dynamics are more important for mass balance changes with high frequencies and short time scales. Thereby is the mass distribution much more efficient in the interior than near the glacier terminus. The geometry adjustment can be described by three overlapping stages: (I) changes in interior thickness drive (II) changes in terminus ice flux, which in turn (III) drive the changes in glacier length. (Roe and Baker, 2014)

> On a constant slope, the temporal evolution of the length response is controlled by a single timescale (as function of glacier geometry, melt factor and lapse rate as already derived by Jóhannesson et. al. (1989)) (Roe and Baker, 2014)



> [...] cannot be matched by any one-stage model of the form of Eqn (2). The glacier has more short-term persistence, and is more damped at high frequencies, than Eqn (2) would imply. (Roe and Baker, 2014)

A one stage model driven by the linearized mass conservation equation
$$
\frac{\mathrm{d}L'}{\mathrm{d}t} + \frac{L'}{\tau} = \alpha T' + \beta P',
$$
is not able to reproduce the characteristic sigmoidal evolution of glacier length in response to a step in climate. Hereby $L'$ is the glacier length perturbation linearized about the mean length $\bar{L}$, $T'$ and $P'$ are the melt season temperature and annual mean precipitation, respectively. The coefficients are given by
$$
\alpha = -\frac{\mu A_{T>0}}{wH},\\
\beta = \frac{A_\text{tot}}{wH},\\
\tau = \frac{wH}{\mu \Gamma\tan\phi A_\text{abl}},
$$
whereby $H$ is the ice thickness, $w$ is the terminus width,  $\tan\phi$ the basal slope, $A_\text{tot}$, $A_{T>0}$ and $A_\text{abl}$ the total area, the melting area (melt season temperature T above freezing), and the ablation area (below ELA), respectively. $\Gamma$ is the atmospheric lapse rate and $\mu$ the melt factor relating ablation to melt season temperature. 

> This timescale must not be interpreted as an e-folding timescale [...] if $\tau$ is estimated from the time it takes for a numerical glacier model to reach $(1 - 1/e)$ of its new equilibrium length in response to a step change, the inferred timescale will be about twice that of the actual one.

### 2. One-stage model

A linearization of the mass conservation equation about the mean length $\bar L$ gives a model for glacier length perturbations $L'$ driven by changes in melt season temperature $T'$ and annual mean precipitation $P'$. For simple geometries, i.e., constant ice thickness $H$, terminus width $w$ and bed slope $\phi$ this yields
$$
\renewcommand{dfrac}[2]{\frac{\mathrm d\, #1}{\mathrm d\, #2}}
\dfrac{L'}{t} + \frac{1}{\tau}L' = \alpha T' + \beta P'.
$$
The coefficients are given by
$$
\alpha = -\frac{\mu A_{T>0}}{wH},\\
\beta = \frac{A_\text{tot}}{wH},\\
\tau = \frac{wH}{\mu \Gamma\tan\phi A_\text{abl}},
$$
whereby $A_\text{tot}$, $A_{T>0}$ and $A_\text{abl}$ are the total area, the melting area (melt season temperature $T$ above freezing), and the ablation area (below ELA), respectively. $\Gamma$ is the atmospheric lapse rate and $\mu$ the melt factor relating ablation to melt season temperature.

This model is advantageous to prior and similar approaches (e.g., [Jóhanneson et al. (1989)][]). The model is derived directly from the mass conservation equation, it distinguishes between temperature and precipitation as forcings and the model parameters are related to glacier geometries. However, changes in the mass balance are immediately converted to a tendency on the terminus length (Roe and Baker, 2014). Hence the name, stemming from the single dynamical stage.

Simple solutions:

- Under no climate anomalies ($T' = 0$, $P' = 0$) an initial length perturbation $L'(0)$ will decay exponentially as $L'(t) = \exp(-t/\tau)$. Hereby, $\tau$ can be seen as the e-folding time scale. This can also be seen in the auto correlation function $\mathrm{ACF}(t) = \exp(-t/\tau)$.
- If a glacier in equilibrium $L'(t=0) = 0$ is subjected to a step change in temperature $\Delta T$ and precipitation $\Delta P$, the solution for the length change is $L'(t) = \Delta L (1-\exp(-t/\tau))$. Hereby, $\Delta L = \tau(\alpha\Delta T + \beta\Delta P)$ is the equilibrium length.
- Applying trends of temperature $T' = \dot T t$ and precipitation $P' = \dot P t$ to a glacier in equilibrium, the length changes as $L'(t) = \tau(\alpha\dot T + \beta\dot P) (t-\tau(1-\exp(-t/\tau)))$. For $t\gg\tau$, this solutions asymptotes to a straight line with slope $\tau(\alpha\dot T + \beta\dot P)$ and intersect $\tau$ on the time axis.
- A discretization into time steps of year $\Delta t = 1\text{ year}$ yields $L'_t = (1-\Delta t/\tau)L'_{t-\Delta t} + \alpha \Delta t T'_t + \beta \Delta t P'_t$. This allows to investigate the model glaciers response to natural random climate variability. Normal distributed and uncorrelated year-to-year fluctuations in temperature $T'$ and precipitation $P'$, with standard deviations of $\sigma_T$ and $\sigma_P$ result in a glacier length fluctuation with variance $\sigma_L^2 = 0.5\tau\Delta t (\alpha^2 \sigma_T^2 + \beta^2 \sigma_P^2)$.

**Intermediate summary:** The time scale $\tau$ occurs in all above solutions, which makes it a central number to estimate a glaciers response to a variety of forcings (step change, climate trend, natural variability). In other words, the glaciers behavior to the above scenarios go hand in hand. A glacier with a higher natural variability will also react stronger to climate change.

- [ ] Power spectrum for discrete length equation, complex frequency response function, ... ???

#### Performance of a one stage model

The one stage model is compared to a flowline model on a constant slope using standard numerical techniques.

The one stage model responds to quickly to a step in accumulation $\Delta P \pm 0.5\ \mathrm{ma^{-1}}$, when compared to the flowline model. While the flowline model shows the classic sigmoidal length evolution, the one stage model shows exponential decay and asymptotic growth, respectively. The e-folding time for the one stage model is around $\tau_\text{one stage} \simeq 6.7\ \mathrm{years}$, while for the flowline model it is more than twice that $\tau_\text{flowline} \approx 15\ \mathrm{years}$. However, after 20 years both models converge and the equilibrium lengths agree to within 5%. Similar behavior is found for a variety of different step functions (cp. [Roe (2011)][]).



