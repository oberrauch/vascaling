### Glacier response to climate perturbations: an accurate linear geometric model

**Gerard H. Roe and Marcia B. Baker (2014), Journal of Glaciology - https://doi.org/10.3189/2014jog14j016**

#### Abstract

Using a numerical **flowline model** solving the **shallow-ice equations** with sliding, the paper investigates the fundamental parameters governing glacier advance and retreat as well as the spectral properties of fluctuations in glacier length. Thereby they find, that

> [...] the time evolution and spectral shape of glacier excursions depend on a single parameter, a time constant determined by the geometrical properties of the glacier. (Roe and Baker, 2014)

The glacier flow as a response to mass balance perturbation can be solved by a three-stage linear model (third order linear differential equation), which captures the following (overlapping) stages:

1. changes in interior thickness
2. changes in terminus flux
3. changes in glacier length

#### Introduction

The authors ... @TODO

#### Summary and Discussion

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