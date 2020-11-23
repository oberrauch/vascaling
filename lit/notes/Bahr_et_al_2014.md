### Glacier volume estimation as an ill-posed inversion

**David B. Bahr, W. Tad Pfeffer, Georg Kaser (2014) - Journal of Gaciology - http://doi.org/10.3189/2014JoG14J062**

**Abstract**:

> Estimating a glacier’s volume by inferring properties at depth (e.g. bed topography or basal slip) from properties observed at the surface (e.g. area and slope) creates a calculation instability that grows exponentially with the size of the glacier. [...] Volume/area scaling inherently filters these short wavelengths and automatically eliminates the instability, [numerical inversions can also filter these wavelength]. but when calculating volume, neither the modeling nor the scaling approach offers a fundamental advantage over the other. Both are significantly limited by the inherently ‘ill-posed’ inversion, and even though both provide stable volume solutions, neither can give unique solutions. (Bahr et al., 2014)

**Introduction:**

Glacier volume is the most fundamental geometric property, yet the least known one for the vast majority of all the world's glaciers and ice caps. Direct ice volume or ice thickness measurements (and all other internal properties such subsurface geometry, internal velocities, sliding rates, ...) are difficult to obtain, hence they are estimated from surface properties (i.e., inversion).

> Inversions of this type over-specify data on one boundary (the glacier surface) and under-specify the conditions on another boundary (the glacier bed). An unbalanced placement of boundary conditions is a notorious way to create unstable solutions [...], and for the system of equations that describe a glacier (continuity, force balance and constitutive), placing all of the boundary conditions at the surface can create both a calculation instability and an ‘ill-posed boundary-value problem’. (Bahr et al., 2014)

The unbalanced application of boundary conditions (all on the surface, non on the bottom) results in an ill-posed problem, the type of numerical solution is thereby irrelevant. A problem is **ill-posed** if the solution is either not unique, not stable or does not exist. Given that the glacier has a volume, a solution must exist. Therefore the following investigation focuses on unstable and/or non-unique solutions, both cases imply the existence of a cloud of possible solutions.

However, the ill-posed nature of the problem agrees with our understanding of glaciers. Glacier with similar surface areas can have vastly different ice volumes, based on their geographic location, climatic conditions, ... (non-uniqueness). Additionally, given the glaciers tendency to smooth out small scale basal features, similar surface topographies can arise from vastly different basal topographies. In other words, small changes at the surface will most likely be due to substantially different conditions on the bed (instabilities).

**Identifying reasonable (regularized) solutions:**

Ice thickness inversion is and will be ill-posed, however there are the following three mitigation factors that render scaling and numerical approaches useful nonetheless:

1. Ill-posed does not mean unsolvable. Given that ill-posed problems are well studied in other fields (engineering, geophysics), progress was made in turning ill-posed problems into sets of well-posed problems. The process is called regularization and places reasonable bounds on the solution by assumes parts of it and thereby constraining the inversion.

2. By the law of large numbers, the sum over a substantially large set of glaciers' volume will be a good estimate since the random errors of each glacier are likely to balance out. However, additional uncertainties introduced by other parameters (e.g., Glen's A-parameter, sliding parameters, ...) must be addressed separately.

3. > The ill-posed instability is dependent on the spatial wavelength or resolution of the solution. [...] defining the spatial wavelength (and related spatial frequency) as the length scale of variation of whatever surface properties are being projected downwards in the inversion for distributions at the bed. (Bahr et al., 2014)
   
Whatever that means... as far as I understand it, the instability depends on the spatial resolution (or more abstract the spatial wavelength and associated frequencies). 
   
> The trick is to eliminate the unstable short frequencies that generate most of the errors. Eliminating or smoothing high frequencies has been the strategy of many glacier inversions, and it is another example of regularization (which assumes something about the structure of the solution at the bed) (e.g. Truffer, 2004; Habermann and others, 2012). (Bahr et al., 2014)
   
No matter the type of inversion, errors grow chaotically and exponentially. Thereby, higher frequencies result in a stronger exponential growth. Hence, substantially large spatial wavelength can produce stable solutions.

**Error analysis for the ill-posed volume problem:**

> Any errors in a specified (or calculated) stress at the surface of a glacier will grow exponentially larger when the resulting stress is estimated at the bed of a glacier (Bahr and others, 1994). [...] deviatoric stress errors from ill-posed calculation grow exponentially as $|| \Delta \hat\sigma_{ij}^{'}(k, z_b) || = || \Delta\hat\sigma_{ij}^{'}(k, z_s) ||\ \mathrm{e}^{kr(n)H}$. [...]  $\hat\sigma_{ij}^{'}(k, z_s)$ is the measurement error of a stress at the surface of the glacier, and $\hat\sigma_{ij}^{'}(k, z_b)$ is the resulting inversion calculation error at the bed. The hat notation, $\hat{}$, indicates that the function has been Fourier transformed in the along-glacier direction, $x$ (but not transformed in the vertical direction, $z$); in other words, along-glacier spatial variations in $x$ have been transformed to the spatial frequency domain, $k$. Spatial frequency is related to spatial wavelength by $k = 2\pi/\lambda$, so Eqn ($\ref{eq:exponential_growth}$) indicates that the low-frequency (and long-wavelength) components of the solution have relatively small errors while the high-frequency (and short-wave-length) components of the solution have exponentially larger errors. (Bahr et al., 2014)

Bahr et al. (1994) show that measurement (or calculation) errors (denoted by $\Delta$) in stress at the surface grow exponentially with depth according to
$$
|| \Delta \hat\tau_{ij}(k, z_b) || = || \Delta\hat\tau_{ij}(k, z_s) ||\ \mathrm{e}^{kr(n)H}, \label{eq:exponential_growth}
$$
whereby $\tau_{ij}$ represents a deviatoric stress component, $z_b$ and $z_s$ are the surface and bedrock elevation, respectively, $H = z_s - z_b$ is the ice thickness and $k = 2\pi/\lambda$ represents the wave number corresponding to the wavelength $\lambda$. The constant $r(n)$ depends on the exponent of Glen's flow law, with a typical value of $r(n) \approx 0.6$ for the common assumption of $n = 3$. The hat ($\hat{\square}$) denotes the Fourier-transformation of the along-glacier axis into frequency domain (the vertical axis remains in time domain).

It is therefore paramount to keep the exponent below one, in order to minimize the error growth below one order of magnitude
$$
\begin{align}
\begin{split}
kr(n)H = \frac{2\pi r(n)H}{\lambda} \ll 1 \\
\frac{3.77\ H}{\lambda} \ll 1 \Leftrightarrow \lambda \gg 4\ H
\end{split}
\end{align}
$$
which is the case if the wavelength is approximately four times the glacier thickness.

> For reasonable errors that grow by no more than an order of magnitude, one ice thickness is the shortest practical horizontal wavelength at the bed, in which case $|| \Delta \hat\sigma_{ij}^{'}(k, H) || / || \Delta \hat\sigma_{ij}^{'}(k, z_s) || = 43$. [...] Practically this means keeping the exponent’s value $\ll 1$: $kr(n)H = \frac{2\pi r(n)H}{\lambda} \ll 1$. For $n = 3$, the constant terms become $2\pi r(n) = 3.77$, and the shortest acceptable wavelength will be approximately four times the glacier thickness. (Bahr et al., 2014)
>

**Thickness errors:**

The parallel flow model describes the surface velocity $u_s$ as follows, whereby $u_b$ describes the basal velocity, $A$ and $n$ are the parameters of Glen's flow law, $\tau_b$ describes the basal stress and $H$ the ice thickness:
$$
u_s = u_b + \frac{2A}{n+1}\tau_b^nH\\
$$
Given that the basal stress $\tau_b$ can be approximated by $\tau_b \simeq \rho g \alpha H$, the surface velocity is proportional to the ice thickness to the power of 4 (assuming $n=3$) or more general
$$
u_s - u_b \propto H^{n+1} \qquad\text{and}\qquad H \propto (u_s - u_b)^\frac{1}{1+n}.
$$
This shows that an error in measured surface velocity will result in a related error in ice thickness. Furthermore, an exponentially large error in calculated basal velocity will produce an exponentially large error in ice thickness. The exponent $1/(1+n)$ may dampen the magnitude, but does not change the exponential nature.

Similarly, the linear relation between basal stress and ice thickness $\tau_b \simeq \rho g \alpha H$, results in a 1:1 error propagation
$$
\tau_b + \Delta\tau_b = \rho g \alpha (H + \Delta H).
$$
Solving for the error in ice thickness $\Delta H$, Fourier transforming to the spatial frequency domain and using the exponential relation between deviatoric stress errors (Eq. ($\ref{eq:exponential_growth}$)) results in:
$$
||\Delta \hat H|| = \frac{1}{\rho g \alpha}||\Delta \hat\tau_{x,z}(k,z_s)||\ \mathrm{e}^{2\pi r(n)H/\lambda}
$$
In words, calculating the basal stress from the surface stress results in an exponentially large error, which in term propagates to the ice thickness estimate (for parallel flow). This result can be generalized, using the relation between shear stress and thickness derived from a dimensional analysis (or equivalently from a stretching transformation). ...

**Ill-posed errors from volume/area scaling:**

Volume/area scaling is derived from a dimensional analysis (or stretching transformation) and relates the ice volume to the glacier's surface area by the following power law
$$
V = c\,S^\gamma.\label{eq:vas}
$$
While the scaling exponent is a constant $\gamma = 1.375$, the scaling constant $c$ is a (normal distributed) random variable. Hence, the solution cannot be unique, and the inversion is per definition ill-posed. However, the solution is stable. More precisely, the volume error $\Delta V$ scales linearly with surface error $\Delta S$ as
$$
\Delta V\approx c\gamma\,S^{\gamma-1}\Delta S, \label{eq:vas_error}
$$
 given that for small errors in surface measurements $\Delta S \ll S$ higher order terms of $\Delta S$ vanish
$$
\begin{align}
\begin{split}
V + \Delta V &= c\,(S+\Delta S)^\gamma \\
	&= c\,\sum_{k=0}^\gamma\binom{\gamma}{k}S^{(\gamma-k)}\Delta S^k\\
	&\approx c\,S^\gamma + c\gamma\,S^{\gamma-1}\Delta S.
\end{split}
\end{align}
$$
The stability arises from the fact, that volume/area scaling uses sufficiently large spatial wavelengths to keep the errors irrelevant.

The dimensional analysis produced a whole suite of equation, relating all fundamental continuum parameters with each other. Especially, the surface area can be expressed as a power law of any glacial continuum parameter $P$ as
$$
S = c_p P^{\gamma_p}.\label{eq:surf_scaling}
$$
This allows to translate the volume/area relation into a volume/stress relation by substituting $S = c_\tau {\tau_{ij}}^{\gamma_\tau}$ into Eq ($\ref{eq:vas}$)
$$
\begin{align}
\begin{split}
V &= cS^\gamma\\
	&= c(c_\tau {\tau_{ij}}^{\gamma_\tau})^\gamma \\
	&= c{c_\tau}^\gamma {\tau_{ij}}^{\gamma_\tau\gamma} \\
	&= c_V\, {\tau_{ij}}^{\gamma_V}. \\
\end{split}
\end{align}
$$

The scaling exponents are related such that $\gamma_\tau = \frac{1}{\gamma-1}$ and hence $\gamma_V = \gamma_\tau\gamma = \frac{\gamma}{\gamma-1}.$ An error in (measured or calculated) stress $\Delta \tau_{ij}$ will generate an error in calculated volume $\Delta $V according to $V + \Delta V = c_V(\tau_{ij} + \Delta \tau_{ij})^{\gamma_V},$ whereby for small errors $\Delta \tau_{ij} \ll \tau_{ij}$ higher order terms of $\Delta \tau_{ij}$ are negligible, hence
$$
\Delta V \approx c_V\gamma_V{\tau_{ij}}^{\gamma_V-1}\Delta \tau_{ij} \qquad\text{and analogously}\qquad
\Delta S \approx c_\tau\gamma_\tau{\tau_{ij}}^{\gamma_\tau-1}\Delta \tau_{ij}.
$$
Using the general scaling relation surface area and any glacial continuum parameter $P$ (Eq. ($\ref{eq:surf_scaling}$)), small perturbations $\Delta P \ll P$ results in the following surface error $\Delta S \approx c_P\gamma_PP^{\gamma_P-1}\Delta P$. Using this statement, multiplying by $P$ and rearranging yields
$$
\frac{\Delta S}{S} \propto \frac{\Delta P}{P} \propto \frac{\Delta L}{L} \propto \frac{\Delta V}{V} \propto \frac{\Delta \tau_{ij}}{\tau_{ij}} \propto \ldots
$$
All parameters and their errors are related the same, which shows the necessity for a common length scale $L$. In other words, the length scale (and all other characteristic quantities for that matter) must be the same for all scaling relations and cannot change during the application. The most appropriate characteristic length scale is the glacier length, which can be derived from the  area/length scaling relation
$$
S = b\,L^{q+1},
$$
whereby the scaling parameter $b \approx 1$ and the scaling exponent $q=0.6$ for glaciers. Since the magnitude of errors depends on the used spatial wavelengths, and the minimum resolvable wavelength is twice the shortest length scale $\lambda = 2L$, solving the area/length relation for $L$ and substituting it yields
$$
\lambda = 2L = 2\left(\frac{S}{b}\right)^\frac{1}{q+1}.\label{eq:wavelength_scaling}
$$
Using the exponential relation between deviatoric stress errors (Eq. ($\ref{eq:exponential_growth}$)) and substituting $H=cS^{\gamma-1}$ yields
$$
\begin{align}
\begin{split}

|| \Delta \hat\sigma_{ij}^{'}(z_b) || &= || \hat\sigma_{ij}^{'}(z_s) ||\ \mathrm{e}^{kr(n)H} \\
	&= || \hat\sigma_{ij}^{'}(z_s) ||\ \mathrm{e}^{r(n)H/L} \\
	&= || \hat\sigma_{ij}^{'}(z_s) ||\ \exp\left({r(n)cS^{\gamma-1}}\left(\frac{b}{S}\right)^\frac{1}{q+1} \right)\\
	&= || \hat\sigma_{ij}^{'}(z_s) ||\ \exp\left(r(n)cb^\frac{1}{q+1}\ S^{\gamma-1-\frac{1}{q+1}}\right) \\
	&= || \hat\sigma_{ij}^{'}(z_s) ||\ \exp\left(m\ S^\eta \right)
\end{split}
\end{align}
$$
The exponents take the form $m=r(n)cb^\frac{1}{q+1} = 0.064$ and $\eta = \gamma - 1 - \frac{1}{q+1} = -0.25$ for both glaciers and ice caps. Applying the McLaurin series expansion for $\exp(1/\Lambda) \approx (1 + 1/\Lambda)$, the fundamental (normalized) wavelength associated with scaling is
$$
\Lambda_S = \frac{1}{mS^\eta} \approx 16\,S^{0.25},
$$
and the corresponding fractional error $\delta_S = mS^\eta = 0.064\, S^{-0.25}$. So for all glaciers bigger than $S>1\ \mathrm{km}^2$, stress errors do not increase significantly during the inversion, since the wavelength is big enough $\Lambda_S \gg 1$ and the relative error tends to zero $\delta_s \to 0$.

These findings can be applied to the volume/area scaling relation, showing that Eq. ($\ref{eq:vas_error}$) is equally valid as Fourier transformed statement
$$
\begin{align}
\begin{split}
||\Delta\hat V|| &= c\gamma S^{\gamma-1}||\Delta\hat S||\mathrm{e}^{mS^\gamma} \\
 &\approx c\gamma S^{\gamma-1}||\Delta\hat S|| \qquad\qquad \text{for } S > 1 \ \mathrm{km^2}.
\end{split}
\end{align}
$$


> At wavelengths relevant to scaling, the stress errors do not increase significantly with the calculation depth, and the volume errors do not increase exponentially with the thickness of the glacier. At the long wavelengths inherently used by volume/area and volume/stress scaling, there is no exponential growth of the errors between the surface and the bed of the glacier. The short and problematic wavelengths have been automatically excluded by volume/area scaling. (Bahr et al., 2014)

**Ill-posed errors from numerical inversions:**

As before, the minimum resolvable wavelength is twice the shortest length scale, which for numerical models is the grid size $\mathrm dx$ and hence $\lambda = 2\mathrm dx$. Using this in Eq. ($\ref{eq:exponential_growth}$) and the relation between wavelength and wavenumber yields $|| \Delta \hat\tau_{ij}(k, z_b) || = || \Delta\hat\tau_{ij}(k, z_s) ||\ \exp(\pi r(n)H/\mathrm dx)$. This shows, that for a fixed grid size the error will grow exponentially. In order to limit the error growth (keep exponent around zero), the grid spacing must be adjusted with the glacier size as
$$
\frac{\pi r(n) H}{\mathrm dx} \approx 0 \qquad\Leftrightarrow\qquad \mathrm dx \gg \Lambda_m\pi r(n) H \simeq 18.8H
$$
Hereby, the normalized wavelength $\Lambda_m = 10$ assure a typical order of magnitude estimate. This indicates that the horizontal grid spacing should be about 20 times the ice thickness, which for all but the largest glaciers results in less than 10 grid points.



**A comparison of numerical inversions and scaling:**

The spatial wavelength associated with scaling growths with glacier size according to $\lambda_s = 2\,S^{0.625}$ (Eq. ($\ref{eq:wavelength_scaling}$)), while for numerical inversions the wavelength growth with glacier size according to $\lambda_m = \Lambda_m\, 2\pi r(n)c\ S^{0.375}$ with $\Lambda_m \approx 10$ as measure for one order of magnitude. Hence, numerical model will have higher resolutions than scaling (even after taking care of short-wavelength components) for all glaciers with $S>1\ \mathrm{km^2}$. To be more precise, "the scaling wavelengths $\lambda_s$ will be quadratically larger than the modeling wavelengths $\lambda_m$" (Bahr et al., 2014).
$$
\begin{align}
\begin{split}
\lambda_m &= \Lambda_m\, 2\pi r(n)c\ \left(\frac{\lambda_s}{2}\right)^{(\gamma-1)(q+1)}\\
&\approx {\lambda_s}^{0.6} \text{ (for glaciers)}\\
&\approx {\lambda_s}^{0.5} \text{ (for ice caps)}
\end{split}
\end{align}
$$
But, even though the modeling wavelength is chosen to minimize the errors arising from the ill-posed nature of the inversion, the errors still grow exponentially. As shown above, the fractional error for modeling is around 10%, whereas the fractional error for scaling is only around 1% (for glaciers with $S\approx1\,000\ \mathrm{km^2}$). Hence scaling is an order of magnitude more accurate than modeling $\varepsilon = \delta_m/\delta_s = 0.10/0.01 = 10$. While further limiting the spatial wavelength in the numerical inversion will decrease the relative error $\varepsilon$ compared to scaling, in turn the accuracy will go down.

Other sources of error may have a significant impact on the thickness and/or volume estimate, whereby numerical models have much more potential sources than scaling ...

**Aggregate errors:**

In essence: while the law of large number and the central limit theorem suggest that the sum over randomly distributed and unbiased errors will yield a well defined mean and variance, nothing in the above analysis suggest anything other than an exponential growth. Hence, not accounting for errors from the ill-posed problem may result in meaningless results...

**Conclusions:**

