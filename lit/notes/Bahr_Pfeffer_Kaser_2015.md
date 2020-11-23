### A review of volume-area scaling of glaciers

**David B. Bahr, W. Tad Pfeffer, and Georg Kaser (2015)**

> Bahr,D. B., W. T. Pfeffer, andG. Kaser (2015), A review of volume-area scaling of glaciers, Rev. Geophys., 53, 95–140, https://doi.org/10.1002/2014RG000470.

#### Abstract:

> The application of the theory is not entirely straightforward, [...] many of the recently published results contain analyses that are in conflict with the theory as originally described by [Bahr et al. (1997)][Bahr et al. (1997)]. In this review we describe the general theory of scaling for glaciers in full three-dimensional detail without simplifications, including an improved derivation [...] (Bahr et al., 2015)

Glacier area volume scaling is based on physical principles governing glacier flow (originally derived by Bahr et al. (1997)). Despite that, the usage in recent studies is often in conflict with the underlying theory.

Beyond the derivation of the scaling laws by two different starting points, the following topics are discussed:

- "common misconceptions of the theory" (Bahr et al., 2015)
- "potential future developments in power law scaling" (Bahr et al., 2015)
- "the relationship between power law scaling and other modeling approaches" (Bahr et al., 2015)
- "the advantages and limitations of scaling techniques" (Bahr et al., 2015)

#### 10. Conclusions:

> Volume-area scaling and related power laws have been derived by two different methods. No simplifications are necessary in either approach. The derivations do not assume plane strain, shallow ice, or steady state conditions. A single closure condition is required, and available data for three possible choices of the closure condition all give the same result. Each of these closure conditions imply the others and are self-consistent. (Bahr et al., 2015)

The volume-area scaling power law has a solid physical basis, the derivation does not rely on any simplifications. Two separate derivations and three different choices for the single closure condition all give the same results and imply self consistency.

> For consistency with the underlying theory, we suggest the following guidelines for practical applications of volume-area scaling [...]:
>
> 1. Fix the scaling exponent to the theoretical constant γ = 1.375. [...]
> 2. Apply volume-area scaling to collections of many glaciers, not to individuals. Treat the resulting set of volumes as a probability distribution. [...]
> 3. Volume-area scaling should only be applied to an individual glacier [...]
> 4. Include time dependence by using response time scaling [...]
> 5. Avoid applying volume-area scaling to glacier complexes. [...]
> 6. Compare scaled volumes to numerically modeled volumes [...] only if the model is using the same
>    closure conditions.
>
> (Bahr et al., 2015)

The above mentioned simple rules ensure a sound application of volume-area scaling theory.

#### 2. Fundamental Dimensional Principles

The glacier volume $V$ can be estimated by the glacier surface $S$ using the power law $V = c S^\gamma$, whereby the scaling constant $c$ is a random variable and the scaling exponent $\gamma$ is a constant for a given geometric class of glaciers (especially valley glaciers vs. radially symmetric ice caps). (Bahr et al., 2015)

> [...] whereby parameters $c$ and $\gamma$ were originally determined as empirical constants [Macheret et al., 1988; Chen and Ohmura,1990; Zhuravlev, 1988]. Additional data were later compiled by Meier and Bahr [1996], Bahr et al. [1997], Cogley [2010], and Grinsted [2013]. (Bahr et al., 2015)

Even though the scaling relation was originally developed empirically, the theoretical (physical) justification was established by [Bahr et al. (1997)][] and [Bahr (1997)][]. This allows to apply the scaling relation under circumstance different from those used to establish it, say different climatic/environmental conditions, significantly smaller and/or larger glaciers, equilibrium vs. non equilibrium conditions.

The main pitfalls of using the volume-area scaling relation are summarized. While the scaling constant $c$ may vary, the scaling exponent $\gamma$ is a constant (in space and time). The scaling relation is determined, validated and should therefore be used solely on samples of glaciers spanning a large range of sizes. Additionally, the scaling relation is not valid for parts or individual branches of glaciers, as well as glacier complexes spanning over flow divides and draining into multiple outlets.

The following paragraphs describe the general nature of scaling relations with examples from fluid dynamics (Reynolds number, Froud number, Euler number) and find analogies to the glacier volume are scaling laws.

> Volume-area scaling and response time scaling are both well-known glaciological relationships, but one cannot be considered separately from the other. Instead, a real glacier is characterized and modeled by both relationships simultaneously, and there is no need to artificially assign time dependence to the volume-area exponent.

The scaling exponent $\gamma$ is a constant and not time dependent. However, the volume-area scaling must be considered in combination with response time scaling.

**"Good" papers:**  The following list shows studies that are using volume-area scaling relations largely consistent with it's original theory in order to:

- estimate aggregate ice volumes [Radić and Hock, 2010; Grinsted, 2013]
- predict future changes in ice mass [Meier et al., 2007; Leclercq et al., 2011; Bahr et al., 2009; Mernild et al., 2013; Radić and Hock, 2011; Radić et al., 2013; Raper and Braithwaite, 2006; Marzeion et al., 2012]
- sensitivity of volume-area scaling to perturbations [Van de Wal and Wild, 2001; Slangen and Van de Wal, 2011]
- volume response time [Jóhannesson et al., 1989; Bahr et al., 1998; Pfeffer et al., 1998; Raper and Braithwaite, 2009].

**"Bad" papers** that are in direct conflict with the underlying theory because of

- includes improper time dependence [Arendt et al., 2006; Huss and Farinotti, 2012; Möller and Schneider, 2010]
- inappropriate scaling parameters [Grinsted, 2013; Hagg et al., 2013; Möller and Schneider, 2010]
- assumptions of steady state and/or shallow ice [Huss and Farinotti, 2012; Möller and Schneider,2010; Meehl et al., 2007]
- applications to individual glaciers without regard to accuracy [Agrawal and Tayal,2013; Schneeberger et al., 2003]
- applications to portions of glaciers or to undifferentiated glacier complexes [Grinsted, 2013; Raper and
  Braithwaite, 2006]
- a general rejection of volume-area scaling based on inappropriate statistical arguments [Haeberli et al., 2007]

**2.1 Basic derivations**

While glacier volume can be scaled with a number of quantities (such as area, thickness, length, velocity, ...), surface area is the most easily measurable quantity. To obtain volume $V$ with dimension $\mathrm{L}^3$, the surface $S$ with dimension $\mathrm{L}^2$ must be multiplied with one additional quantity with dimension $\mathrm{L}$. For this glacier thickness is an obvious choice, whereby the centerline thickness $h$ must be adjusted with a shape factor $F$ to account for lateral drag which scales as a ratio of width to thickness $F\propto\frac{w}{h}$. Assuming that the surface area is a product of length $l$ and width $w$ and given that glacier width scales with glacier length as $w\propto l^q$ with $q\approx 0.6$, the glacier volume scales with surface area as
$$
\begin{align}
\begin{split}
V\ \mathrm{[L^3]} &\propto \text{Surface}\ \mathrm{[L^2]}\cdot\text{Thickness}\ \mathrm{[L]}&\qquad \big|\ \text{account for lateral drag}\\
	&\propto S\cdot Fh &\qquad \big|\ \text{shape factor scales with witdth/thickness ratio}\\
	&\propto wl\cdot \frac{w}{h}h &\qquad \big|\ \text{simplify}\\
	&\propto w^2l &\qquad \big|\ \text{width scales with length to the power of }q\\
	&\propto l^{2q+1} &\qquad \big|\ \text{substitute }S\propto wl = l^{q+1}\\
	&\propto S^\frac{2q+1}{q+1} = S^\gamma = S^{1.375} \\
	
\end{split}
\end{align}
$$
This abbreviated derivation relies on the physical intuition for the choice of $Fh$ as additional length dimension and on the geometric closure condition $q=0.6$ (which applies on average to all glaciers). In the following thorough derivation shows that the choice of $h$ is inevitable (through directional analysis) and that the result is consistent with other closure conditions based on mass balance scaling or the accumulation area ratio.

**2.2 An outline of the completer derivation**

As the heading may suggest, this paragraph contains an outline of the following derivation(s) and briefly described the content of each of the following sections:

3. Examines the relevant continuum equations to find all the relevant variables
4. Uses the Buckingham Pi Theorem to derive the full set of dimensionless variables. The additionally performed directional analysis refines the dimensionless parameters and assures that the final results are independent of the underlying flow equations.
5. Applying a stretching symmetry transformation to the continuum equations yields the same results as in Section 4, using a completely independent pathway. Additionally it shows the source of each dimensionless parameter allowing to trace scaling relations to their origin in the continuum equation. This section is of special importance to numerical modelers, since it outlines the derivation of scaling relations using a different (simplified) set of equations. This is crucial to allow for a "apples-to-apples" comparison when comparing models using modified equations (say shallow ice approximation, plane strain, ...) to volume-area scaling.
6. "Summarizes the full set of dimensionless parameters and outlines how to choose and use characteristic values." (Bahr et al., 2015)
7. Shows a complete derivation of the volume area scaling exponent $\gamma$ and constant $c$, and shows that only one closure condition is necessary.

#### 3. Continuum Equations

The set of continuum equations describing the ice flow of glaciers account for mass balance and mass conservation via the continuity equation and force balance and momentum conservation via the equation of motion. Both equations are linked by the constitutive equation, relating force to deformation or more precisely stress to strain rate. The inclusion of an equation accounting for energy conservation adds additional non-dimensional parameters, which do not alter the scaling relations. Hence, it is omitted...

The continuity equation describes the conservation of mass as
$$
\renewcommand{dpar}[2]{\frac{\partial#1}{\partial#2}}
-\dpar{\rho}{t} = \nabla\cdot(\rho\vec{u}) + \dot\mu,\label{eq:continuity}
$$
whereby $\vec{u}$ describes the velocity field and $\dot\mu$ rate of mass gain per unit volume. If ice is assumed to be incompressible, Equation ($\ref{eq:continuity}$) can be integrated over glacier thickness $h$ (as difference between surface and bed elevation $h = h_s - h_b$) yielding
$$
\renewcommand{dd}{\mathrm{d}}
\dpar{h}{t} = \dot b - \dpar{}{x}\int_{h_b}^{h_s}u_x\dd z - \dpar{}{x}\int_{h_b}^{h_s}u_y\dd z.\label{eq:continuity-integrated}
$$
Hereby, the mass balance rate $\dot b$ explicitly includes basal, englacial, and surface terms (even though the first two are often negligible).

The equation of motion includes conservation of momentum as
$$
\renewcommand{dD}[2]{\frac{\mathrm{D}#1}{\mathrm{D}#2}}
\rho\dD{u_i}{t} = \dpar{\sigma_{ii}}{x_i} + \dpar{\sigma_{ij}}{x_j} + \dpar{\sigma_{ik}}{x_k} + \rho g_i,\label{eq:motion}
$$
whereby $i,j,k \in \{x,y,z\},\ i\neq j\neq k$. $\mathrm{D}/\mathrm{D}t$ is the material derivative, which is often assumed to be zero since acceleration is negligible for glaciers and Equation ($\ref{eq:motion}$) becomes the static stress equilibrium.

Using Glen's flow law as constitutive equation to relate stress to strain rate is the only assumption made during the scaling analysis. For $i,j\in \{x,y,z\}$ holds
$$
\dot \varepsilon_{ij} = A\tau_E^{(n-1)}\tau_{ij}.
$$
The effective stress $\tau_E$ is the second invariant of the deviatoric stress tensor and is defined as 
$$
\tau_E^2 = \frac{1}{2}[\tau^2_{xx} + \tau^2_{yy} + \tau^2_{zz}] + \tau^2_{xy} + \tau^2_{xz} + \tau^2_{yz},
$$
whereby the deviatoric stress tensor describes the deviations of the stress from an isotropic state as
$$
\tau_{ii} = \sigma_{ii} - \frac{1}{3}(\sigma_{xx} + \sigma_{yy} + \sigma_{zz}), \\
\tau_{ij} = \sigma_{ij}.
$$
The strain rates are defined as
$$
\dot \varepsilon_{ij} = \frac{1}{2}\left(\dpar{u_i}{x_j} + \dpar{u_j}{x_i}\right).
$$
The specific geometries, surface mass balance, sliding velocity, ... are defined by boundary conditions. Since those boundary conditions do not alter the dimensions of the relevant variables, they have no effect on the following dimensional analysis.

#### 4. Dimensional Analysis and Buckingham Pi

**4.1 Dimensional Analysis**

The continuum equations are based on the following 18 fundamental variables: density $\rho$, time $t$, length along all three dimensions $\vec{x}$, velocity in all three dimensions $\vec{u}$, gravitational acceleration in all three dimensions $\vec{g}$, the six unique components of the symmetric stress tensor $\sigma_{ij}$ for $i,j \in \{x,y,z\}$ and the ice creep parameter $A$. The ice thickness $h$ is indirectly included in the $z$ coordinate, and the mass balance rate $\dot b$ is expressed by the more general vertical velocity $u_z$. The following table lists all variables and their corresponding fundamental dimensions (as exponents e.g., $\rho = \mathbf{[M/L]}^3$).

|                     | $\rho$ | $t$  | $\vec x$ | $\vec u$ | $\vec g$ | $\sigma_{ij}$ |  $A$   |
| :------------------ | :----: | :--: | :------: | :------: | :------: | :-----------: | :----: |
| Length $\mathbf{L}$ |  $-3$  | $0$  |   $1$    |   $1$    |   $1$    |     $-1$      |  $n$   |
| Time $\mathbf{T}$   |  $0$   | $1$  |   $0$    |   $-1$   |   $-2$   |     $-2$      | $2n-1$ |
| Mass $\mathbf{M}$   |  $1$   | $0$  |   $0$    |   $0$    |   $0$    |      $1$      |  $-n$  |

Following the Buckingham Pi Theorem, we have 18 variables with 3 fundamental dimensions, hence there must be 15 dimensionless parameters describing the glacier flow. By selecting three (since we have three fundamental dimensions) variables that span all fundamental dimensions and are linearly independent of each other (here $\rho$, $t$ and $x$), the dimensionless parameter for every other variable $\theta$ can be constructed as
$$
\Pi_\theta = \rho^at^bx^c\theta
$$
for some exponents $a$, $b$ and $c$. As an example, consider the vertical velocity $u_x$ with dimensions $\mathbf{L/T}$. The corresponding $\Pi$-group must be the following
$$
\begin{align}
\begin{split}
	\Pi_{u_x} &= \rho^at^bx^cu_x\\
		&= \rho^0t^{1}x^{-1}u_x\\
		&= \frac{tu_x}{x}
\end{split}
\end{align}
$$
All other non dimensional parameters are constructed analogously, additional new and different parameters can be constructed by substitution e.g.,
$$
\frac{\Pi_{u_x}}{\Pi_{u_y}} = \frac{\frac{tu_x}{x}}{\frac{tu_y}{x}} = \frac{u_x}{u_y}.
$$
This new parameter is obviously not independent of the others, and according to the Buckingham PI theorem in this case there can only be 15 independent parameters (no matter how they are constructed). The final set of $\Pi$-groups obtained after some simplifications and substitutions is
$$
\begin{align}
\begin{split}
\Pi_1 &= \frac{x}{y} \qquad\qquad \Pi_2 &= \frac{x}{z}\\
\Pi_3 &= \frac{u_x}{u_y} \qquad\qquad \Pi_4 &= \frac{u_x}{u_z} \qquad\qquad \Pi_5 &= \frac{tu_x}{x}\\
\Pi_6 &= \frac{g_x}{g_z} \qquad\qquad \Pi_7 &= \frac{g_y}{g_z} \qquad\qquad \Pi_8 &= \frac{t^2g_z}{z}\\
\Pi_9 &= \frac{A\rho^n g_x^n z^{n+1}}{u_x}\\
\Pi_{10} &= \frac{\rho g_x z}{\sigma_{xz}} \qquad\qquad \Pi_{11} &= \frac{\sigma_{xz}}{\sigma_{xx}} \qquad\qquad \Pi_{12} &= \frac{\sigma_{xz}}{\sigma_{xy}}\\
\Pi_{13} &= \frac{\sigma_{xz}}{\sigma_{yy}} \qquad\qquad \Pi_{14} &= \frac{\sigma_{xz}}{\sigma_{yz}} \qquad\qquad \Pi_{15} &= \frac{\sigma_{xz}}{\sigma_{zz}}
\end{split}
\end{align}
$$
The stresses might be considered redundant since they are implied by the constitutive equation, hence the set of non dimensional parameters could be reduced to $\Pi_1$ through $\Pi_9$. The addition of an energy conservation equation introduces new variables and dimensions and adds new dimensionless parameters, however it does not change the above derivations. The same holds true for the introduction of basal sliding laws and other boundary equations.

**4.2 Directional Analysis**

The directional analysis first introduced by Donal B. Siano ([Siano (1985a)][] and [Siano (1985b)][]) adds a coordinate direction or orientation to each variable with fundamental dimension of length in form of a unit vector pointing in the given direction. The variable $x$ with dimension $\mathbf{L}$ and orientation $1_x$ can be represented by the vector $\mathbf{L}1_x$ The orientation of areas is described by their respective normal vector, so for example the product $xy$ has dimension $\mathbf{L}^2$ and orientation $1_x1_y = 1_z$. Generally, the products of orientations are defined by a *Klein* group with products $1_i1_j = 1_k$ for $i\neq j$, identity $1_i1_i = 1$ and inverse $1_i^{-1} = 1_i$. The extension to three dimensions shows that volumes have no directional orientation, given that $V = xyz$ with dimension $\mathbf{L}^3$ and orientation $1_x1_y1_z = 1_z1_z = 1$.

Applying this concept to the above derived dimensionless parameters, the resulting orientation of $\Pi_1$ is along the vertical axis, since
$$
\Pi_1 = \frac{x}{y} = \frac{\mathbf{L}1_x}{\mathbf{L}1_y} = \frac{1_x}{1_y} = {1_x}{1_y} = 1_z.
$$
This is obviously not directionally consistent, which calls for the introduction of dimensionless angles. Analogous to an area, the orientation of an angle is normal to the plane in which it lies. So an angle $\phi_{xy}$ in the $xy$ plane (and all its odd powers ${\phi_{xy}}^k$ for $k \in \mathbb{N}_\text{odd}$) have orientation $1_z$. Introducing such an angle in the above equation via
$$
\Pi_1{\phi_{xy}}^k = \frac{x}{y} \qquad
\left[1_z = \frac{\mathbf{L}1_x}{\mathbf{L}1_y}\right]
$$
ensures dimensional and orientational consistency (as long as $k$ is an odd integer). The variables $x$ and $y$ give magnitudes along their respective axis, hence forming the two legs of a right-angled triangle. Consequentially, the ratio $x/y$ equals the tangent $\tan \phi$, which can be expanded as a series with only odd powers of $\phi$. The value of $k$ cannot be determined by the dimensional or directional analysis, but must be based on physical knowledge (here $k=1$ is the obvious choice).

Angles are introduced to all dimensionless parameters, except $\Pi_5$, $\Pi_8$, $\Pi_9$, $\Pi_{10}$, resulting in
$$
\begin{align}
\begin{split}
&\Pi_1{\phi_{xy}}^k = \frac{x}{y} \qquad\qquad &\Pi_2{\phi_{xz}}^k = \frac{x}{z}\\
&\Pi_3{\phi_{xy}}^k = \frac{u_x}{u_y} \qquad\qquad &\Pi_4{\phi_{xz}}^k = \frac{u_x}{u_z} \qquad\qquad &\Pi_5 = \frac{tu_x}{x}\\
&\Pi_6{\phi_{xz}}^k = \frac{g_x}{g_z} \qquad\qquad &\Pi_7{\phi_{xz}}^k = \frac{g_y}{g_z} \qquad\qquad &\Pi_8 = \frac{t^2g_z}{z}\\
&\Pi_9 = \frac{A\rho^n g_x^n z^{n+1}}{u_x}\\
&\Pi_{10} = \frac{\rho g_x z}{\sigma_{xz}} \qquad\qquad &\Pi_{11}{\phi_{xz}}^k = \frac{\sigma_{xz}}{\sigma_{xx}} \qquad\qquad &\Pi_{12}{\phi_{yz}}^k = \frac{\sigma_{xz}}{\sigma_{xy}}\\
&\Pi_{13}{\phi_{xz}}^k = \frac{\sigma_{xz}}{\sigma_{yy}} \qquad\qquad &\Pi_{14}{\phi_{xy}}^k = \frac{\sigma_{xz}}{\sigma_{yz}} \qquad\qquad &\Pi_{15}{\phi_{xz}}^k = \frac{\sigma_{xz}}{\sigma_{zz}}
\end{split}
\end{align}
$$
This allows to rewrite the dimensionless parameter $\Pi_9$ and solving for $u_x$ as
$$
u_x = \frac{{\Pi_6}^n}{\Pi_9}A(\rho g_z\phi_{xz})^nz^{n+1}
$$
 using the directional consistent form of $\Pi_6$, which introduces the glacier slope $\phi_{xz}$. This equation is equivalent to the equation for the surface velocity $u_s$ for parallel flow
$$
u_s = \frac{2A}{n+1}\tau_b^nH = \frac{2A}{n+1}(\rho g \alpha)^nH^{n+1}.
$$
Hereby, $H$ describes the ice thickness, and $\tau_b \simeq \rho g \alpha H$ the basal/driving stress (cf. Cuffey and Paterson, Chap. 8.3.1).

**4.3 Discussion**

...

#### 5. Stretching analysis of continuum equations

**5.1 The stretching transformation and continuity equation**

The stretching transformation entails a multiplication of each variable $X$ with a scalar factor $\lambda^{k_X}$, whereby $\lambda$ is the same constant for all variables and $k_X$ is a constant stretching exponent specific to $X$. The rescaled quantity is indicated by a bar, hence $\overline X = \lambda^{k_X}X$. Notice that the scaling factor can be represented as dimensionless number by $\lambda^{k_X} = \overline X / X$. Thereby, the equations must be preserved for the stretched variables. This forces some constraints to the scaling constants, for example the topographic variables $h$, $h_s$ and $h_b$ must always rescale by the same amount, since 
$$
\begin{equation}
\begin{split}
h = h_s - h_b \quad &\Rightarrow \quad \overline h = \overline h_s - \overline h_b, \\
\lambda^{k_h} h = \lambda^{k_{h_s}} h_s - \lambda^{k_{h_b}} hb \quad&\Leftrightarrow\quad k_h = k_{h_s} = k_{h_b}.
\end{split}
\end{equation}
$$
Transforming a derivative results in the subtraction of the associated exponents, for example
$$
\dpar{\overline h}{\overline t} = \lambda^{k_h-k_t}\dpar{h}{t}.
$$
Consequentially, transforming an integral results in the addition of the associated exponents, for example
$$
\dpar{}{\overline x}\int^{\overline h_s}_{\overline h_b} \overline u_x \dd\overline z = \lambda^{k_{u_x}-k_{x}+k_{h}}\dpar{}{x}\int^{h_s}_{h_b} u_x \dd z.
$$
Given that $z$ is a placeholder (or dummy variable) for the integration limits, it must scale identically to them. Hence, the scaling exponent $k_z$ is equal to $k_{h} = k_{h_s} = k_{h_b}$ and is therefor not part of the equation. The stretched form of the integrated continuity equation Eq. ($\ref{eq:continuity-integrated}$) reads
$$
\lambda^{k_h-k_t}\dpar{h}{t} = \lambda^{k_{\dot b}}\dot b - \lambda^{k_{u_x}-k_{x}+k_{h}}\dpar{}{x}\int^{h_s}_{h_b} u_x \dd z - \lambda^{k_{u_y}-k_{y}+k_{h}}\dpar{}{y}\int^{h_s}_{h_b} u_y \dd z.
$$
As discussed before, the stretched form must be equal to the original equation. This is only the case, if all $\lambda$-terms can be factored out, hence
$$
\lambda^{k_h - k_t} = \lambda^{k_{\dot b}} = \lambda^{k_{u_x} - k_x + k_h} = \lambda^{k_{u_y} - k_y + k_h}, \\
k_h - k_t = k_{\dot b} = k_{u_x} - k_x + k_h = k_{u_y} - k_y + k_h.
$$
With this, we know how the variables of the balance equation must be scaled in respect to each other. Representing the scaling factors as non-dimensional parameters (e.g., $\lambda^{\dot b} = \overline {\dot b}/{\dot b}$) gives six pairings of ratios, and separating stretched from original variables results in the following six dimensionless numbers
$$
\begin{align}
\begin{split}
	\Pi_{1,c} &= \frac{\dot b t}{h} \\
	\Pi_{2,c} &= \frac{u_x t}{x} \\
	\Pi_{3,c} &= \frac{u_y t}{y} \\
	\Pi_{4,c} &= \frac{\dot b x}{u_x h} \\
	\Pi_{5,c} &= \frac{\dot b y}{u_y h}\\
	\Pi_{6,c} &= \frac{u_y x}{u_x y}
\end{split}
\end{align}
$$
The index $c$ indicates that it was derived from the continuity equation. Additionally, the three trivial relationships $\Pi_{7,c} = \frac{h}{h_s}$, $\Pi_{8,c} = \frac{h}{h_b}$ and $\Pi_{9,c} = \frac{h_s}{h_b}$ can be derived (again, the variable $z$ does not appear outside the integral). Stretching symmetries do not guarantee a set of unique parameters. In our case only $\Pi_{1,c}$, $\Pi_{2,c}$, $\Pi_{3,c}$ and $\Pi_{7,c}$ are unique, which can easily seen above e.g., $\Pi_{4,c} = \Pi_{1,c}/\Pi_{2,c}$.

- [ ] Stretching symmetry of differential form of the continuity equation Eq. ($\ref{eq:continuity}$)

**5.2 Stretching analysis of the equation of motion**

The stretching analysis of the equation of motion Eq. ($\ref{eq:motion}$) follows the same steps as before:

1. Rescale every single variable e.g., $\overline x_i = \lambda^{k_{x_i}}x_i$

2. Establish the rescaled equation (expand the material derivative)
   $$
   \dpar{\overline \sigma_{ii}}{\overline x_i} + \dpar{\overline \sigma_{ij}}{\overline x_j} + \dpar{\overline \sigma_{ik}}{\overline x_k} + \overline \rho \overline g_i = \overline \rho \dD{\overline u_i}{\overline t}
   $$

3. Substitute all stretched variables with the respective product of original variable and scaling factor, whereby derivations are represented as subtractions of the corresponding scaling exponents
   $$
   \lambda^{k_{\sigma_{ii}} - k_{x_i}}\dpar{\sigma_{ii}}{x_i} + \ldots + \lambda^{k_{\rho} + k_{g_i}}\rho g_i = \lambda^{k_{\rho} + k_{u_i} - k_{t}}\rho\dpar{u_{i}}{t} + \lambda^{k_{\rho} + k_{u_x} + k_{u_i} - k_{x}}\rho u_x\dpar{u_{i}}{x} + \ldots
   $$
   
4. The stretched equation must be equal to the original equation, which is the case if all $\lambda$-terms can be factored out
   $$
   \lambda^{k_{\sigma_{ii}} - k_{x_i}} = \ldots = \lambda^{k_{\rho} + k_{g_i}} = \lambda^{k_{\rho} + k_{u_i} - k_{t}} = \lambda^{k_{\rho} + k_{u_x} + k_{u_i} - k_{x}} = \ldots
   $$
   
5. Represent the scaling factors as non-dimensional parameters (e.g., $\lambda^{k_\rho} = \overline \rho / \rho$) resulting in dimensionless parameters

The resulting set of unique parameters is for all $i,j,k \in \{x,y,z\}$
$$
\begin{align}
\begin{split}
	\Pi_{1,\sigma} &= \frac{x_i}{x_j}\frac{\sigma_{ij}}{\sigma_{jk}} \\
	\Pi_{2,\sigma} &= \frac{\rho g_i x_j}{\sigma_{ij}} \\
	\Pi_{3,\sigma} &= \frac{t g_i}{u_i}\\
	\Pi_{4,\sigma} &= \frac{t u_i}{x_i}
\end{split}
\end{align}
$$
whereby, the index $\sigma$ denotes the derivation from the equation of motion. From the combination of $\Pi_{1,\sigma}$ and $\Pi_{2,\sigma}$ yields
$$
\Pi_{5,\sigma} = \frac{x_i}{x_j}\frac{g_j}{g_i}.
$$
The following section shows the derivation of geometric similarity from the constitutive equation, hence also gravity scales equally in all directions. Additionally, since the material derivative has the same dimensional structure as the continuity equation. Hence, $\Pi_{4,\sigma}$ is equal to $\Pi_{c}$—rendering it superfluous—and therefore the statement of continuity is helpful for our intuition but unnecessary for a dimensional analysis.

**Stretching analysis of the constitutive equation**

After substituting the definition of effective stress and deviatoric stress into the flow law and expanding the terms it reads
$$
\varepsilon_{ij} = A \left[\frac{1}{3}\left(\sigma_{xx}^2+\sigma_{yy}^2+\sigma_{zz}^2-\sigma_{xx}\sigma_{yy} - \sigma_{xx}\sigma_{zz} - \sigma_{yy}\sigma_{zz}
\right) +\sigma_{xy}^2 +\sigma_{yz}^2 \sigma_{xz}^2\right]^{\frac{n-1}{2}}\tau_{ij}.
$$
To avoid added complexities the assumption $n=3$ is made, the full expansion can be derived using a Taylor series.

And again, the stretching analysis of the constitutive equation follows the same steps ... yielding a total of 286 dimensionless parameters. As before, most parameters are not unique, i.e. can be constructed using other parameters. There are only two groups of unique parameters, whereby the index $G$ indicates the origin in Glen's flow law.

Most terms in the above equations involve only stresses and no velocities, and each of those terms gives a non-dimensional parameter of the form
$$
\Pi_{1,G} = \frac{\sigma_{ij}}{\sigma_{jk}},
$$
whereby none, some, of all of $i$, $j$ and $k$ can be the same. All stresses must scale identically, describing the similarity of forces in static equilibrium. The other terms include stresses and velocities, resulting in the dimensionless parameter of the form
$$
\Pi_{2,G} = \frac{u_i}{Ax_j\sigma_{ij}^n}.
$$
All other parameters follow from combinations of $\Pi_{1,G}$ and $\Pi_{2,G}$, most important the relations of geometric similarity and static similarity
$$
\Pi_{3,G} = \frac{x_i}{x_j} \qquad\qquad \Pi_{4,G} = \frac{u_i}{u_j}.
$$
If the glacier geometry is stretched in one direction, it must stretch equally in the other two dimensions. The same holds true for velocities and forces.

#### 6. Using dimensionless parameters

**6.1 Summary of dimensionless parameters by physical origins**

The stretching transformation gives the same set of 15 dimensionless parameters as derived using the Buckingham Pi theorem, with the added benefit of knowing the respective physical origin.

Geometric similarity, kinematic similarity and static similarity follow from the constitutive equation $\Pi_1'$, $\Pi_2'$ and $\Pi_3'$, respectively. Gravitational similarity follows from the equation of motion combined with geometric similarity as $\Pi_4'$. The material derivative gives a relation about dynamics and gravity $\Pi_5'$ as well as a response time relationship $\Pi_6'$ (can also be derived from the equation of continuity). The equation of motion gives a stress relationship $\Pi_7'$ and a velocity relationship $\Pi_8'$ when combined with the constitutive equation.
$$
\begin{equation}
\begin{split}
\Pi_1' = \frac{x_i}{x_j}\\
\Pi_2' = \frac{u_i}{u_j}\\
\Pi_3' = \frac{\sigma_{ij}}{\sigma_{jk}}\\
\Pi_4' = \frac{g_i}{g_j}\\
\Pi_5' = \frac{t^2g_i}{x_i}\\
\Pi_6' = \frac{tu_i}{x_i}\\
\Pi_7' = \frac{x_i\rho g_i}{\sigma_{ij}}\\
\Pi_8' = \frac{A(\rho g_i)^nx_j^{n+1}}{u_i}\\
\end{split}
\end{equation}
$$
**6.2 Characteristic values**

The choice of characteristic quantities may seem somewhat arbitrary, especially since there is no single correct choice. In a nutshell, "all of the variables in the dimensionless relations are characteristics and should be assigned values relevant to the particular problem" (Bahr et al., 2015). In this case, values describing the average length, average width and average thickness will give the best estimate of the total ice volume. Furthermore, the actual surface area is the obvious choice for the characteristic area.

**6.3 Choosing dimensionless parameters**

Similarly as above, the dimensionless parameters are chosen with relevance to the particular problem. However, to avoid mistakes a Buckingham Pi analysis could be performed. Thereout arises a particularity. The problem of volume-area scaling does not include time as a relevant quantity. However, volume-area scaling does not assume any steady-state conditions and can be applied to transient state glaciers. In this case, the additional response time scaling must be applied additionally.

#### 7. Volume-area scaling from dimensionless parameters

Even though volume area scaling is the most common application of dimensional analysis in glaciology, the same approach can be used to determine relationships between any other continuum variables (e.g., thickness and horizontal velocity). In more detail, [Bahr (1997)][] showed that any scaling relationship can be derived from the volume-area scaling.

**7.1 Volume-area scaling for glaciers**

> A volume-area scaling relationship can be derived by combining the dimensionless parameters given above and by selecting appropriate characteristic values.

Combining the dimensionless parameters $\Pi_2$, $\Pi_4$ and $\Pi_9$ yields
$$
\frac{\Pi_4\Pi_9}{\Pi_2} = \frac{A\rho^n g_x^n z^{n+2}}{u_z x},
$$
which is valid for any orientation. The expression is independent of any angle, which make a slope closure unnecessary. The following characteristic values and closure conditions are chosen:

- glacier length $l$ as characteristic length $x$
- terminus mass balance $\dot b$ as characteristic vertical velocity $u_z$
- average glacier width $w$ as characteristic width
- average ice thickness $h$ as characteristic thickness $z$
- width/length closure $w = c_q\,l^q$ with $q\approx0.6$ suggested by data ("the value of $q$ reflects the relationship between the glacier itself and the underlying topography" Bahr et al., (2015))
- mass balance/length closure $\dot b = c_m\,l^m$ with $m\approx2$ suggested by data

Substitution gives
$$
\begin{align}
\begin{split}
\frac{\Pi_4\Pi_9}{\Pi_2} &= \frac{A\rho^n g_x^n h^{n+2}}{\dot b\, l}\\
 &= \frac{A\rho^n g_x^n h^{n+2}}{c_m\, l^{m+1}}\\
\end{split}
\end{align}
$$
and solving for average ice thickness $h$ gives
$$
h = \left(\frac{\Pi_4\Pi_9}{\Pi_2}\frac{c_m}{A\rho^n g_x^n}\right)^\frac{1}{n+2}\, l^\frac{m+1}{n+2}.
$$


Volume can be estimated as the product of characteristic ...

**7.2 Volume-area scaling for ice caps**

**7.3 Alternative closure conditions**



#### 8. Discussion

- (8.1) Volume-Area Scaling is derived without assumptions of plane strain, steady state, perfect plasticity, or shallow-ice approximations, rather uses the full set of continuum equations for the derivation of the scaling parameters. Two different and independent techniques (dimensional analysis and stretching symmetries) arrive at the same result. The value for the scaling exponent is completely determined by physics, up to the choice of a single closure condition, which is generated from data. However, all three possible choices for glaciers (width-length scaling exponent, mass balance scaling exponent, equilibrium value for the accumulation area ratio) yield the same results and are internally consistent (both in theory and data).

- (8.2) 

- (8.3)

- (8.5) Volume-area scaling applies to populations of glaciers, and not to individual glaciers, since the worldwide mean value of the scaling constant $c=0.034$ km^3-2γ^ is determined from the (roughly) normal distribution may not apply to a specific glacier. The fractional error in glacier volume estimate is directly proportional to the fractional error in scaling constant
  $$
  \frac{\Delta V}{V} = \frac{\Delta c}{c}.
  $$
  Considering the normal distribution and a on standard deviation error of the scaling constant $c$, a volume error of ~34% is easily possible. On the other hand, the normal distribution and the law of large numbers makes the volume area scaling method more robust for the estimation of the aggregate volume of an ensemble of glaciers. Even if the error of $c$ for individual glaciers is large, the mean value is a reasonable choice for al glaciers in the summation.

  The scaling constant $c$ can be eliminated altogether when estimating global (or regional) changes in
  aggregated ice volume, given that the fractional change in glaciers volume $p_V$ and today's volume $V$ are independent (today's glacier knows nothing about tomorrow's change in climate).

- 

#### References:

[Bahr et al. (1997)]: http://doi.org/10.1029/97JB01696	"Bahr, D. B., M. F. Meier, and S. D. Peckham (1997), The physical basis of glacier volume-area scaling, J. Geophys. Res., 102(B9), 20'355–20'362"
[Bahr (1997)]: http://doi.org/10.1029/97WR00824	"Bahr, D. B. (1997), Global distribution of glacier properties: A stochastic scaling paradigm, Water Resour. Res., 33(7), 1669-1679"

[Siano (1985a)]: http://doi.org/10.1016/0016-0032(85)90031-6	"Siano, D. B. (1985a), Orientational analysis—A supplement to dimensional analysis—I, J. Franklin Inst. Eng. Appl. Math., 320(6), 267–283."
[Siano (1985b) ]: http://doi.org/10.1016/0016-0032(85)90032-8	"Orientational analysis, tensor analysis and the group properties of the SI supplementary units—II, J. Franklin Inst. Eng. Appl. Math., 320(6), 285–302"

