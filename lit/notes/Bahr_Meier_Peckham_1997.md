### The physical basis of glacier volume-area scaling

**David B. Bahr, Mark F. Meier, Scott D. Peckham (1997) - Journal of Geophysical Research - https://doi.org/10.1029/97JB01696**

**Abstract**:

Ice volumes are known for only a few of the roughly 160,000 glaciers worldwide but are important components of many climate and sea level studies which require water flux estimates. A scaling analysis of the mass and momentum conservation equations shows that glacier volumes can be related by a power law to more easily observed glacier surface areas. The relationship requires four closure choices for the scaling behavior of glacier widths, slopes, side drag and mass balance. Reasonable closures predict a volume-area scaling exponent which is consistent with observations, giving a physical and practical basis for estimating ice volumes. Glacier volume is insensitive to perturbations in the mass balance scaling, but changes in average accumulation area ratios reflect significant changes in the scaling of both mass balance and ice volume.

**Volume Size Scaling:**

In the simplest representation, a glaciers volume can be computed as product of its characteristic width $[w]$, length $[x]$ and ice thickness $[h]$. Hence, the surface area equals $[S] \propto [w][x]$. Contrary to simple geometric objects (like spheres, cubes, ...) , width, length, and thickness of glaciers do not scale identically. This result in the following relationship
$$
\begin{align}
\begin{split}
  [V] &\propto [w][x][h]\\
  	&\propto [S][h]\\
  	&\propto [S]^\gamma
\end{split}
\label{eq:volume_surface_scaling}
\end{align}
$$
whereby the value of $\gamma \simeq 1.36$ is suggested from field data. Given that valley glaciers are not one-dimensional surfaces with widths scaling as constants (e.g., parabolic cross section with $[w] \propto [h]^{1/2}$), it is necessary to include a relation between glacier width and length which allows for all possibilities. Assuming $[w] \propto [x]^q$, where $q$ is a constant, Equation ($\ref{eq:volume_surface_scaling}$) becomes
$$
\begin{align}
\begin{split}
  [V] &\propto [x]^{1+q}[h]\\
  [S] &\propto [x]^{1+q}
\end{split}
\label{eq:volume_surface_width}
\end{align}
$$
**Scaling Analysis:**

Combining the relations above (Eqn. ($\ref{eq:volume_surface_scaling}$) and Eqn. ($\ref{eq:volume_surface_width}$)) we get
$$
\begin{align}
\begin{split}
  [V] &\propto [S]^\gamma \\
	[w][x][h] &\propto [x]^\gamma[w]^\gamma \\
\end{split}
\end{align}
$$
and solving for the ice thickness yields $[h] \propto [x]^\theta$ where $\theta = (q+1)(\gamma-1)$. This implies that the ice thickness $[h]$ must be a power law in glacier length $[x]$, which can be expected from a similarity analysis (see Buckingham Pi theorem).

> The constant $\theta$ is derived explicitly from a scaling analysis of the partial differential equations describing glacier dynamics. (Bahr, Meier and Peckham, 1997).

The scaling analysis is based on ~~the shallow ice~~ equations for mass conservation (continuity), momentum conservation (ice flow) and Glen's constitutive ice creep law, hence
$$
\begin{align}
\begin{split}
	[\dot{b}] \propto \frac{[u_x][h]}{[x]}
\end{split}
\qquad\qquad
\begin{split}
	[u_x] = [A][\rho]^n[g\sin\alpha]^n[h]^{n+1}[F]^n
\end{split}
\end{align}
$$
whereby $\dot b$ is the point mass balance, $u_x$ is the down-glacier component of the velocity, $A$ is Glen's ice creep parameter, $g$ is the gravitational acceleration, $\alpha$ is the slope and $F$ is a shape factor accounting for ice drag on the side walls and the glacier bed. Combining those two relations to eliminate $u_x$ yields
$$
[h]^{n+2} \propto \frac{[\dot b][x]}{[\sin\alpha]^n[F]^n}.\label{eq:thickness_scaling}
$$
To express the ice thickness as function of glacier length, all terms in Eqn. ($\ref{eq:thickness_scaling}$) must be functions of glacier length. Using power laws in the form of
$$
[\sin\alpha] \propto [x]^{-r}, \qquad\qquad [F] \propto [x]^{-f}, \qquad\qquad [\dot b] \propto [x]^{m},
$$
for some constants $r$, $f$ and $m$ and substituting them into Eqn. ($\ref{eq:thickness_scaling}$) yields
$$
[h] \propto [x]^\tfrac{1+m+n(f+r)}{n+2}
$$
and consequentially
$$
[V] \propto [S]^{1+\tfrac{1+m+n(f+r)}{(q+1)(n+2)}}.
$$
**Choices for closure:**

The **width closure** is based on a study of 24'000 Eurasian glaciers and 5'400 Alps glaciers, suggesting that $[w] \propto [x]^{0.6}$, hence $q=0.6$ .

@TODO...

**Closure summary:**

> Although many closures are possible, for valley glaciers $q\simeq 0.6$ (width exponent) and $m\simeq 2$ (mass balance exponent) are both suggested by data. The value $f\simeq 0$ for side drag may also be appropriate. There is no supporting data for the slope closure, but of the reasonable values discussed, only the steep slope scaling $r=0$ is consistent with the other closure choices for $q$, $m$ and $f$. [...] For valley glaciers the closure choices predict $\gamma = 1.375$ [...] in excellent agreement with the available data. (Bahr, Meier and Peckham, 1997).

