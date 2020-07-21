### Width and length scaling of glaciers

**David B. Bahr (1997) - Journal of Glaciology - https://doi.org/10.3189/s0022143000035164**

**Abstract**:

> An analysis of hundreds of mountain and valley glaciers in the former Soviet Union and the Alps shows that the characteristic glacier widths scale as characteristic glacier lengths raised to an exponent of 0.6. [...] the linear width-thickness realtionship is not incosnistent width parabolic cross-sections [...] (Bahr, 1997)

#### 2. Observations

The analysis of the length-width behavior is based on 303 Eurasian glaciers and 112 Alpine glaciers, withouth longitudinal discontinuities. The attributes of mean widt, maximum length, total area ablation area and mean depth are obtained from a digital inventory of almost 24.5 thousand Eurasian glaciers and ice caps and almost 5.5 thousand glaciers from the European Alps selected based on given accuracy.

#### 3. Analysis

The measurements taken from the inventories of the geometric parameters (mean widt, maximum length, total area ablation area, ablation length and mean depth) are assumed to be the characteristic values for the corresponding glacier, since the excat choice of the characteristic values is not critical (given a scaling constant). Following the Buckingham Pi theorem the characteristic width $[w]$ can be related to the characteristic length $[x]$ using a power law, whereby $q$ is the scaling exponent.
$$
[w] = [x]^q
$$
Most observations suggest a parabolic cross section of the glacier bed (or any other simple polinomial of the order $p$), hence the characteristic width $[w]$ scales with the characteristic depth or thickness $[h]$  as
$$
[w] = \sqrt[p]{[h]} = [h]^{\frac{1}{p}}.
$$
Using the relation between characteristic thickness $[h]$ and characteristic length $[x]$ from [Bahr et al. (1997)][] one obtaines another width-length scaling in form of a power law:
$$
\begin{equation}
\begin{split}
[w] = [x]^{\frac{(m+1)}{p(n+2)}}
\qquad\ldots\qquad\text{for steep slopes}
\\{}
[w] = [x]^{\frac{(m+n+1)}{2p(n+1)}} \qquad\ldots\qquad\text{for shallow slopes}
\end{split}\label{eq:polinomial-shape}
\end{equation}
$$
The exponent $n$ comes from the ice creep relation $\dot{\varepsilon} = A\tau^n$ (Glen's law), and the exponent $m$ comes from the mass balance rate discribed as a function of distance along the glacier $\dot b(x) \simeq -c_mx^m + c_0$.

**Area-length data**

The scaling exponent $q$ is obtained by a linear regression on the log-log plot of mean glacier width and maximum glacier length. Given that the observations of mean width are highly subjective, surface area will be used instead. The characteristic surface area is proportional to the product of characteristic width and characteristic length, hence
$$
[S] = [w][x] = [x]^{1+q}.
$$
This results in a scaling exponent $q=0.61$ for the Eurasian data set (with $R^2 = 0.81$) and $q=0.69$ for the Alpine data set (with $R^2 = 0.89$).

**Volume-area data**

The scaling exponent is definitely positive ($q>0$) and the best fit arises from $q \in [0.6,0.7]$, but even $q=1$ gives reasonable appearing fits (same scaling exponent as for ice caps). Therefore, the scaling exponent $q$ is determined independently using the glacier surface and the glacier volume, given that
$$
[V]\propto[S][h]\propto[S]^\gamma.
$$
The mean thickness and volume data from the inventory seems to be calculated rather than measured for a large number of glaciers, since the regression analysis shows an irritating lack of noise. However, enough independent studies suggest a scaling exponent in the range of $\gamma \in [1.36, 1.40]$. If the creep scaling exponent $n\simeq 3$ is assumed to be know, the width scaling exponent $q$ can be determined for any selected value of the mass balance exponent $m$ using
$$
\begin{equation}
\begin{split}
q = \frac{m+1}{(\gamma-1)(n+2)}-1
\qquad\ldots\qquad\text{for steep slopes, and}
\\{}
q = \frac{m+n+1}{2(\gamma-1)(n+1)}-1
\qquad\ldots\qquad\text{for shallow slopes.}
\end{split}\label{eq:q-from-volume}
\end{equation}
$$
The scaling exponent $m$ can be estimated by measurements of accumulation-area ratio, which can be shown to be

It can be shown that the accumulation area ration is given by
$$
\text{AAR} = \left(\frac{1}{m+1}\right)^\frac{1}{m}.
$$
With a measurement average of 0.577 and 0.580 for the Eurasian and Alpine glaciers, respectively, the exponent $m$ is estimated to be $m=1.99$ and $m=2.04$.

> For $n = 3$ and for the observed values of $\gamma \approx 1.36$ and $m\approx2$ (from the mean AAR), Equations [($\ref{eq:q-from-volume}$)] predict $q = 0.6$ for steep slopes and $q = 1.0$ for shallow slopes [...] good agreement with the area-length data estimates of $q = 0.56$ (Eurasia ) and $q = 0.69$ (Alps ).

**Discussion**

The discussion is mainly concerned with the implicated linear relation between characteristic glacier width and characteristic glacier thickness for steep slopes (as follows from Equations ($\ref{eq:polinomial-shape}$) for $q=0.6$, $n=3$ and $m=2$). 

> Characteristic values can be derived from relationships between variables, but the far more general relationship between variables cannot be derived from the reduced information contained in the characteristic numbers.
>
> The only reasonable explanation, therefore, is that the characteristic *glacier* width and thickness are different from the characteristic width and thickness of the *glacier valley*.

**Conclusion**

> The Eurasian and Alps inventories show that characteristic valley glacier widths are related to characteristic valley glacier lengths by $[w]\propto[x]^q$ with $q\approx0.6$ [...] derived from observations of length-area scaling and [...] supported by observations of volume-area scaling [...] consistend with an accumulation-area ratio of slightly less than two-thirds (as frequently assumed) and comparatively steep (rather than shallow) surface slopes.
>
> [...] implies that characteristic glacier width is linearly related to characteristic glacier thickness [... suggesting] that glaciated and glacierized valleys should be V-shaped [...] changes in channel width with longitudinal distance and variations in the glacier surface profile, glacier and valley with a fairly general class of shape can have parabolic cross-sections but a linear relationship $[w]\propto[h]$.

