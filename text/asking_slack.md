# Asking Ben:

1. Scaling parameters (constant and exponent) for the volume/area and volume/length scaling, respectively.
2. Response times or better time scales for length and area adjustment.
3. Iterative process of finding the best historic start area.

## Scaling parameters and response times:

Hallo Ben. Ich bin gerade dabei dabei das Scaling-Model aus eurem Paper ([Marzeion et. al., 2012](https://www.the-cryosphere.net/6/1295/2012/tc-6-1295-2012.html)) im OGGM Framework zu implementieren. Dabei bin ich beim *Dynamischen Model*, sprich dem *area/volume* und *length/volume scaling*, auf einige Unklarheiten gestoßen. Vielleicht kannst du mir dabei helfen?!

Danke inzwischen und Beste Grüße!

### Start length

Laut eurem Paper berechnet ihr die Startlänge wie folgt (2.1.5, Eq. 9, S. 1299):
$$
L_\text{start} = c_L \cdot (A_\text{start})^q
$$
Diese Gleichung ist aber nicht stimmig mit den gegebenen Einheiten, nämlich:
$$
\mathrm[m \stackrel{?}{=} m^{3-q} \cdot (m^2)^q = m^{3-q+2q} = m^{3+q} \neq m]
$$
Laut [Radić et. al., 2008](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/E7C46372D4056A5120E3BF3B77B3FDFE/S002214300020871Xa.pdf/analysis_of_scaling_methods_in_deriving_future_volume_evolutions_of_valley_glaciers.pdf) gilt aber $V = c_L \cdot L^q$, was im Umkehrschluss bedeuten würde
$$
L_\text{start} = (\frac{1}{c_L} V_\text{start})^{1/q}.
$$
Auf welchen Überlegungen basiert die von euch benutze Formel?!

Auch wenn man das Problem mit den Einheiten vernachlässigt, ist auch der Wert der Anfangslänge zu klein. Am Beispiele Hintereisferner, mit einer Anfangsfläche von $A_\text{start} = 8.036\ \mathrm{km^2}$ und einem entsprechendem Anfangsvolumen von $V_\text{start} = 0.596\ \mathrm{km^3}$ ergeben sich folgenden Werte:
$$
L_\text{start, Marzeion} = c_L \cdot (A_\text{start})^q = 1.78\ \mathrm{km}
$$
laut der Formel von [Marzeion et. al., 2012](https://www.the-cryosphere.net/6/1295/2012/tc-6-1295-2012.html), und
$$
L_\text{start, Radić} = (\frac{1}{c_L} V_\text{start})^{1/q} = 4.89\ \mathrm{km}
$$
laut [Radić et. al., 2008](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/E7C46372D4056A5120E3BF3B77B3FDFE/S002214300020871Xa.pdf/analysis_of_scaling_methods_in_deriving_future_volume_evolutions_of_valley_glaciers.pdf). Die im RGI vermerkte Länge des Hintereisferners in 2003 beläuft sich auf $7.18\ \mathrm{km}$. Demzufolge unterschätzen beide Formeln die effektive Länge, wenn auch die letztere deutlich bessere Werte liefert. (Zur Berechnung der Länge wurden folgende Parameter verwendet: $c_L = 0.018\ \mathrm{[km^{3-q}]}$ und $q = 2.2$.)

**Anmerkung: Unterschiedliche Einheiten ergeben unterschiedliche Werte.** Die Nichtlinearität der Gleichung erfordert unterschiedliche Scaling Konstanten je nach verwendeter Einheit. [Radić et. al., 2008](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/E7C46372D4056A5120E3BF3B77B3FDFE/S002214300020871Xa.pdf/analysis_of_scaling_methods_in_deriving_future_volume_evolutions_of_valley_glaciers.pdf) arbeiten mit SI-Einheiten gegeben, sprich $c_L = 4.5507\ \mathrm{[m^{3-q}]}$ und $c_A = 0.191\ \mathrm{[m^{3-2q}]}$. Dies entspricht den Konstanten in Marzeion et. al., 2012 mit  $c_L = 0.018\ \mathrm{[km^{3-q}]}$ und $c_A = 0.034\ \mathrm{[km^{3-2q}]}$, wobei Länge bzw. Fläche in $[\mathrm{km}]$ bzw. $[\mathrm{km^2}]$ gegeben sind.

```python
def scale_volume(area, SI=True):
    """Compute volume using the following area/volume scaling:
    	V = c_a * A^gamma
    If SI is True, the area is supplied in squared meters [m2],
    and in squared kilometers otherwise [km2]."""
    if SI:
        # scaling constant for area in m2
        ca = 0.191
    else:
        # scaling constant for area in km2
        ca = 0.034
    # scaling exponent
    gamma = 1.375
    # compute volume by scaling law
    return ca * area ** gamma

# area and volume in km2 and km3, respectively
area_km2 = 8.036
volume_km3 = scale_volume(area_km2, SI=False)
print('Area of {:.4} km2 corresponds a volume of {:.4} km3'.format(area_km2, volume_km3))

# area and volume in m2 and m3, respectively
area_m2 = 8036000.0
volume_m3 = scale_volume(area_m2)
print('Area of {:.4} m2 corresponds a volume of {:.4} m3'.format(area_m2, volume_m3))
```

Das obige Code Snippet produziert folgenden Output, und $0.5969\ \mathrm{km^3} \approx 0.5963\cdot{}10^9  \ \mathrm{m^3} $.

```python
$> Area of 8.036 km2 corresponds a volume of 0.5969 km3
$> Area of 8.036e+06 m2 corresponds a volume of 5.963e+08 m3
```

### Time scales or response times

Die Zeitskala für Längenänderungen ist gegeben durch das Verhältnis von Gletschervolumen und klimatologischem festen Niederschlag, $\tau_L = V/P(t^*)^\text{solid}_\text{i, clim}$. Dabei verstehe ich nicht genau welcher Wert sich hinter $P(t^*)^\text{solid}_\text{i, clim}$ verbirgt. Bis jetzt habe ich angenommen, es handelt sich um die Mittelwert des festen Jahresniederschlags über die 31-jährige Klimaperiode zentriert um $t^*$. Stimmt das so weit?!

Die nächste Frage bezieht sich auf die Einheiten und damit auch auf die Größenordnung. Gilt die angegeben Formel nur als *best guess* und die Einheiten werden vernachlässigt?!
$$
\tau_L = \frac{V}{P(t^*)^\text{solid}_\text{i, clim}} \quad \stackrel{?}{=} \quad
\mathrm{\left[\frac{km^3}{\frac{mm}{m^2\ a}} = \frac{10^9m^3}{\frac{10^{-3}m}{m^2\ a}} = 10^{12}\ m^4\ a \right]}
$$

Mir ist auch nicht klar, welcher Term in welcher Einheit verwendet werden soll... Und alle *Trial and Error* Versuche ergeben entweder zu geringe Längen-Änderungen (von wenigen Metern auf 100 Jahren) oder numerische Instabilität. Laut [Jóhannesson et. al., 1989](https://www.cambridge.org/core/journals/journal-of-glaciology/article/timescale-for-adjustment-of-glaciers-to-changes-in-mass-balance/B34E832DC10EADCF6DF4220695F34B9A) sollten die Zeitskalen in der Größenordnung zwischen 10 un 100 Jahren liegen. Wenn wir das am Beispiel Hintereisferner durchspielen, ein Volumen von $0.59\ \mathrm{km^2}$ annehmen sowie eine mittleren jährlichen Niederschlag von $1942\ \mathrm{mm\ w.e.}$, ergibt sich
$$
\tau_L = \frac{V}{P(t^*)^\text{solid}_\text{i, clim}} = \frac{0.59 \cdot 10^9\ \mathrm{m^3}}{1.942\ \mathrm{\frac{m}{m^2\ a}}} = 0.307\cdot 10^9\ m^4\ a.
$$
Dieser Wert scheint mir viel zu groß. Aber ich finde meinen Denkfehler nicht...

## 3. Iterative start area seeking process

Hello Ben and all others interested!

**WHAT AM I DOING:**

I'm still working on the implementation of the 'original' volume/area scaling model. So far, the basic structure seems to work and the model produces some sensible results. However, I have issues estimating the glacier’s surface area `A_0` at the beginning of the model integration (i.e. 1850, since I'm working with HistAlp cliamte data). As described in Ben's sea level rise paper, the initial surface area "is estimated by iteratively seeking that surface area in the starting year of the integration that will result in the measured surface area in the year of the measurement."

**WHERE IS THE PROBLEM:**

The main problem is: for some glaciers the model is not able to produce a surface area closed to the measured in 2003 for any initial surface `A_0`. Not even  It seems like the climate forcing is much more important than the initial glacier geometry. Additionally, the larger the initial area the smaller the final glacier area after around 150 years of model evolution. I prepared some plots for different glaciers, to illustrate my problem.

![](../plots/start_area/hintereisferner_overturn.pdf)

The plots show the modeled evolution of glacier area for Hintereisferner over the HistAlp climate period, given different initial values for `A_0`. The dashed line represents the initial surface area which results in todays (or better 2003's) measured surface area, while the dash-dotted line shows the model evolution for a glacier initialized with a historic measurement from 1850. Plots for other (reference) glaciers can be found in the `../plots/start_area/` directory.

**ASK FOR HELP:**

All equations, units, ... are triple checked at this point, while it is still more than possible that I have some bugs in my code. So I'm happy for any suggestions... Thank's so far.

### Answers:

My response/explanation to the first reactions on Slack contains not much useful information. [#science, 2019 April 8th](https://oggm.slack.com/archives/C0L9KNKJN/p1554719026000400)

> ... As [@nchampollion](https://oggm.slack.com/team/U5EV0PW0P) [pointed out], small glaciers seem to work better, even if I have not tried enough different glaciers to find a pattern. But in combination with [@matthias](https://oggm.slack.com/team/U9SV1M0MU)' comment, I now assume that there is a problem with my way of computing the change in terminus elevation
>
> The terminus elevation is just a linear function of glacier length and the initial (i.e. measured) average slope. But I think it is wrong to assume a constant slope, independent of the glacier size. I assume a larger (and therefore longer) glacier to have a lower average slope, especially in alpine terrain considering steep mountains and flat(er) valleys. This could also explain the *turnover*, given too low and too high terminus elevations for too large and too small glaciers, respectively...
> About your point [@julia.eis](https://oggm.slack.com/team/U2KF3F2KE): I didn't occur to me, but as I'm thinking about it now it is strange that the same glacier area evolves differently under the same climate. My explanation is the influence of the time scales on the area and length changes. The time scales depend on the glacier volume which is not necessarily equal for equal areas during the model run. But I definitively have to look into that...

Some good comments came from [@David Parkes](https://oggm.slack.com/team/UBM45M2CX), raising the assumption that the system could be underdamped (i.e., that the calculated response times are too long).

> Hi [@Moritz](https://oggm.slack.com/team/U7FCK7CV7) - this is an implementation of the model from Ben's 2012 paper, yes? As far as I can see, nothing from the graphs suggests there is necessarily a **problem** with the implementation, **but if** there is, my first assumption would be that it's with the **response time**. Because this model has mass balance immediately update the glacier's volume annually, but then adjust the area and length gradually towards an expected equilibrium volume/area ratio over a number of years, the crossover points can be due to the same glacier area with different volumes, where the bulk of the annual area change comes not from that year's mass balance but from the 'stored' response to past years' volume changes. To see quite a significant oscillation and the areas not converging over the timescales you're using suggests the system is **underdamped**, i.e. the response time used to calculate how much the area/length correct towards that expected for the volume is too long - of course there's no physical reason a glacier 'should' be optimally damped but it does raise questions and I know the calculation of response times is a valid lever to tune. Whether that's an error in implementation or in the model I don't know.
>
> The points on the glacier geometry and terminus are entirely reasonable, and are limitations of the model, but the errors due to this should relate to the equilibrium a glacier is tending towards at any given time, rather than explain the 'bounce' you're seeing with large glaciers becoming too small and small glaciers becoming too large.

Here is **[ben.marzeion](https://app.slack.com/team/U0KQG20JU)**'s comment, also he points in the direction of time scales (or response times).

> ... I also think that the different behaviors (connected with getting smaller glaciers when starting with larger ones) is caused by the **time scaling** (i.e., larger disequilibrium between volume and area, and then very sudden disappearence of a glacier that is stretched thin).
>
> Also in my original model, this iteration wasn't always successful, but it was rather the small glaciers that caused problems. Because of this, it was only a very small percentage of ice area in a region where the iteration failed, and I didn't bother to look deeper...

After that I looked into the codebase (Ben original Matlab code) and found what causes the problem.

> As suspected by many of you, the problem arises from a wrong (or different) computation of the terminus elevation. The glacier is initialized with a value for the surface area, from which the volume and length are computed via given scaling laws. After that point my model differs from Ben's code:
>
> While Ben still uses the minimal and maximal elevation from the RGI entity, I scaled the terminus elevation proportional to the new glacier length. This results in strongly positive/negative mass balances for smaller/larger glaciers in the first couple year, which explains the *overturning*.
>
> Now it seems to work just fine, even though I think my first approach (or at least the thinking behind it) is more consistent. Maybe [@ben.marzeion](https://oggm.slack.com/team/U0KQG20JU) still knows what their reasoning was, because obviously it work[s better].

Here is the same plot as above (model evolution for the Hintereisferner) but with the "fixed" code.

![](../plots/start_area/hintereisferner.pdf)

My explanation why my approach is physically more consistent,

> Let's assume the new initial area to be bigger than the RGI value. Given a larger surface are, the glacier gets a new larger volume and a new longer length (through the scaling laws). And a longer glacier will reach further down valley, which is why I scaled the terminus elevation accordingly (Eq. 8 of your paper). Otherwise I'd end up with a bigger glacier confined to the same "elevation band" as the RGI glacier, right?! Or am I missing something?!

and Ben's response:

> No, you are not missing anything - but I thought I do that, too?
>
> Oh no: I agree, I don't do that. But **that is not deliberate, it's a bug**!
>
> This is terrible!
>
> And wow, that's a strong effect. If you do it correctly (which creat[e]s the strange results) and **remove the time scaling** (i.e., adjust area and volume **instantaneously**), does it get better? I think it should, because it should remove the overshoots.

So, I also looked into that (i.e., setting the response time to 1 year):

> Looks like your right. If the length and area adjustments occur instant, the "overturning" problem does not occur. In fact, the model behaves almost the same as without an initial terminus elevation adjustment but with response times in place. Just a bit more noisy.

![](../plots/start_area/hintereisferner_instant.pdf)

And finally, here is [@fmaussion](https://app.slack.com/team/U0KQ07VJT)'s final  

> My suggestion for now will be for Moritz to implement the algorithm as you have it for your 2012 paper (edited). We already made some changes, which means that OGGM-VAS will not be 1 to 1 same as your model, but close.
>
> The other question will require more testing, maybe from you Ben if you have time. Because the convergence plots in the "correct" version are very troubling, especially because all curves converge too fast to a certain glacier (i.e. it reminds me of the "real" OGGM)

## Next question

https://oggm.slack.com/archives/C0L9KNKJN/p1561549522002800