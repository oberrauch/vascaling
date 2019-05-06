# Ask Ben about scaling parameters

Hallo Ben. Ich bin gerade dabei dabei das Scaling-Model aus eurem Paper ([Marzeion et. al., 2012](https://www.the-cryosphere.net/6/1295/2012/tc-6-1295-2012.html)) im OGGM Framework zu implementieren. Dabei bin ich beim *Dynamischen Model*, sprich dem *area/volume* und *length/volume scaling*, auf einige Unklarheiten gestoßen. Vielleicht kannst du mir dabei helfen?!

Danke inzwischen und Beste Grüße!

## Start length

Laut eurem Paper berechnet ihr die Startlänge wie folgt (2.1.5, Eq. 9, S. 1299):
$$
L_\text{start} = c_L \cdot (A_\text{start})^q
$$
Diese Gleichung ist aber nicht stimmig mit den gegebenen Einheiten, nämlich:
$$
\mathrm[m \stackrel{?}{=} m^{3-q} \cdot (m^2)^q = m^{3-q+2q} = m^{3+q} \neq m]
$$
Laut [Radić et. al., 2008](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/E7C46372D4056A5120E3BF3B77B3FDFE/S002214300020871Xa.pdf/analysis_of_scaling_methods_in_deriving_future_volume_evolutions_of_valley_glaciers.pdf) gilt aber $V = c_L \cdot L^q​$, was im Umkehrschluss bedeuten würde
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







## Time scales

Die Zeitskala für Längenänderungen ist gegeben durch das Verhältnis von Gletschervolumen und        klimatologischem festen Niederschlag, $\tau_L = V/P(t^*)^\text{solid}_\text{i, clim}$. Dabei verstehe ich nicht genau welcher Wert sich hinter $P(t^*)^\text{solid}_\text{i, clim}$ verbirgt. Bis jetzt habe ich angenommen, es handelt sich um die Mittelwert des festen Jahresniederschlags über die 31-jährige Klimaperiode zentriert um $t^*$. Stimmt das so weit?!

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

## First results

Anbei noch die ersten Plots, welche die Entwicklung der Länge, der Fläche und des Volumens über die HistAlp Klimaperiode der beiden Modelle vergleichen. Beide Modelle laufen mit der Outline von 2003 und ohne SpinUp, dementsprechend handelt es sich mehr um qualitative Plots. Das grundlegende Entwicklung scheint durchaus vergleichbar, die absoluten Werte und dabei vor allem die Amplituden sind in meiner Implementation jedoch  um mehrere Größenordnungen zu klein. Dabei habe ich mich entschieden, die Zeitskalen für die Längen- ($\tau_A$) sowie die Flächenänderung ($\tau_L$) konstant auf 10 Jahre zu fixieren...

Die Ergebnisse des Scaling-Models sind jeweils in blau gezeichnet und beziehen sich auf die linke y-Achse, während die Ergebnisse des OGGMs in orange dargestellt sind und sich auf die rechts y-Achse beziehen.

![Length comparison](../plots/length.png) 

 ![Area comparison](../plots/area.png) 

![Volume comparison](../plots/volume.png) 





