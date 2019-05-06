# Poster Concept

Title: Testing the importance of explicit glacier dynamics for future glacier evolution in the Alps

## Introduction

Even though the *OGGM* is a rather simple model (1.5-dimensional flowlines, shallow ice approximation), the original glacier model used by [Marzeion et. al., 2012]() is even more basic. While the OGGM used the shallow ice approximation to drive the model glacier (and thereby incorporates some physics), the VAS model relies solely on a volume/area and volume/length scaling principle.

- [ ] **What** do I want to do?

  This is my Master Thesis, so I personally don't want to do anything but would much rather go skiing. 😅 But despite that, I'm interested in the behaviour of glacier models. Currently I'm implementing the original volume/area scaling model used by Ben (at least I'm trying to). 

- [ ] **Why** do I do it?

  A more complex model comes with a higher computational cost. For detailed analysis this may be necessary, but not all scientific questions call for such a high degree of accuracy. In addition, alpine glaciers wont advance much in the coming decades - meaning, not much ice physics going on (on a larger scale).

- [ ] **What** do I want to achieve? (Research question)

  In a nutshell: how far can I dumb the dynamic model down, while still producing reasonable (comparable, physical, ... choose one) results on a regional and/or global scale. This includes an closer investigation of the following topics:

  - implementation of glacier model(s) with different levels of complexity in the OGGM framework

  - strength and weaknesses of the VAS model approach
  - regional runs for the alpine region (past and future)

## Methodology

- [ ] **How** do I do it?

  - The mass balance model is quite similar to the one used by the OGGM. The (yearly) mass balance is computed as the difference of solid precipitation fallen onto the glacier surface (mass input) and positive melting temperature at the glacier terminus (energy input) scaled by the temperature sensitivity parameter $\mu^*$. 

    It uses monthly temperature and precipitation data as input, which are scaled to the given glacier elevation. The main difference to the OGGM model is, that it computes a single (average) value of melting temperature and solid precipitation for the entire glacier.

  - The scaling model ...

  - Confusogram of the day...

  - Finding the glacier surface area at the beginning of the model integration ...

- [ ] ...

## Preliminary results

- [ ] Compare mass balance models

- [ ] Compare length/area/volume over HistAlp climate period

  ![](../plots/length.png)

  ![](../plots/area.png)

  ![](../plots/volume.png)

  

- [ ] ...

## Discussion - i.e. problems I'm facing currently

- [ ]  Start area
- [ ] Time scales
- [ ] ...

## Conclusion

Not so much to conclude here, since it's just a concept so far...

