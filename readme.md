# Master thesis 

This directory/repository contains all files related to my Master's Thesis "Testing the importance of explicit glacier dynamics for mountain glacier change projections".

For this work I re-implemented the volume/area scaling model by [Marzeion et al. (2012)][Marzeion et al. (2012)] into the Open Global Glacier Model ([OGGM](https://oggm.org/)) framework. This allowed to perform a set of different experiments comparing the OGGM flowline model to the scaling model under controlled boundary conditions.

The model implementation can be found under [github.com/OGGM/oggm-vas](https://github.com/OGGM/oggm-vas).

## Content

- `code`: Contains all the code (mainly python scripts) used to run the model and/or plot results. The code for the volume/area scaling model itself has a dedicated OGGM repository.
- `data`: The directory contains produced results stored to file(s).
- `defensio`: LaTex environment containing the presentation for the Master's Thesis defence.
- `flowchart`: Latex environment containing different TikZ schematics.
- `icu-meeting`: LaTex environment containing the presentation for the ICU meeting in October 2020.
- `lit`: Contains relevant papers.
- `matlab_model_marzeion`: The original MatLab model code.
- `notebooks`: Contains all Jupyter Notebooks.
- `plots`: Contains all produced plots and graphics.
- `poster`: LaTex environment for the poster presented at the 2019 AGM in Innsbruck.
- `text`: Contains all write-ups (mainly markdown files).
- `thesis`: LaTex environment containing the thesis.
- `working_directories`: Contains all OGGM working directories.

## References:

[Marzeion et al. (2012)]: https://doi.org/10.5194/tc-6-1295-2012	"Marzeion, B., Jarosch, A. H., and Hofer, M.: Past and future sea-level change from the surface mass balance of glaciers, The Cryosphere, 6, 1295–1322, https://doi.org/10.5194/tc-6-1295-2012, 2012."