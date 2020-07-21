% -----------------------------------------------------------------------------
% Description of LaTeX Template for a Science Thesis
%
% File:    thesis_template_v2.0.tar
%
% Version: 2.0 (December 2009)
%
% Author:  Alexander Gohm
%          Institut of Meteorology and Geophysics
%          University of Innsbruck
%          alexander.gohm@uibk.ac.at
% -----------------------------------------------------------------------------

(1) CONTENT -------------------------------------------------------------------
This README.txt file is a description of a LaTeX template for a science thesis.
All files of this template are stored in the tar-file thesis_template_v2.0.tar.
On a Linux/Unix platform use the following command to unpack this file:

tar -xvf thesis_template_v2.0.tar

This will create a directory "thesis" that contains the following files:
README.txt             ... This file
thesis.tex             ... Main LaTeX file of the thesis template
thesis.dvi             ... DVI file of the thesis template
thesis.pdf             ... PDF file of the thesis template

The following tex-files are "included" (called) in thesis.tex:
titlepage_master.tex   ... Title page for Master's thesis
titlepage_diploma.tex  ... Title page for Diplomar thesis (optionally called)
dedication.tex         ... Dedication
preface.tex            ... Preface
abstract.tex           ... Abstract
chapter1.tex           ... Chapter 1 (Introduction)
chapter2.tex           ... Chapter 2 (Methodology)
chapter3.tex           ... Chapter 3 (Results)
discussion.tex         ... Discussion
conclusions.tex        ... Conclusion
appendix.tex           ... Appendix A
thanks.tex             ... Acknowledgments
curricv.tex            ... Curriculum Vitae
epilogue.tex           ... Epilogue

Here are some additional files that are called in thesis.tex:
figure_topo.eps        ... EPS figure of topography of target area
figure_theta.eps       ... EPS figure of time series of theta
mybibfile.bib          ... BibTeX file containing citations
ametsoc.bst            ... Bibliography style file
journals-ams.tex       ... Short cuts for AMS journals


(2) COMPILE LATEX FILE --------------------------------------------------------
To compile the main latex file thesis.tex (i.e., to create a dvi file) use the
following four commands:

latex thesis
bibtex thesis
latex thesis
latex thesis

This will create and thus overwrite the existing file thesis.dvi.
To convert the DVI file to a PDF file use the following command:

dvipdfm thesis.dvi

This will create and thus overwrite the existing file thesis.pdf


(3) MODIFY TEMPLATE -----------------------------------------------------------
Use thesis.tex as well as the included tex-files as template for your thesis.
Some hints how to create your own thesis based on this template can be found in
the files thesis.pdf and thesis.tex. Use a text editor or a latex editor
(such as kile on LINUX) to open and modify the LaTeX source code.
