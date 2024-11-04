# STATA-SIMUL

This project relates to files for a simulation study to better understand bias of test-negative design COVID-19 vaccine effectiveness studies at priamry care level that can be due to SARS-CoV-2 self-testing.
These files are associated with the article: The potential bias introduced into COVID-19 vaccine effectiveness studies at primary care level due to the availability of SARS-CoV-2 tests in the general population, which will be submitted to preprint on 8 November 2024.

## Description

In the above-mentioned article, we outline our hypothesis around the mechanism of a potential bias in the COVID-19 VE estimates obtained from primary care TND studies. We describe a simulation study, to help understand the magnitude and direction of potential bias and suggest analytical solutions to eliminate this bias.

These are the files behind the simulation study, including data generation and estimation of bias.

Please see the article for details on methods.

## Getting Started

### Prerequisites

All files created in Stata 16. You will need good computational power and a reasonably high file storage capacity for running the scripts.

### Executing the scripts
Please put all files in the same working directory. You run the "simulation self test datagen.do" first, in which "simulation self test datagen_include.do" is nested. Then you can run "simulation self test estimated bias.do" for the estimated bias.

## Authors

Esther Kissling and Gilles Desve. Please contact e.kissling@epiconcept.fr in case of questions.

## Version History

v3.3
New version 1 November 2024
Updated estimated bias scripts and review of within-script comments.

v3.2
New version 2 August 2024
Restructured 12 September 2024 by Gilles Desvé  => v3.2 
use of binomial for binary variables instead of uniform

## License

This project is licensed under CC BY-NC.

## Acknowledgments

Many thanks to Charlotte Laniece and Baltazar Nunes for their numerous inputs into these scripts.
