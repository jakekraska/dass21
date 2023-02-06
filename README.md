# Graded Response Model Analysis and Computer Adaptive Test Simulation of the DASS-21

This repository includes the data and code required for analysis in the **Graded Response Model Analysis and Computer Adaptive Test Simulation of the DASS-211** manuscript. This manuscript is currently under peer review.

## Abstract

### Background

The Depression Anxiety Stress Scale 21 (DASS-21) is a mental health screening tool with conflicting research on its factor structure. No studies have yet attempted to develop a Computerised Adaptive Test (CAT) version of it. 

### Objective

This study calibrated items for, and simulated, a DASS-21 CAT using a non-clinical sample. 

### Methods

An evaluation sample (n = 580) was used to evaluate the DASS-21 scales via confirmatory factor analysis, Mokken, and Graded Response Modeling. A CAT was simulated with a validation sample (n = 248) and a simulated sample (n = 10,000) to confirm the generalisability of the model developed. 

### Results

A bifactor model, also known as the “quadripartite” model (one general factor with three specific factors) in the context of the DASS-21, displayed good fit. All scales displayed acceptable fit with the Graded Response Model. Simulation of three unidimensional (depression, anxiety, and stress) CATs resulted in an average 17% to 48% reduction in items administered when a reliability of .80 was acceptable. 

### Conclusion

The current study clarifies previous conflicting findings regarding the DASS-21 factor structure and suggests the quadripartite model for the DASS-21 items fits best. IRT modelling suggests that the items measure their respective constructs best between 0 and 3 theta (mild to moderate severity).

## Authors

Jake Kraska wrote the final manuscript and all code for this analysis. Jake Kraska, Karen Bell and Shane Costello designed the study, and collected the data. All authors are affiliated with the Faculty of Education, Monash University, Australia. Jake Kraska is also affiliated with the Krongold Clinic, Monash University.

## Student Project Contributions

Nicole Robinson participated in data analysis of the DASS-21 Stress factor. Paige Macdonald participated in the data analysis of the DASS-21 Anxiety factor. Suganya Chandrasekaran participated in the data analysis of the DASS-21 Depression factor.

## Correspondence

All correspondence regarding this analysis and manuscript can be directed to [Jake Kraska](mailto:jake.kraska@gmail.com)

## Installation and Running

* Clone this repo to your local machine using https://github.com/jakekraska/dass21
* Open analysis.R in your preferred R IDE
* Set the working directory to the source file within the IDE or using `setwd(dirname(rstudioapi::getActiveDocumentContext($path)`
* Install the necessary packages using the `install.packages("PACKAGENAME")` command
    * car version 3.0-10
    * plyr version 1.8.6
    * ggplot2 version 3.3.3
    * psych version 2.1.3
    * knitr version 1.32
    * lavaan version 0.6-8
    * mirt version 1.33.8
    * mirtCAT version 1.10.4
    * mokken version 3.0.6
    * dplyr version 1.0.5
    * tidyr version 1.1.3
    * latticeExtra version 0.6-29
    * irtDemo version 0.1.4
    * parallel version 4.0.5
* The code is commented heavily and minimal SEM, IRT and CAT knowledge should allow understanding of the output. Minimal R knowledge is required to understand the code.