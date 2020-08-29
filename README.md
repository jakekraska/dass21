# Computer Adaptive Testing with the DASS-21

This repository includes the data and code required for analysis in the **Computer Adaptive Testing with the DASS-21** manuscript. This manuscript is currently under peer review.

## Abstract

**Background:** The Depression Anxiety Stress Scale 21 (DASS-21) is a mental health screening tool. Initial investigations into whether the DASS-21 items fit into Item Response Theory (IRT) models have yet attempted to develop a Computerised Adaptive Test (CAT) version of it. **Objective:** To calibrate items for, and simulate, a DASS-21 CAT. Method: An evaluation sample (n = 580) was used to evaluate the DASS-21 scales via CFA, Mokken, and IRT (Graded Response Model). A validation sample (n = 248) and a simulated sample (n = 10,000) was utilised in an IRT analysis and CAT to confirm the generalisability of the model developed during the evaluation phase. **Results:** Confirmatory Factor Analysis comparisons demonstrated that the data fit well with a bifactor model, also known as the “quadripartite” model as demonstrated by Szabó (2010), Gomez (2013) and Henry and Crawford (2005). All scales met the assumptions of IRT and displayed acceptable fit with the Graded Response Model. Multidimensional IRT analysis was attempted but failed to converge, likely due to sample size, thus leaving Multidimensional CAT unattainable. Simulation of three unidimensional DASS-21 CATs resulted in an average 17% to 48% reduction in items administered when a reliability of .8 was acceptable. **Conclusions:** The current study supports the quadripartite model for the DASS-21 items. IRT modelling suggests that the items measure their respective constructs best between 0 and 3 theta, which could be considered mild to moderate. Limitations and potential future research are discussed.

## Authors and Contributions

Jake Kraska wrote the final manuscript and all code for this analysis.

Jake Kraska, Karen Bell and Shane Costello designed the study, collected the data and supervised the other contributors in their Graduate Diploma of Professional Psychology.

Nicole Robinson participated in data analysis of the DASS-21 Stress factor and the first draft of the manuscript. Paige Macdonald participated in the data analysis of the DASS-21 Anxiety factor and the first draft of the manuscript. Suganya Chandrasekaran participated in the data analysis of the DASS-21 Depression factor and the first draft of the manuscript.

All authors and contributors are affiliated with the Faculty of Education, Monash University, Australia. Jake Kraska is also affiliated with the Krongold Clinic, Monash University.

## Correspondence

All correspondence regarding this analysis and manuscript can be directed to [Jake Kraska](mailto:jake.kraska@monash.edu)

## Installation and Running

* Clone this repo to your local machine using https://github.com/jakekraska/dass21
* Open analysis.R in your preferred R IDE
* Set the working directory to the source file within the IDE or using `setwd(dirname(rstudioapi::getActiveDocumentContext($path`
* Install the necessary packages using the `install.packages("PACKAGENAME")` command
    * plyr version 1.8.6
    * ggplot2 version 3.3.2
    * psych version 2.0.7
    * knitr version 1.29
    * lavaan version 0.6-7
    * mirt version 1.32.1
    * mirtCAT version 1.10
    * mokken version 3.0.2
    * dplyr version 1.0.1
    * tidyr version 1.1.1
    * latticeExtra version 0.6-29
    * irtDemo version 0.1.4
    * parallel version 4.0.2
* The code is commented heavily and minimal SEM, IRT and CAT knowledge should allow understanding of the output. Minimal R knowledge is required to understand the code.