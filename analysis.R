# CFA Analysis, IRT Analysis and CAT Simulation of the DASS21
# Original code developed by Jake Kraska (Monash University)
# Built on R 4.0.2 with R Studio 1.3.1073
# requires file "data" - a file with 828 observations of DASS21, SEX, AGE, ZIP, EMPLOY and QUALIFICATIONS
# set the working directory to source file directory manually

#---                                        ---#
################# Load Packages #################
#---                                        ---#

require(car) # version 3.0-10
require(plyr) # version 1.8.6
require(ggplot2) # version 3.3.3
require(psych) # version 2.1.3
require(knitr) # version 1.32
require(lavaan) # version 0.6-8
require(mirt) # version 1.33.8
require(mirtCAT) # version 1.10.4
require(mokken) # version 3.0.6
require(dplyr) # version 1.0.5
require(tidyr) # version 1.1.3
require(latticeExtra) # version 0.6-29
require(irtDemo) # version 0.1.4
require(parallel) # version 4.0.5

#---                                       ---#
################# Set Options #################
#---                                       ---#

set.seed(2019) # set the seed so that randomised outcomes are the same across analyses
options(max.print = 3000) # enable the console to print more lines

#---                                     ---#
################# Functions  #################
#---                                     ---#

calculate_sem <- function(rxx,sd) {sd * sqrt(1-rxx)} # creates a function that calculates sem

cat_simulation <- function(model, sem, pattern = NULL) {
  
  if (is.null(pattern) | !exists("pattern")) { # if there is no provided pattern, then simulate a theta distribution
    n <- 10000 # simulate n participants
  } else { # if there is a provided pattern, then get theta distribution
    n <- nrow(pattern)
  } # simulate n
  
  
  if (is.null(pattern) | !exists("pattern")) { # if there is no provided pattern, then simulate a theta distribution
    Theta <- matrix(rnorm(n = n, mean = 0, sd = 1)) 
    } else { # if there is a provided pattern, then get theta distribution
      Theta <- matrix(fscores(object = model, method = "EAP", response.pattern = pattern, append_response.pattern = FALSE)[,1])
      } # Simulate theta
  
  if (is.null(pattern) | !exists("pattern")) { 
    pattern <- generate_pattern(mo = model, Theta = Theta) 
    } else { # if there is a provided pattern, then generate theta and create appropriate object for simulation
      pattern <- data.matrix(pattern)
      attr(pattern, "Theta") <- matrix(fscores(object = model, method = "EAP", response.pattern = pattern, append_response.pattern = FALSE)[,1])
      } ### Generate response pattern 
  
  cat_design <- list(min_SEM = sem) ### the design factors of the simulation
  cl <- makeCluster(detectCores()) ### for parallel processing
  cat <- mirtCAT(mo = model, local_pattern = pattern, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = cat_design, cl = cl) # run
  stopCluster(cl) # stop the cluster
  est_Theta <- laply(cat, function(y) y$thetas) # list of theta
  average_items <- mean(plyr::laply(cat, function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(Theta[,1], est_Theta) ### correlation between true theta and estimated theta
  bias <- mean(Theta[,1] - est_Theta) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((Theta[,1] - est_Theta)^2))
  results <- data.frame("sem" = round(sem, 3), "n" = n, "averageItems" = round(average_items, 3), 
                        "correlation" = round(correlation, 3), "bias" = round(bias, 3), "rmsd" = round(rmsd, 3)) # put the results into a dataframe
  return(results)
}

#---                                          ---#
################# Simulation Data #################
#---                                          ---#

reliability <- c(seq(from = 1.0, to = .50, by = -.1)) # create levels of reliability
sem <- calculate_sem(reliability,1) # calculate SEM associated with each level of reliability
Theta_groups <- seq(-3,3,.6) # create theta groups for comparison of performance

#---                                      ---#
################# Import Data #################
#---                                      ---#

data <- read.csv("data.csv", stringsAsFactors = FALSE) # load the data

#---                                       ---#
################# Split Data #################
#---                                      ---#

sample.split <- 0.70 # set the percentage of the data that will be used for analysis
ind <- sample(c(rep(TRUE,ceiling(nrow(data)*sample.split)),rep(FALSE,floor(nrow(data)*(1-sample.split)))))
data <- cbind(data,ind) # identify which rows will be kept for analysis

#---                                          ---#
################# Rename Columns #################
#---                                          ---#

names(data)[names(data) == "SEX"] <- "Gender" #change the SEX column to Gender
names(data)[names(data) == "EMPLOY1"] <- "Employment" # Change the EMPLOY1 column to Employment
names(data)[names(data) == "QUAL_GRP"] <- "Qualification" # Change the QUAL_GRP column to Qualification
names(data)[names(data) == "AGE"] <- "Age" # Change the AGE column to Age
names(data)[names(data) == "ind"] <- "Sample" # Change the ind column to Sample

#---                                             ---#
################# Recode Data Types #################
#---                                            ---#

data$Gender <- as.character(data$Gender) # recode the variable as type "character"
data$Employment <- as.character(data$Employment) # recode the variable as type "character"
data$Qualification <- as.character(data$Qualification) # recode the variable as type "character"
data$Sample <- as.character(data$Sample) # recode the variable as type "character"

#---                                                      ---#
################# Remove Unnecessary Columns #################
#---                                                     ---#

data$ï..ID <- NULL # remove variable
data$V1 <- NULL # remove variable
data$ZIP <- NULL # remove variable

## Recode Demographics ##

data$Gender <- recode(data$Gender, 
                      "1" = "Male",
                      "2" = "Female")

data$Employment <- recode(data$Employment,
                          "1" = "Student", # Full Time Student
                          "2" = "Employed", # Part Time or Full Time Employment
                          "3" = "Unemployed", # Not Currently In Workforce
                          "4" = "Retired") # Retired

data$Qualification <- recode(data$Qualification,
                             "1" = "Secondary School Incomplete", # high school incomplete 
                             "2" = "TAFE", # TAFE
                             "3" = "Secondary", # completed high school
                             "4" = "Undergraduate", # bachelors
                             "5" = "Honours", # fourth year programs, postgraduate diploma
                             "6" = "Masters/PhD") # PhD, masters

data$Sample <- recode(data$Sample,
                      "TRUE" = "evaluation",
                      "FALSE" = "validation")

data$AgeGroup <- cut(data$Age, breaks = c(18,30,40,50,60,70,80,90), right = FALSE)

data$AgeGroup <- recode(data$AgeGroup,
                        "[18,30)" = "18 to 29",
                        "[30,40)" = "30 to 39",
                        "[40,50)" = "40 to 49",
                        "[50,60)" = "50 to 59",
                        "[60,70)" = "60 to 69",
                        "[70,80)" = "70 to 79",
                        "[80,90)" = "80 to 90")

#---                                            ---#
################# Calculate Totals #################
#---                                            ---#

# calculate a total score

data <- mutate(
  data,
  dep = DEP1 + DEP2 + DEP3 + DEP4 + DEP5 + DEP6 + DEP7,
  anx = ANX1 + ANX2 + ANX3 + ANX4 + ANX5 + ANX6 + ANX7,
  str = STR1 + STR2 + STR3 + STR4 + STR5 + STR6 + STR7,
  tot = dep + anx + str
)

#---                                                ---#
################# Setup Items and Data #################
#---                                                ---#

# create variables that contain the item names

all.items <- names(select(data, contains("DEP", ignore.case = FALSE), contains("ANX", ignore.case = FALSE), contains("STR", ignore.case = FALSE)))
dep.items <- names(select(data, contains("DEP", ignore.case = FALSE)))
anx.items <- names(select(data, contains("ANX", ignore.case = FALSE)))
str.items <- names(select(data, contains("STR", ignore.case = FALSE)))

# create data frames that have data for each scale

all.data <- select(data, all_of(all.items))
dep.data <- select(data, all_of(dep.items))
anx.data <- select(data, all_of(anx.items))
str.data <- select(data, all_of(str.items))

# create data frames that have data for each scale but only for the evaluation data

all.eval <- select(filter(data, Sample == "evaluation"), all_of(all.items))
dep.eval <- select(filter(data, Sample == "evaluation"), all_of(dep.items))
anx.eval <- select(filter(data, Sample == "evaluation"), all_of(anx.items))
str.eval <- select(filter(data, Sample == "evaluation"), all_of(str.items))

# create data frames that have data for each scale but only for the validation data

all.vali <- select(filter(data, Sample == "validation"), all_of(all.items))
dep.vali <- select(filter(data, Sample == "validation"), all_of(dep.items))
anx.vali <- select(filter(data, Sample == "validation"), all_of(anx.items))
str.vali <- select(filter(data, Sample == "validation"), all_of(str.items))

#---                                        ---#
################# ICC Examples #################
#---                                        ---#

## ICC by difficulty ##

# create an ICC plot that shows different item characteristic curves based on difficulty

items.by.difficulty <- cbind(c(1,1,1,1,1),c(-2,-1,0,1,2))

twopl <- function(a, b, theta){
  1 / (1 + exp(-a * (theta - b)))}
theta <- seq(-4,4,.1)

items <- list()
for (i in 1:nrow(items.by.difficulty)) {
  items[[i]] <- twopl(a = items.by.difficulty[i,1], b = items.by.difficulty[i,2], theta = theta)
}

items <- data.frame(theta,as.data.frame(items))
colnames(items) <- c("theta","-2 theta","-1 theta","0 theta","1 theta","2 theta")
longer.format <- gather(items,item,measurement,2:6)

ggplot(longer.format, aes(theta, measurement, colour=item)) + 
  geom_line() + theme_bw() +
  xlab(expression(theta)) + 
  ylab(expression(P(theta))) + 
  theme(text = element_text(size=16),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.title=element_blank())

rm(items.by.difficulty,items,i,theta,longer.format,twopl)

## ICC by Discrimination ##

# create an ICC plot that shows different item characteristic curves based on discrimination

items.by.discrimination <- cbind(c(.1, .5, 1, 5),c(0, 0, 0, 0))

twopl <- function(a, b, theta){
  1 / (1 + exp(-a * (theta - b)))}
theta <- seq(-4,4,.1)

items <- list()
for (i in 1:nrow(items.by.discrimination)) {
  items[[i]] <- twopl(a = items.by.discrimination[i,1], b = items.by.discrimination[i,2], theta = theta)
}

items <- data.frame(theta,as.data.frame(items))
colnames(items) <- c("theta",".1",".5 ","1 ","5")
longer.format <- gather(items,item,measurement,2:5)

ggplot(longer.format, aes(theta, measurement, colour=item)) + 
  geom_line() + theme_bw() +
  xlab(expression(theta)) + 
  ylab(expression(P(theta))) +
  theme(text = element_text(size=16),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.title=element_blank())

rm(items.by.discrimination,items,i,theta,longer.format,twopl)

## ICC with Difficulty Line ##

# create an ICC plot that shows a 50% probability - i.e. the item difficulty

items.by.difficulty <- cbind(c(1),c(-1))

twopl <- function(a, b, theta){
  1 / (1 + exp(-a * (theta - b)))}
theta <- seq(-4,4,.1)

items <- list()
for (i in 1:nrow(items.by.difficulty)) {
  items[[i]] <- twopl(a = items.by.difficulty[i,1], b = items.by.difficulty[i,2], theta = theta)
}

items <- data.frame(theta,as.data.frame(items))
colnames(items) <- c("theta","item1")
longer.format <- gather(items,item,measurement,item1)

ggplot(longer.format, aes(theta, measurement, colour=item)) + 
  geom_line() + theme_bw() +
  xlab(expression(theta)) + 
  ylab(expression(P(theta))) + 
  geom_vline(aes(xintercept = items.by.difficulty[1,2])) +
  geom_hline(aes(yintercept = 0.5)) + theme_bw() + 
  theme(text = element_text(size=16),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.title=element_blank())

rm(items.by.difficulty,items,i,theta,longer.format,twopl)

#---                                          ---#
#################  Demographics #################
#---                                          ---#

# Plot for Gender
table(data$Gender)
table(data[ind,"Gender"])
table(data[!ind,"Gender"])

ggplot(data) +
  geom_bar(aes(Gender, fill = Sample), position = position_dodge(width = 0.7)) +
  geom_text(
    stat = "count",aes(Gender, fill = Sample, label = ..count..),
    vjust = -0.5,colour = "black",position = position_dodge(width = .7),size = 3) +
  theme(
    axis.text.x = element_text(colour = "grey20",size = 12,angle = 90,hjust = 1,vjust = .5,face = "plain"),
    axis.text.y = element_text(colour = "grey20",size = 12,angle = 0,hjust = 1,vjust = .5,face = "plain"),
    axis.title.x = element_text(colour = "grey20",size = 12,angle = 0,hjust = .5,vjust = -1,face = "plain"),
    axis.title.y = element_text(colour = "grey20",size = 12,angle = 90,hjust = .5,vjust = .5,face = "plain")) +
  ylab("Count") +
  scale_y_continuous(breaks = seq(0,505,50), limits = c(0,505))
ggsave("plots/Figure 1 - Gender.png")

# Tables and Plot for Age

table(data$AgeGroup)
table(data[ind,"AgeGroup"])
table(data[!ind,"AgeGroup"])
psych::describe(data$Age) %>%
  kable(digits=2, format="pandoc", caption=paste0("Age for Whole Sample (n = ", nrow(data),")"))
psych::describe(data[ind,"Age"]) %>%
  kable(digits=2, format="pandoc", caption=paste0("Age for Evaluation Sample (n = ", sample.split * nrow(data),")"))
psych::describe(data[!ind,"Age"]) %>%
  kable(digits=2, format="pandoc", caption=paste0("Age for Validation Sample (n = ", round((1 - sample.split) * nrow(data)),")"))

ggplot(data) +
  geom_bar(mapping = aes(AgeGroup, fill = Sample),
           position = position_dodge(width = 0.7)) +
  geom_text(
    stat = "count",aes(AgeGroup, fill = Sample, label = ..count..),
    vjust = -0.5,colour = "black",position = position_dodge(width = .7),size = 3) +
  xlab("Age Group") +
  scale_x_discrete(
    labels = c("[18,30)" = "18-29","[30,40)" = "30-39","[40,50)" = "40-49",
               "[50,60)" = "50-59","[60,70)" = "60-69","[70,80)" = "70-79","[80,90)" = "80+")) +
  theme(
    axis.text.x = element_text(colour = "grey20",size = 12,angle = 90,hjust = 1,vjust = .5,face = "plain"),
    axis.text.y = element_text(colour = "grey20",size = 12,angle = 0,hjust = 1,vjust = .5,face = "plain"),
    axis.title.x = element_text(colour = "grey20",size = 12,angle = 0,hjust = .5,vjust = 0,face = "plain"),
    axis.title.y = element_text(colour = "grey20",size = 12,angle = 90,hjust = .5,vjust = .5,face = "plain")) +
  ylab("Count") +
  scale_y_continuous(breaks = seq(0, 250, by = 25), limits = c(0, 250))
ggsave("plots/Figure 2 - Age.png")

# Employment Plot
table(data$Employment)
table(data[ind,"Employment"])
table(data[!ind,"Employment"])

ggplot(data) +
  geom_bar(mapping = aes(Employment, fill = Sample),position = position_dodge(width = 0.7)) +
  geom_text(stat = "count",aes(Employment, fill = Sample, label = ..count..),
            vjust = -0.5,colour = "black",position = position_dodge(width = .7),size = 3) +
  xlab("Employment Category") +
  scale_x_discrete(
    labels = c(
      "Full Time Student" = "Student",
      "Not Currently In Workforce" = "Unemployed",
      "Part Time or Full Time Employment" = "Employed",
      "Retired" = "Retired")) +
  theme(
    axis.text.x = element_text(colour = "grey20",size = 12,angle = 90,hjust = 1,vjust = .5,face = "plain"),
    axis.text.y = element_text(colour = "grey20",size = 12,angle = 0,hjust = 1,vjust = .5,face = "plain"),
    axis.title.x = element_text(colour = "grey20",size = 12,angle = 0,hjust = .5,vjust = 0,face = "plain"),
    axis.title.y = element_text(colour = "grey20",size = 12,angle = 90,hjust = .5,vjust = .5,face = "plain")) +
  ylab("Count") +
  scale_y_continuous(breaks = seq(0, 500, by = 50), limits = c(0, 500))
ggsave("plots/Figure 3 - Employment.png")

# Qualification Plot
table(data$Qualification)
table(data[ind,"Qualification"])
table(data[!ind,"Qualification"])

ggplot(data) +
  geom_bar(mapping =  aes(Qualification, fill = Sample),position = position_dodge(width = 0.7)) +
  geom_text(stat = "count",aes(Qualification, fill = Sample, label = ..count..),
            vjust = -0.5,colour = "black",position = position_dodge(width = .7),size = 3) +
  xlab("Qualification Category") +
  scale_x_discrete(
    labels = c(
      "Bachelor Degree" = "Bachelors",
      "Did not complete high school" = "Secondary Incomplete",
      "Completed High School" = "Secondary School",
      "Honours/Postgraduate Diploma" = "Honours/PG Diploma",
      "Masters/PhD" = "Postgraduate",
      "TAFE Certificate/Diploma/Trade" = "TAFE")) +
  theme(
    axis.text.x = element_text(colour = "grey20",size = 12,angle = 90,hjust = 1,vjust = .5,face = "plain"),
    axis.text.y = element_text(colour = "grey20",size = 12,angle = 0,hjust = 1,vjust = .5,face = "plain"),
    axis.title.x = element_text(colour = "grey20",size = 12,angle = 0,hjust = .5,vjust = 0,face = "plain"),
    axis.title.y = element_text(colour = "grey20",size = 12,angle = 90,hjust = .5,vjust = .5,face = "plain")) +
  ylab("Count") +
  scale_y_continuous(breaks = seq(0, 200, by = 25), limits = c(0, 200))
ggsave("plots/Figure 4 - Qualification.png")

# Scale Scores

psych::describe(filter(data, Sample == "evaluation")$dep)
psych::describe(filter(data, Sample == "evaluation")$anx)
psych::describe(filter(data, Sample == "evaluation")$str)

select(filter(data, Sample == "evaluation"), dep, anx, str) %>%
  as_tibble() %>%
  gather(key = "variable", value = "value") %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot(aes(color = variable)) +
  scale_x_discrete(name = "Scale") +
  scale_y_continuous(name = "Scale Score", limits = c(0,21), breaks = seq(0,21,7)) +
  theme(legend.position = "none")
ggsave("plots/Figure 5 - Scale Scores.png")

select(filter(data, Sample == "evaluation"), AgeGroup, dep, anx, str) %>%
  as_tibble() %>%
  gather(key = "variable", value = "value", -AgeGroup) %>%
  ggplot(aes(x = AgeGroup, y = value)) +
  geom_boxplot(aes(color = variable)) +
  scale_x_discrete(name = "Age") +
  scale_y_continuous(name = "Scale Score", limits = c(0,21), breaks = seq(0,21,7)) +
  scale_color_discrete(name = "Scale") +
  theme(legend.position = "bottom")
ggsave("plots/Figure 6 - Scale Scores by Age.png")

  
#---                                               ---#
################# Sample Independence #################
#---                                               ---#

data$Sample <- as.factor(data$Sample)

leveneTest(tot ~ Sample, data = data)
t.test(formula = tot ~ Sample, data = data) # t test to test differences in total score between samples

leveneTest(dep ~ Sample, data = data)
t.test(formula = dep ~ Sample, data = data) # t test to test differences in depression score between samples

leveneTest(anx ~ Sample, data = data)
t.test(formula = anx ~ Sample, data = data) # t test to test differences in anxiety score between samples

leveneTest(str ~ Sample, data = data)
t.test(formula = str ~ Sample, data = data) # t test to test differences in stress score between samples

leveneTest(Age ~ Sample, data = data)
t.test(formula = Age ~ Sample, data = data) # t test to test differences in age between samples

#---                                        ---#
################# Reliability #################
#---                                        ---#

psych::alpha(all.eval) # calculate cronbachs alpha for total scale
psych::alpha(dep.eval) # calculate cronbachs alpha for depression scale
psych::alpha(anx.eval) # calculate cronbachs alpha for anxiety scale
psych::alpha(str.eval) # calculate cronbachs alpha for stress scale

#---                                                ---#
################# CFA Model Comparisons #################
#---                                                ---#

# One factor model - Model A
model.a <- paste("tot =~ ", paste(all.items, collapse = "+"), sep = "")

# Test one factor model - Model A1
model.a1.fit <- cfa(model.a, data = all.eval, std.lv = TRUE, ordered = c(all.items))
fitMeasures(model.a1.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Three factor model - Model B
model.b <- paste(
  paste("dep =~ ", paste(dep.items, collapse = "+"), sep = ""),
  paste("anx =~ ", paste(anx.items, collapse = "+"), sep = ""),
  paste("str =~ ", paste(str.items, collapse = "+"), sep = ""),
  sep = "\n "
)

# Test three factor with covariances - Model B1
model.b1.fit <- cfa(model.b, data = all.eval, std.lv = TRUE, ordered=c(all.items))  
fitMeasures(model.b1.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")) 

# Test three factor with orthogonal latent variables - Model B2
model.b2.fit <- cfa(model.b, data = all.eval, orthogonal = TRUE, std.lv = TRUE, ordered=c(all.items))  
fitMeasures(model.b2.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")) 

# Two factor model - Model C
model.c <- paste(
  paste("dep =~ ", paste(dep.items, collapse = "+"), sep = ""),
  paste("dis =~ ", paste(anx.items, collapse = "+"), "+", paste(str.items, collapse = "+"), sep = ""),
  sep = "\n "
)

# Test two factor with covariances - Model C1
model.c1.fit <- cfa(model.c, data = all.eval, std.lv = TRUE, ordered=c(all.items))
fitMeasures(model.c1.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Test two factor with orthogonal latent variables - Model C2
model.c2.fit <- cfa(model.c, data = all.eval, orthogonal = TRUE, std.lv = TRUE, ordered=c(all.items))
fitMeasures(model.c2.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Tripartite model with orthongonal negative affect only - Model D2
model.d1 <- paste(
  # two domain specific factors
  paste("dep =~ ", paste(dep.items, collapse = "+"), sep = ""),
  paste("anx =~ ", paste(anx.items, collapse = "+"), sep = ""),
  # one negative affect factor
  paste("na =~ ", paste(all.items, collapse = "+"), sep = ""),
  # orthogonal factors
  paste("na ~~ 0*dep", sep=""),
  paste("na ~~ 0*anx", sep=""),
  # seperate all by lines
  sep = "\n "
)

# Tripartite model with all orthogonal factors - Model D2
model.d2 <- paste(
  # two domain specific factors
  paste("dep =~ ", paste(dep.items, collapse = "+"), sep = ""),
  paste("anx =~ ", paste(anx.items, collapse = "+"), sep = ""),
  # one negative affect factor
  paste("na =~ ", paste(all.items, collapse = "+"), sep = ""),
  # orthogonal factors
  paste("na ~~ 0*dep", sep=""),
  paste("na ~~ 0*anx", sep=""),
  paste("dep ~~ 0*anx", sep=""),
  # seperate all by lines
  sep = "\n "
)

# Test the Tripartite model with orthogonal negative affect - Model D1
model.d1.fit <- cfa(model.d1, data = all.eval, std.lv = TRUE, ordered = c(all.items))
fitMeasures(model.d1.fit , c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Test the Tripartite model with all orthogonal factors = Model D2
model.d2.fit <- cfa(model.d2, data = all.eval, std.lv = TRUE, ordered = c(all.items))
fitMeasures(model.d2.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Quadripartite model with orthongonal negative affect - Model E1
model.e1 <- paste(
  # three domain specific factors
  paste("dep =~ ",paste(dep.items, collapse = "+"), sep = ""),
  paste("anx =~ ",paste(anx.items, collapse = "+"), sep = ""),
  paste("str =~ ",paste(str.items, collapse = "+"), sep = ""),
  # one negative affect factor
  paste("na =~ ",paste(all.items, collapse = "+"), sep = ""),
  # orthogonal factors
  paste("na ~~ 0*dep", sep=""),
  paste("na ~~ 0*anx", sep=""),
  paste("na ~~ 0*str", sep=""),
  # variances to 1
  paste("na ~~ 1*na", sep=""),
  paste("dep ~~ 1*dep", sep =""),
  paste("anx ~~ 1*anx", sep=""),
  paste("str ~~ 1*str", sep=""),
  # add new lines
  sep = "\n "
)

# Quadripartite model with all orthogonal factors - Model E2
model.e2 <- paste(
  # three domain specific factors
  paste("dep =~ ",paste(dep.items, collapse = "+"), sep = ""),
  paste("anx =~ ",paste(anx.items, collapse = "+"), sep = ""),
  paste("str =~ ",paste(str.items, collapse = "+"), sep = ""),
  # one negative affect factor
  paste("na =~ ",paste(all.items, collapse = "+"), sep = ""),
  # orthogonal factors
  paste("na ~~ 0*dep", sep=""),
  paste("na ~~ 0*anx", sep=""),
  paste("na ~~ 0*str", sep=""),
  paste("dep ~~ 0*anx", sep=""),
  paste("dep ~~ 0*str", sep=""),
  paste("anx ~~ 0*str", sep=""),
  # variances to 1
  paste("na ~~ 1*na", sep=""),
  paste("dep ~~ 1*dep", sep=""),
  paste("anx ~~ 1*anx", sep=""),
  paste("str ~~ 1*str", sep=""),
  # add new lines
  sep = "\n "
)

# Test the quadripartite model with orthogonal negative affect
model.e1.fit <- cfa(model.e1, data = all.eval, std.lv = TRUE, ordered = c(all.items))
fitMeasures(model.e1.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Test the quadripartite model with all orthogonal factors
model.e2.fit <- cfa(model.e2, data = all.eval, std.lv = TRUE, ordered = c(all.items))
fitMeasures(model.e2.fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Compare the results of the quadripartite models
anova(model.e1.fit, model.e2.fit)

# Compare the model fit for the eval sample and the vali sample
model.e2.fit.vali <- cfa(model.e2, data = all.vali, std.lv = TRUE, ordered = c(all.items))
fitMeasures(model.e2.fit.vali, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Parameter Estimates (this model used because it is more parsimonious)
parameterEstimates(model.e2.fit, standardized = TRUE) %>%
  filter(op == "=~") %>%
  dplyr::select("Trait" = lhs, op, "Item" = rhs, "SE" = se, "Z" = z, "p-value" = pvalue, "Std. Beta" = std.all) %>%
  kable(digits = 3, format="pandoc", caption="Factor Loadings for DASS-21")

rm(model.a, model.b, model.c, model.d1, model.d2, model.e1)
rm(model.a1.fit, model.b1.fit, model.b2.fit, model.c1.fit, model.c2.fit, model.d1.fit, model.d2.fit, model.e1.fit)

#---                                              ---#
################# Mokken Total Scale #################
#---                                             ---#

all.mokken <- coefH(all.eval) # conduct mokken analysis for total scale
all.mokken$Hi # check the Loevinger's H values (>.3 is ok, >.5 is excellent)
all.mokken$H # check the Loevinger's H values (>.3 is ok, >.5 is excellent)
all.aisp <- aisp(all.eval) # conduct the automatic item selection for the total scale
all.aisp  # print results of automatic item selection

#---                                ---#
################# MIRT #################
#---                                ---#

bfactor(data = all.eval, model = c(rep(1,7),rep(2,7),rep(3,7))) # fails even with 2000 EM cycles, need more data

#---                                           ---#
################# DIF Preparation #################
#---                                          ---#

# load data for dif
dif.data <- filter(data, Sample == "evaluation")

# 
dif.data$Gender <- dif.data$Gender

# Collapse employment categories
dif.data$Employment <- recode(dif.data$Employment,
                     "Student" = "Unemployed",
                     "Employed" = "Employed",
                     "Unemployed" = "Unemployed",
                     "Retired" = "Unemployed")

#
dif.data$Qualification <- recode(dif.data$Qualification,
                        "Secondary School Incomplete" = "Non Tertiary",
                        "TAFE" = "Tertiary",
                        "Secondary" = "Non Tertiary",
                        "Undergraduate" = "Tertiary", 
                        "Honours" = "Tertiary",
                        "Masters/PhD" = "Tertiary")

dif.data$AgeGroup <- recode(dif.data$AgeGroup,
                        "18 to 29" = "18 to 29",
                        "30 to 39" = "30 to 39",
                        "40 to 49" = "40 to 49",
                        "50 to 59" = "50 to 59",
                        "60 to 69" = "60+",
                        "70 to 79" = "60+",
                        "80 to 90" = "60+")


#---                                                     ---#
################# Unidimensional Depression #################
#---                                                     ---#

## Depression  Item Descriptives ##

psych::describe(dep.eval) %>%
  kable(digits = 2, format = "pandoc", caption = paste0("Depression for Evaluation Sample (n = ", round(sample.split * nrow(data)),")"))
  
## Depression CFA ##

dep.model <- paste("dep =~ ",paste(dep.items, collapse = "+"), sep = "")
dep.fit <- cfa(dep.model, data = dep.eval, std.lv = TRUE, ordered = c(dep.items))
fitMeasures(dep.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) 

# Parameter Estimates
parameterEstimates(dep.fit, standardized = TRUE) %>%
  filter(op == "=~") %>%
  dplyr::select(Item = rhs, SE = se, Z = z, "p-value" = pvalue, "Std. Beta" = std.all) %>%
  kable(digits = 3, format="pandoc", caption="Factor Loadings for Depression")

## Depression Mokken ##

dep.mokken <- coefH(dep.eval)
dep.mokken$Hi #check the Loevinger's H values (>.3 is ok, >.5 is excellent)
dep.mokken$H #check the Loevinger's H values (>.3 is ok, >.5 is excellent)
dep.aisp <- aisp(dep.eval)
dep.aisp  

## Depression Item Response Theory ##

dep.grm <- mirt(dep.eval, model = 1, type = "graded")
itemfit(dep.grm) # check item fit, should be > .05, or at minimum .01
M2(dep.grm, type = "C2") 
marginal_rxx(dep.grm)
dep.pars <- coef(dep.grm, simplify = TRUE)
dep.pars

plot(dep.grm, type = "trace", xlim = c(-6, 6), main = "")
plot(dep.grm, type = "info", xlim = c(-6, 6), main = "Depression Test Information")
plot(dep.grm, type = "rxx", xlim = c(-6, 6), main = "")

## Depression Local Independence ##

residuals(dep.grm, type = "Q3", digits = 2) # positive correlations over .2 that are problematic 

## Depression Differential Item Functioning ##

for (i in 1:7) { print(table(dif.data$Gender, dep.eval[,i]))}

dep.mg <- multipleGroup(data = dep.eval, model = 1, group = dif.data$Gender, invariance = "free_mean")
DIF(MGmodel = dep.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = dep.mg)

for (i in 1:7) { print(table(dif.data$Employment, dep.eval[,i]))}

dep.mg <- multipleGroup(data = dep.eval, model = 1, group = dif.data$Employment, invariance = "free_mean")
DIF(MGmodel = dep.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = dep.mg)

for (i in 1:7) { print(table(dif.data$Qualification, dep.eval[,i]))}

dep.mg <- multipleGroup(data = dep.eval, model = 1, group = dif.data$Qualification, invariance = "free_mean")
DIF(MGmodel = dep.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = dep.mg)

## Depression CAT Simulation (Simulated Participants) ##

# CAT with simulated participants
dep.cat <- list()
for (i in 1:length(sem)) { dep.cat[[i]] <- cat_simulation(model = dep.grm, sem = sem[i]) }
dep.cat <- bind_rows(dep.cat)
dep.cat$reliability <- reliability
dep.cat

# Performance of simulated participants by group
dep.Theta <- matrix(rnorm(n = 10000, mean = 0, sd = 1)) # Simulate theta for 10000 people
dep.pattern <- generate_pattern(mo = dep.grm, Theta = dep.Theta) ### Generate response pattern
dep.cat.design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
dep.cat.group <- mirtCAT(mo = dep.grm, local_pattern = dep.pattern, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = dep.cat.design, cl = cl) # run
stopCluster(cl) # stop the cluster
dep.est.Theta <- laply(dep.cat.group, function(y) y$thetas) # list of theta
dep.sim.groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(dep.Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(dep.Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(dep.Theta > Theta_groups[i-1] & dep.Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(dep.cat.group[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(dep.Theta[rows,1], dep.est.Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(dep.Theta[rows,1] - dep.est.Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((dep.Theta[rows,1] - dep.est.Theta[rows])^2))
  nParticipants <- length(rows)
  dep.sim.groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                  "averageItems" = round(average_items, 3), "correlation" = round(correlation, 3), 
                                  "bias" = round(bias, 3), "rmsd" = round(rmsd, 3))
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
dep.sim.groups <- bind_rows(dep.sim.groups)
dep.sim.groups

## Depression CAT Simulation (Real Participants) ##

# CAT with real participants
dep.cat.real <- list()
for (i in 1:length(sem)) { dep.cat.real[[i]] <- cat_simulation(model = dep.grm, sem = sem[i], pattern =  dep.vali) }
dep.cat.real <- bind_rows(dep.cat.real)
dep.cat.real$reliability <- reliability
dep.cat.real

# Performance of real participants by group
dep.Theta.real <- matrix(fscores(object = dep.grm, method = "EAP", response.pattern = dep.vali, append_response.pattern = FALSE)[,1])
dep.pattern.real <- data.matrix(dep.vali)
attr(dep.pattern.real, "Theta") <- matrix(fscores(object = dep.grm, method = "EAP", response.pattern = dep.vali, append_response.pattern = FALSE)[,1])
dep.cat.design.real <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
dep.cat.group.real <- mirtCAT(mo = dep.grm, local_pattern = dep.pattern.real, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = dep.cat.design.real, cl = cl) # run
stopCluster(cl) # stop the cluster
dep.est.Theta.real <- laply(dep.cat.group.real, function(y) y$thetas) # list of theta
dep.sim.groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(dep.Theta.real < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(dep.Theta.real > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(dep.Theta.real > Theta_groups[i-1] & dep.Theta.real < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(dep.cat.group.real[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(dep.Theta.real[rows,1], dep.est.Theta.real[rows]) ### correlation between true theta and estimated theta
  bias <- mean(dep.Theta.real[rows,1] - dep.est.Theta.real[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((dep.Theta.real[rows,1] - dep.est.Theta.real[rows])^2))
  nParticipants <- length(rows)
  dep.sim.groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                   "averageItems" = average_items, "correlation" = correlation, 
                                   "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
dep.sim.groups <- bind_rows(dep.sim.groups)
dep.sim.groups

#---                                                  ---#
################# Unidimensional Anxiety #################
#---                                                  ---#

## Anxiety  Item Descriptives ##

psych::describe(anx.eval) %>%
  kable(digits = 2, format = "pandoc", caption = paste0("Anxiety for Evaluation Sample (n = ", round(sample.split * nrow(data)),")"))

## Anxiety CFA ##

anx.model <- paste("anx =~ ",paste(anx.items, collapse = "+"), sep = "")
anx.fit <- cfa(anx.model, data = anx.eval, std.lv = TRUE, ordered = c(anx.items))
fitMeasures(anx.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) 

# Parameter Estimates
parameterEstimates(anx.fit, standardized = TRUE) %>%
  filter(op == "=~") %>%
  dplyr::select(Item = rhs, SE = se, Z = z, "p-value" = pvalue, "Std. Beta" = std.all) %>%
  kable(digits = 3, format = "pandoc", caption = "Factor Loadings for Anxiety")

## Anxiety Mokken ##

anx.mokken <- coefH(anx.eval)
anx.mokken$Hi #check the Loevinger's H values (>.3 is ok, >.5 is excellent)
anx.mokken$H #check the Loevinger's H values (>.3 is ok, >.5 is excellent)
anx.aisp <- aisp(anx.eval)
anx.aisp  

## Anxiety Item Response Theory ##

anx.grm <- mirt(anx.eval, model = 1, type = "graded")
itemfit(anx.grm) # check item fit, should be > .05, or at minimum .01
M2(anx.grm, type = "C2") 
marginal_rxx(anx.grm)

plot(anx.grm, type = "trace", xlim = c(-6, 6), main = "")
plot(anx.grm, type = "info", xlim = c(-6, 6), main = "Anxiety Test Information")
plot(anx.grm, type = "rxx", xlim = c(-6, 6), main = "")

## Anxiety Local Independence ##

residuals(anx.grm, type = "Q3", digits = 2) # positive correlations over .2 that are problematic 

## Anxiety Differential Item Functioning ##

for (i in 1:7) { print(table(dif.data$Gender, anx.eval[,i]))}

anx.mg <- multipleGroup(data = anx.eval, model = 1, group = dif.data$Gender, invariance = "free_mean")
DIF(MGmodel = anx.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = anx.mg)

for (i in 1:7) { print(table(dif.data$Employment, anx.eval[,i]))}

anx.mg <- multipleGroup(data = anx.eval, model = 1, group = dif.data$Employment, invariance = "free_mean")
DIF(MGmodel = anx.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = anx.mg)

for (i in 1:7) { print(table(dif.data$Qualification, anx.eval[,i]))}

anx.mg <- multipleGroup(data = anx.eval, model = 1, group = dif.data$Qualification, invariance = "free_mean")
DIF(MGmodel = anx.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = anx.mg)

## Anxiety CAT Simulation (Simulated Participants) ##

# CAT with simulated participants
anx.cat <- list()
for (i in 1:length(sem)) { anx.cat[[i]] <- cat_simulation(model = anx.grm, sem = sem[i]) }
anx.cat <- bind_rows(anx.cat)
anx.cat$reliability <- reliability
anx.cat

# Performance of simulated participants by group
anx.Theta <- matrix(rnorm(n = 10000, mean = 0, sd = 1)) # Simulate theta for 10000 people
anx.pattern <- generate_pattern(mo = anx.grm, Theta = anx.Theta) ### Generate response pattern
anx.cat.design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
anx.cat.group <- mirtCAT(mo = anx.grm, local_pattern = anx.pattern, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = anx.cat.design, cl = cl) # run
stopCluster(cl) # stop the cluster
anx.est.Theta <- laply(anx.cat.group, function(y) y$thetas) # list of theta
anx.sim.groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(anx.Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(anx.Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(anx.Theta > Theta_groups[i-1] & anx.Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(anx.cat.group[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(anx.Theta[rows,1], anx.est.Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(anx.Theta[rows,1] - anx.est.Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((anx.Theta[rows,1] - anx.est.Theta[rows])^2))
  nParticipants <- length(rows)
  anx.sim.groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                   "averageItems" = round(average_items, 3), "correlation" = round(correlation, 3), 
                                   "bias" = round(bias, 3), "rmsd" = round(rmsd, 3))
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
anx.sim.groups <- bind_rows(anx.sim.groups)
anx.sim.groups

## Anxiety CAT Simulation (Real Participants) ##

# CAT with real participants
anx.cat.real <- list()
for (i in 1:length(sem)) { anx.cat.real[[i]] <- cat_simulation(model = anx.grm, sem = sem[i], pattern =  anx.vali) }
anx.cat.real <- bind_rows(anx.cat.real)
anx.cat.real$reliability <- reliability
anx.cat.real

# Performance of real participants by group
anx.Theta.real <- matrix(fscores(object = anx.grm, method = "EAP", response.pattern = anx.vali, append_response.pattern = FALSE)[,1])
anx.pattern.real <- data.matrix(anx.vali)
attr(anx.pattern.real, "Theta") <- matrix(fscores(object = anx.grm, method = "EAP", response.pattern = anx.vali, append_response.pattern = FALSE)[,1])
anx.cat.design.real <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
anx.cat.group.real <- mirtCAT(mo = anx.grm, local_pattern = anx.pattern.real, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = anx.cat.design.real, cl = cl) # run
stopCluster(cl) # stop the cluster
anx.est.Theta.real <- laply(anx.cat.group.real, function(y) y$thetas) # list of theta
anx.sim.groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(anx.Theta.real < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(anx.Theta.real > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(anx.Theta.real > Theta_groups[i-1] & anx.Theta.real < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(anx.cat.group.real[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(anx.Theta.real[rows,1], anx.est.Theta.real[rows]) ### correlation between true theta and estimated theta
  bias <- mean(anx.Theta.real[rows,1] - anx.est.Theta.real[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((anx.Theta.real[rows,1] - anx.est.Theta.real[rows])^2))
  nParticipants <- length(rows)
  anx.sim.groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                   "averageItems" = average_items, "correlation" = correlation, 
                                   "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
anx.sim.groups <- bind_rows(anx.sim.groups)
anx.sim.groups

#---                                                 ---#
################# Unidimensional Stress #################
#---                                                ---#

## Stress  Item Descriptives ##

psych::describe(str.eval) %>%
  kable(digits = 2, format = "pandoc", caption = paste0("Stress for Evaluation Sample (n = ", round(sample.split * nrow(data)),")"))

## Stress CFA ##

str.model <- paste("str =~ ",paste(str.items, collapse = "+"), sep = "")
str.fit <- cfa(str.model, data = str.eval, std.lv = TRUE, ordered = c(str.items))
fitMeasures(str.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) 

# Parameter Estimates
parameterEstimates(str.fit, standardized = TRUE) %>%
  filter(op == "=~") %>%
  dplyr::select(Item = rhs, SE = se, Z = z, "p-value" = pvalue, "Std. Beta" = std.all) %>%
  kable(digits = 3, format = "pandoc", caption = "Factor Loadings for Stress")

## Stress Mokken ##

str.mokken <- coefH(str.eval)
str.mokken$Hi #check the Loevinger's H values (>.3 is ok, >.5 is excellent)
str.mokken$H #check the Loevinger's H values (>.3 is ok, >.5 is excellent)
str.aisp <- aisp(str.eval)
str.aisp  

## Stress Item Response Theory ##

str.grm <- mirt(str.eval, model = 1, type = "graded")
itemfit(str.grm) # check item fit, should be > .05, or at minimum .01
M2(str.grm, type = "C2") 
marginal_rxx(str.grm)
str.pars <- coef(str.grm, simplify = TRUE)
str.pars

plot(str.grm, type = "trace", xlim = c(-6, 6), main = "")
plot(str.grm, type = "info", xlim = c(-6, 6), main = "Stress Test Information")
plot(str.grm, type = "rxx", xlim = c(-6, 6), main = "")

## Stress Local Independence ##

residuals(str.grm, type = "Q3", digits = 2) # positive correlations over .2 that are problematic 

## Stress Differential Item Functioning ##

for (i in 1:7) { print(table(dif.data$Gender, str.eval[,i]))}

str.mg <- multipleGroup(data = str.eval, model = 1, group = dif.data$Gender, invariance = "free_mean")
DIF(MGmodel = str.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = str.mg)

for (i in 1:7) { print(table(dif.data$Employment, str.eval[,i]))}

str.mg <- multipleGroup(data = str.eval, model = 1, group = dif.data$Employment, invariance = "free_mean")
DIF(MGmodel = str.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = str.mg)

for (i in 1:7) { print(table(dif.data$Qualification, str.eval[,i]))}

str.mg <- multipleGroup(data = str.eval, model = 1, group = dif.data$Qualification, invariance = "free_mean")
DIF(MGmodel = str.mg, which.par = c("a1", "d1", "d2", "d3"), scheme = "add", plotdif = TRUE)
DTF(mod = str.mg)

## Stress CAT Simulation (Simulated Participants) ##

# CAT with simulated participants
str.cat <- list()
for (i in 1:length(sem)) { str.cat[[i]] <- cat_simulation(model = str.grm, sem = sem[i]) }
str.cat <- bind_rows(str.cat)
str.cat$reliability <- reliability
str.cat

# Performance of simulated participants by group
str.Theta <- matrix(rnorm(n = 10000, mean = 0, sd = 1)) # Simulate theta for 10000 people
str.pattern <- generate_pattern(mo = str.grm, Theta = str.Theta) ### Generate response pattern
str.cat.design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
str.cat.group <- mirtCAT(mo = str.grm, local_pattern = str.pattern, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = str.cat.design, cl = cl) # run
stopCluster(cl) # stop the cluster
str.est.Theta <- laply(str.cat.group, function(y) y$thetas) # list of theta
str.sim.groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(str.Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(str.Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(str.Theta > Theta_groups[i-1] & str.Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(str.cat.group[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(str.Theta[rows,1], str.est.Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(str.Theta[rows,1] - str.est.Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((str.Theta[rows,1] - str.est.Theta[rows])^2))
  nParticipants <- length(rows)
  str.sim.groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                   "averageItems" = round(average_items, 3), "correlation" = round(correlation, 3), 
                                   "bias" = round(bias, 3), "rmsd" = round(rmsd, 3))
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
str.sim.groups <- bind_rows(str.sim.groups)
str.sim.groups

## Stress CAT Simulation (Real Participants) ##

# CAT with real participants
str.cat.real <- list()
for (i in 1:length(sem)) { str.cat.real[[i]] <- cat_simulation(model = str.grm, sem = sem[i], pattern =  str.vali) }
str.cat.real <- bind_rows(str.cat.real)
str.cat.real$reliability <- reliability
str.cat.real

# Performance of real participants by group
str.Theta.real <- matrix(fscores(object = str.grm, method = "EAP", response.pattern = str.vali, append_response.pattern = FALSE)[,1])
str.pattern.real <- data.matrix(str.vali)
attr(str.pattern.real, "Theta") <- matrix(fscores(object = str.grm, method = "EAP", response.pattern = str.vali, append_response.pattern = FALSE)[,1])
str.cat.design.real <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
str.cat.group.real <- mirtCAT(mo = str.grm, local_pattern = str.pattern.real, start_item = "MPWI", method = "EAP", criteria = "MPWI", design = str.cat.design.real, cl = cl) # run
stopCluster(cl) # stop the cluster
str.est.Theta.real <- laply(str.cat.group.real, function(y) y$thetas) # list of theta
str.sim.groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(str.Theta.real < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(str.Theta.real > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(str.Theta.real > Theta_groups[i-1] & str.Theta.real < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(str.cat.group.real[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(str.Theta.real[rows,1], str.est.Theta.real[rows]) ### correlation between true theta and estimated theta
  bias <- mean(str.Theta.real[rows,1] - str.est.Theta.real[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((str.Theta.real[rows,1] - str.est.Theta.real[rows])^2))
  nParticipants <- length(rows)
  str.sim.groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                   "averageItems" = average_items, "correlation" = correlation, 
                                   "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
str.sim.groups <- bind_rows(str.sim.groups)
str.sim.groups