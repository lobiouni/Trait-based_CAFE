###############################################################
## Manuscript: Zooplankton community trait mean value rather 
## than diversity determines ecosystem functioning 
##
## Date: 08-July-2024
## Updated: 10-July-2024
##
## Author: Lorena Pinheiro-Silva
##
## Encoding: UTF-8
###############################################################

# To make reproducible outputs
set.seed(1234) 

# To avoid scientific notaion
options(scipen=999)

# Color pallet
cols=c("palegreen", "limegreen","forestgreen","darkgreen") 

##=============================================
## Installing and loading packages 
##=============================================

if(!require(data.table)){install.packages("data.table");library(data.table)}
if(!require(tidyverse)){install.packages("tidyverse");library(tidyverse)}
if(!require(magrittr)){install.packages("magrittr");library(dataMaid)}
if(!require(vegan)){install.packages("vegan");library(vegan)}
if(!require(picante)){install.packages("picante");library(picante)}
if(!require(FD)){install.packages("FD");library(FD)}
if(!require(dataMaid)){install.packages("dataMaid");library(dataMaid)}
if(!require(conover.test)){install.packages("conover.test");library(conover.test)}
if(!require(Rmisc)){install.packages("Rmisc");library(Rmisc)}
if(!require(lmPerm)){install.packages("lmPerm");library(lmPerm)}
if(!require(ggthemes)){install.packages("ggthemes");library(ggthemes)}
if(!require(patchwork)){install.packages("patchwork");library(patchwork)}
if(!require(rcompanion)){install.packages("rcompanion");library(rcompanion)}
if(!require(ggpubr)){install.packages("ggpubr");library(ggpubr)}


##=============================================
## Importing data 
##=============================================

data<-read.csv("https://raw.githubusercontent.com/lobiouni/Trait-based_CAFE/main/Trait-based_CAFE_data.csv")

# Checking dataset structure
str(data)

# Fix any nominal column mistakenly read in as numeric
data<-data%>%
  dplyr::mutate_at(c("P", "PL", "Repl"), as.factor)

# Fix any numeric column mistakenly read in as nominal
data<-data%>%
  dplyr::mutate_at(c("S"), as.numeric)

# Rechecking dataset structure
str(data)

##=============================================
## Exploratory data Analysis
##=============================================

# use MakeDataReport with HTML as output
makeDataReport(data, output = "html", replace = TRUE)

# Statistic descriptors for RUE
data %>%
  dplyr::group_by(P, PL) %>%
  dplyr::summarize(mean = mean(RUEzp_D65, na.rm = TRUE),
                   sd = sd(RUEzp_D65, na.rm = TRUE),
                   min = min(RUEzp_D65, na.rm = TRUE),
                   mx = min(RUEzp_D65, na.rm = TRUE))

# Statistic descriptors for SD
data %>%
  dplyr::group_by(P, PL) %>%
  dplyr::summarize(mean = mean(SDshn_dens, na.rm = TRUE),
                   sd = sd(SDshn_dens, na.rm = TRUE),
                   min = min(SDshn_dens, na.rm = TRUE),
                   mx = min(SDshn_dens, na.rm = TRUE))

# Statistic descriptors for CAS
data %>%
  dplyr::group_by(P, PL) %>%
  dplyr::summarize(mean = mean(CAS_dens, na.rm = TRUE),
                   sd = sd(CAS_dens, na.rm = TRUE),
                   min = min(CAS_dens, na.rm = TRUE),
                   mx = min(CAS_dens, na.rm = TRUE))

# Statistic descriptors for S
data %>%
  dplyr::group_by(P, PL) %>%
  dplyr::summarize(mean = mean(S, na.rm = TRUE),
                   sd = sd(S, na.rm = TRUE),
                   min = min(S, na.rm = TRUE),
                   mx = min(S, na.rm = TRUE))

##=============================================
## ANOVAs using permutation tests
##=============================================
## CAS ----

# ANOVA using permutation tests
aovp1 = aovp(data$CAS_dens ~ data$P * data$PL, seqs = T)
summary.aovp(aovp1)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$CAS_dens, data$P, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$CAS_dens, data$PL, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$CAS_dens, data$interaction, method = "bonferroni")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted
)

# Converting "comparisons" to a factor
conover_result_df<- conover_result_df%>%
  dplyr::mutate_at(c("comparisons"), as.factor)

# Generate the CLD using cldList
cld_result <- cldList(P_adjusted ~ comparisons, 
                      data = conover_result_df,
                      threshold = 0.05,
                      remove.zero = FALSE)

# Reordering the letters to match with the data.frame
str(cld_result)

cld_result<- cld_result %>% 
  dplyr::mutate_at(c("Group"), as.factor)%>%
  dplyr::mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  dplyr::arrange(Group)

# Preparing the data to make the plot
df1 <- data%>%summarySE(measurevar="CAS_dens", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df1$Comp <- cld_result$Letter

# Cleaning the data.frame
df1<-df1%>%
  dplyr::mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  dplyr::mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  dplyr::mutate_at(c("P", "PL", "Comp"), as.factor)%>%
  dplyr::mutate_at(c("CAS_dens", "sd", "se", "ci"), as.numeric)

# Figure 1a
A = ggplot(df1, aes(colour=P, y=CAS_dens, x=P)) + 
  geom_point(position=position_dodge(0.9), size=3) + 
  geom_errorbar(aes(ymin=CAS_dens-se, ymax=CAS_dens+se), width=.2,
                position=position_dodge(0.9)) +
  geom_jitter(data = data, aes(y=CAS_dens), alpha=0.3, position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  theme_few() +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="black", size = 14, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 14, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +
  geom_text(aes(label=Comp, y = (CAS_dens + se) + 0.08), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = "",
       y = "CAS")+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))

## SD ----

# ANOVA using permutation tests
aovp2 = aovp(data$SDshn_dens ~ data$P * data$PL)
summary.aovp(aovp2)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$SDshn_dens, data$P, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$SDshn_dens, data$PL, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$SDshn_dens, data$interaction, method = "bonferroni")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted
)

# Converting "comparisons" to a factor
conover_result_df<- conover_result_df%>%
  dplyr::mutate_at(c("comparisons"), as.factor)

# Generate the CLD using cldList
cld_result <- cldList(P_adjusted ~ comparisons, 
                      data = conover_result_df,
                      threshold = 0.05,
                      remove.zero = FALSE)

# Reordering the letters to match with the data.frame
str(cld_result)

cld_result<- cld_result %>% 
  dplyr::mutate_at(c("Group"), as.factor)%>%
  dplyr::mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  dplyr::arrange(Group)

# Preparing the data to make the plot
df2 <- data%>%summarySE(measurevar="SDshn_dens", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df2$Comp <- cld_result$Letter

# Cleaning the data.frame
df2<-df2%>%
  dplyr::mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  dplyr::mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  dplyr::mutate_at(c("P", "PL", "Comp"), as.factor)%>%
  dplyr::mutate_at(c("SDshn_dens", "sd", "se", "ci"), as.numeric)

# Figure 1b
B = ggplot(df2, aes(colour=P, y=SDshn_dens, x=P)) + 
  geom_point(position=position_dodge(0.9), size=3) + 
  geom_errorbar(aes(ymin=SDshn_dens-se, ymax=SDshn_dens+se), width=.2,
                position=position_dodge(0.9)) +
  geom_jitter(data = data, aes(y=SDshn_dens), alpha=0.3, position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  theme_few() +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="black", size = 14, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 14, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +
  geom_text(aes(label=Comp, y = (SDshn_dens + se) + 0.15), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = "",
       y = "SD")+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))


## RICHNESS ----

# ANOVA using permutation tests
aovp3 = aovp(data$S ~ data$P * data$PL)
summary.aovp(aovp3)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$S, data$P, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$S, data$PL, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$S, data$interaction, method = "bonferroni")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted
)

# Converting "comparisons" to a factor
conover_result_df<- conover_result_df%>%
  dplyr::mutate_at(c("comparisons"), as.factor)

# Generate the CLD using cldList
cld_result <- cldList(P_adjusted ~ comparisons, 
                      data = conover_result_df,
                      threshold = 0.05,
                      remove.zero = FALSE)

# Reordering the letters to match with the data.frame
str(cld_result)

cld_result<- cld_result %>% 
  dplyr::mutate_at(c("Group"), as.factor)%>%
  dplyr::mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  dplyr::arrange(Group)

# Preparing the data to make the plot
df3 <- data%>%summarySE(measurevar="S", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df3$Comp <- cld_result$Letter

# Cleaning the data.frame
df3<-df3%>%
  dplyr::mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  dplyr::mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  dplyr::mutate_at(c("P", "PL", "Comp"), as.factor)%>%
  dplyr::mutate_at(c("S", "sd", "se", "ci"), as.numeric)

# Figure 1c
C = ggplot(df3, aes(colour=P, y=S, x=P)) + 
  geom_point(position=position_dodge(0.9), size=3) + 
  geom_errorbar(aes(ymin=S-se, ymax=S+se), width=.2,
                position=position_dodge(0.9)) +
  geom_jitter(data = data, aes(y=S), alpha=0.3, position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  theme_few() +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="black", size = 14, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 14, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = -1.5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +
  geom_text(aes(label=Comp, y = (S + se) + 0.6), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = "P-addition levels",
       y = "S")+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))

## RUE ----

# ANOVA using permutation tests
aovp4 = aovp(log(data$RUEzp_D65) ~ data$P * data$PL)
summary.aovp(aovp4)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)
str(data)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$RUEzp_D65, data$P, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$RUEzp_D65, data$PL, method = "bonferroni")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$RUEzp_D65, data$interaction, method = "bonferroni")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted
)

# Converting "comparisons" to a factor
conover_result_df<- conover_result_df%>%
  dplyr::mutate_at(c("comparisons"), as.factor)

# Generate the CLD using cldList
cld_result <- cldList(P_adjusted ~ comparisons, 
                      data = conover_result_df,
                      threshold = 0.05,
                      remove.zero = FALSE)

# Reordering the letters to match with the data.frame
str(cld_result)

cld_result<- cld_result %>% 
  dplyr::mutate_at(c("Group"), as.factor)%>%
  dplyr::mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  dplyr::arrange(Group)

# Preparing the data to make the plot
df4 <- data%>%summarySE(measurevar="RUEzp_D65", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df4$Comp <- cld_result$Letter

# Cleaning the data.frame
df4<-df4%>%
  dplyr::mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  dplyr::mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  dplyr::mutate_at(c("P", "PL", "Comp"), as.factor)%>%
  dplyr::mutate_at(c("RUEzp_D65", "sd", "se", "ci"), as.numeric)

# Figure 1d
D = ggplot(df4, aes(colour=P, y=log(RUEzp_D65), x=P)) + 
  geom_point(position=position_dodge(0.9), size=3) + 
  geom_errorbar(aes(ymin=log(RUEzp_D65-se), ymax=log(RUEzp_D65+se)), width=.2,
                position=position_dodge(0.9)) +
  geom_jitter(data = data, aes(y=log(RUEzp_D65)), alpha=0.3, position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  theme_few() +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="black", size = 14, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 14, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = -1.5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = 0.5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +  
  geom_text(aes(label=Comp, y = log(RUEzp_D65 + se)+0.4), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = "P-addition levels",
       y = expression(paste("ln RUE"[ZP])))+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))

## Figure 1----
Fig1<- (A | B) / (C | D)

Fig1<-Fig1+plot_annotation(tag_levels = c("a", "b", "c", "d"), tag_suffix = ')')& 
  theme(plot.tag = element_text(face="bold"))
Fig1

tiff(filename= "Fig1.tif",
     height=5800,
     width=7600,
     units="px",
     res=800,
     compression="lzw")
plot(Fig1)
dev.off()

png(filename= "Fig1.png",
    height=5800,
    width=7600,
    units="px",
    res=800)
plot(Fig1)
dev.off()

##=============================================
## variation partitioning
##=============================================
#### mod1: RUE ~ Community structure (SD, CAS, and S), P-addition and HH (Figure 2a) ####

mod1 = varpart(log(data$RUEzp_D65), 
               cbind(data$SDshn_dens, data$CAS_dens, data$S), 
               data$P,
               data$PL)
plot(mod1, 
     Xnames = c("CSS&D", "P", "HH"), 
     cex = 1.1, 
     id.size = 0.8, 
     adj=0.5,
     digits = 1,
     cutoff =-Inf)

# Testing the model
rda.mod1 <- rda(log(data$RUEzp_D65) ~ 
                  cbind(data$SDshn_dens, data$CAS_dens, data$S)+ 
                  data$P+
                  data$PL)
anova(rda.mod1, step=200, perm.max=999) 

# Testing for the pure effect of Biodiversity
rda.mod2 <- rda(log(data$RUEzp_D65) ~ 
                  cbind(data$SDshn_dens, data$CAS_dens, data$S)+ 
                  Condition(data$P)+
                  Condition(data$PL))
anova(rda.mod2, step=200, perm.max=999) 

# Testing for the pure effect of P-addition
rda.mod3 <- rda(log(data$RUEzp_D65) ~ 
                  Condition(cbind(data$SDshn_dens, data$CAS_dens, data$S))+ 
                  data$P+
                  Condition(data$PL))
anova(rda.mod3, step=200, perm.max=999) 

# Testing for the pure effect of HH
rda.mod4 <- rda(log(data$RUEzp_D65) ~ 
                  Condition(cbind(data$SDshn_dens, data$CAS_dens, data$S))+ 
                  Condition(data$P)+
                  data$PL)
anova(rda.mod4, step=200, perm.max=999) 

#### mod2: RUE ~ SD, CAS and S (Figure 2b) ####

mod2 = varpart(log(data$RUEzp_D65), 
               data$SDshn_dens,
               data$CAS_dens,
               data$S)
plot(mod2, 
     Xnames = c("SD", "CAS", "S"), 
     cex = 1.3, 
     id.size = 0.8, 
     digits = 1,
     cutoff =-Inf)

# Testing the model
rda.mod5 <- rda(log(data$RUEzp_D65) ~ 
                  data$CAS_dens+ 
                  data$SDshn_dens+ 
                  data$S)
anova(rda.mod5, step=200, perm.max=999) 

# Testing for the pure effect of SD
rda.mod7 <- rda(log(data$RUEzp_D65) ~ 
                  data$SDshn_dens +
                  Condition(data$CAS_dens) + 
                  Condition(data$S))
anova(rda.mod7, step=200, perm.max=999)

# Testing for the pure effect of CAS
rda.mod6 <- rda(log(data$RUEzp_D65) ~ 
                  Condition(data$SDshn_dens) + 
                  data$CAS_dens+
                  Condition(data$S))
anova(rda.mod6, step=200, perm.max=999) 

# Testing for the pure effect of S
rda.mod8 <- rda(log(data$RUEzp_D65) ~ 
                  Condition(data$SDshn_dens) + 
                  Condition(data$CAS_dens) + 
                  data$S)
anova(rda.mod8, step=200, perm.max=999)


#### mod3: RUE ~ SD, CAS and P-addition (Figure 2c) ####

mod3 = varpart(log(data$RUEzp_D65), 
               data$SDshn_dens, 
               data$CAS_dens,
               data$P)

par(mfrow=c(1,1))
plot(mod3, 
     Xnames = c("SD", "CAS", "P"), 
     cex = 1.3, 
     id.size = 1, 
     digits = 1,
     cutoff =-Inf)

# Testing the model
rda.mod9 <- rda(log(data$RUEzp_D65) ~ 
                  data$SDshn_dens+ 
                  data$CAS_dens+
                  data$P)
anova(rda.mod9, step=200, perm.max=999) 

# Testing for the pure effect of SD
rda.mod10 <- rda(log(data$RUEzp_D65) ~ 
                   data$SDshn_dens+ 
                   Condition(data$CAS_dens)+
                   Condition(data$P))
anova(rda.mod10, step=200, perm.max=999) 

# Testing for the pure effect of CAS
rda.mod11 <- rda(log(data$RUEzp_D65) ~ 
                   Condition(data$SDshn_dens)+ 
                   data$CAS_dens+
                   Condition(data$P))
anova(rda.mod11, step=200, perm.max=999) 

# Testing for the pure effect of P-addition
rda.mod12 <- rda(log(data$RUEzp_D65) ~ 
                   Condition(data$SDshn_dens)+ 
                   Condition(data$CAS_dens)+
                   data$P)
anova(rda.mod12, step=200, perm.max=999) 






#### mod4: RUE ~ SD, CAS and habitat heterogeneity (Figure 2d) ####

mod4 =varpart(log(data$RUEzp_D65), 
              data$SDshn_dens, 
              data$CAS_dens,
              data$PL)

par(mfrow=c(1,1))
plot(mod4, 
     Xnames = c("SD", "CAS", "HH"), 
     cex = 1.3, 
     id.size = 0.8, 
     digits = 1,
     cutoff =-Inf)

# Testing the model
rda.mod13 <- rda(log(data$RUEzp_D65) ~ 
                   data$SDshn_dens+ 
                   data$CAS_dens+
                   data$PL)
anova(rda.mod13, step=200, perm.max=999) 

# Testing for the pure effect of SD
rda.mod14 <- rda(log(data$RUEzp_D65) ~ 
                   data$SDshn_dens + 
                   Condition(data$CAS_dens) + 
                   Condition(data$PL))
anova(rda.mod14, step=200, perm.max=999) 

# Testing for the pure effect of CAS
rda.mod15 <- rda(log(data$RUEzp_D65) ~ 
                   Condition(data$SDshn_dens) + 
                   data$CAS_dens + 
                   Condition(data$PL))
anova(rda.mod15, step=200, perm.max=999) 

# Testing for the pure effect of habitat heterogeneity (HH)
rda.mod16 <- rda(log(data$RUEzp_D65) ~ 
                   Condition(data$SDshn_dens) + 
                   Condition(data$CAS_dens) + 
                   data$PL)
anova(rda.mod16, step=200, perm.max=999) 

##=============================================
## Correlation plot (Figure S1)
##=============================================
A = ggscatter(
  data, x = "SDshn_dens", y = "CAS_dens",
  color = "P",
  alpha = 0.7,
  size=3,
  rug=T) +
  scale_x_continuous(limits=c(-1,4.5), expand=c(0,0))+
  scale_y_continuous(limits=c(0,3.5), expand=c(0,0))+
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  geom_smooth(method = "lm", se = T, color = "#CDC9A5", fill = "#EEEED1", size = 1, alpha = 0.4, linetype = 1) +
  stat_cor(label.x=0.1, method = "spearman") +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(colour="black", size = 13, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 13, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank(), # Background underneath legend keys
        legend.background = element_blank()) +
  labs(title = NULL, y= "CAS (mm)", x = "SD")

B = ggscatter(
  data, x = "S", y = "CAS_dens",
  color = "P",
  alpha = 0.7,
  size=3,
  rug=T) +
  scale_y_continuous(limits=c(0,3.5), expand=c(0,0))+
  scale_x_continuous(limits=c(4,16.5), expand=c(0,0))+
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  geom_smooth(method = "lm", se = T, color = "#CDC9A5", fill = "#EEEED1", size = 1, alpha = 0.4, linetype = 1) +
  stat_cor(label.x=4.3, method = "spearman") +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(colour="black", size = 13, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 13, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank(), # Background underneath legend keys
        legend.background = element_blank()) +
  labs(title = NULL, x= "S", y = "CAS (mm)")

C = ggscatter(
  data, x = "S", y = "SDshn_dens",
  color = "P",
  alpha = 0.7,
  size=3,
  rug=T) +
  scale_y_continuous(limits=c(-1,4.5), expand=c(0,0))+
  scale_x_continuous(limits=c(4,16.5), expand=c(0,0))+
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  geom_smooth(method = "lm", se = T, color = "#CDC9A5", fill = "#EEEED1", size = 1, alpha = 0.4, linetype = 1) +
  stat_cor(label.x=4.3, method = "spearman") +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(colour="black", size = 13, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 13, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank(), # Background underneath legend keys
        legend.background = element_blank()) +
  labs(title = NULL, x= "S", y = "SD")

FigS1<- A / B / C

FigS1<-FigS1+plot_annotation(tag_levels = c("a", "b", "c"), tag_suffix = ')')& 
  theme(plot.tag = element_text(face="bold"))
FigS1

tiff(filename= "FigS1.tif",
     height=7000,
     width=4400,
     units="px",
     res=800,
     compression="lzw")
print(FigS1)
dev.off()

png(filename= "FigS1.png",
    height=7000,
    width=4400,
    units="px",
    res=800)
plot(FigS1)
dev.off()

