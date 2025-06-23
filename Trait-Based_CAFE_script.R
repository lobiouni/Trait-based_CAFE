###############################################################
## Manuscript: Quantifying the contribution of community trait 
## mean and diversity to ecosystem functioning 
##
## Date: 08-July-2024
## Last updated: 23-June-2025
##
## Author: Lorena Pinheiro-Silva
##
## Encoding: UTF-8
###############################################################

# To make reproducible outputs
set.seed(1234) 

# To avoid scientific notation
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

## Figure 1----

Fig1<-ggplot(body_size, aes(x = BS)) +
  geom_density(fill = "#f9eae1", color = "#cc8b86") +
  labs(title = NULL, x = "Body Size (mm)", y = "Density") +
  #theme_minimal(base_size = 12)+
  #theme_few(base_size = 12) +
  theme_bw()+
  theme(text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour="black", linetype = "solid", linewidth=0.5),
        panel.background = element_rect(fill = "white"),
        plot.margin=unit(c(0.3, 0.3, 0.3, 0.3),"cm"),
        axis.text.x = element_text(colour="black", size = 9, angle = 0, vjus = 1),
        axis.text.y = element_text(colour="black", size = 9, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = - 1.5, size = 12, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 12, margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  scale_x_continuous(
    limits = c(0.2, 2.0),  
    breaks = seq(0.2, 2.0, by = 0.2))

Fig1


jpeg(filename = "Fig1.jpeg",
     width = 2.95,  
     height = 2.36, 
     units = "in",
     res = 600,    
     quality = 100)
plot(Fig1)  
dev.off()
##=============================================
## ANOVAs using permutation tests
##=============================================
## CAS ----

# ANOVA using permutation tests
aovp1 = aovp(data$CAS_dens ~ data$P * data$PL, seqs = T)
summary.aovp(aovp1)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons - P factor
conover_result<-conover.test(data$CAS_dens, data$P, method = "BH")

# Pairwise comparisons - PL factor
conover_result<-conover.test(data$CAS_dens, data$PL, method = "BH")

# Pairwise comparisons - Interaction
conover_result<-conover.test(data$CAS_dens, data$interaction, method = "BH")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted)

# Round relevant columns
conover_result_df$P_adjusted <- round(conover_result_df$P_adjusted, 3)

# Converting "comparisons" to a factor
conover_result_df<- conover_result_df%>%
  dplyr::mutate_at(c("comparisons"), as.factor)

# Extract significant comparisons (p < 0.05)
sig_comparisons <- conover_result_df[conover_result_df$P_adjusted < 0.05, ]

# Generate the CLD using cldList
cld_result <- cldList(P_adjusted ~ comparisons, 
                      data = conover_result_df,
                      threshold = 0.05,
                      remove.zero = FALSE)

# Reordering the letters to match with the data.frame
str(cld_result)

cld_result<- cld_result %>% 
  mutate_at(c("Group"), as.factor)%>%
  mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  arrange(Group)

# Preparing the data to make the plot
df1 <- data%>%summarySE(measurevar="CAS_dens", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df1$Comp <- cld_result$Letter

# Cleaning the data.frame
df1<-df1%>%
  mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  mutate_at(c("P", "PL"), as.factor)

# Figure 2a
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
        axis.text.x = element_text(colour="black", size = 12, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 12, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +
  geom_text(aes(label=Comp, y = (CAS_dens + se) + 0.03), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = "",
       y = expression(paste("CWM"[BS], "(mm)")))+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))
A

## SD ----

# ANOVA using permutation tests
aovp2 = aovp(data$SDshn_dens ~ data$P * data$PL)
summary.aovp(aovp2)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$SDshn_dens, data$P, method = "BH")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$SDshn_dens, data$PL, method = "BH")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$SDshn_dens, data$interaction, method = "BH")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted)

# Round relevant columns
conover_result_df$P_adjusted <- round(conover_result_df$P_adjusted, 3)

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
  mutate_at(c("Group"), as.factor)%>%
  mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  arrange(Group)

# Preparing the data to make the plot
df2 <- data%>%summarySE(measurevar="SDshn_dens", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df2$Comp <- cld_result$Letter

# Cleaning the data.frame
df2<-df2%>%
  mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  mutate_at(c("P", "PL"), as.factor)

# Figure 2b
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
        axis.text.x = element_text(colour="black", size = 12, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 12, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +
  geom_text(aes(label=Comp, y = (SDshn_dens + se) + 0.2), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = "",
       y = "SD")+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))
B

## RICHNESS ----

# ANOVA using permutation tests
aovp3 = aovp(data$S ~ data$P * data$PL)
summary.aovp(aovp3)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$S, data$P, method = "BH")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$S, data$PL, method = "BH")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$S, data$interaction, method = "BH")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted)

# Round relevant columns
conover_result_df$P_adjusted <- round(conover_result_df$P_adjusted, 3)

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
  mutate_at(c("Group"), as.factor)%>%
  mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  arrange(Group)

# Preparing the data to make the plot
df3 <- data%>%summarySE(measurevar="S", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df3$Comp <- cld_result$Letter

# Cleaning the data.frame
df3<-df3%>%
  mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  mutate_at(c("P", "PL"), as.factor)

# Figure 2c
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
        axis.text.x = element_text(colour="black", size = 12, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 12, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = -1.5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +
  geom_text(aes(label=Comp, y = (S + se) + 0.8), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = expression("P-addition levels (µg/L)"), 
       y = "S")+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))
C

## RUE ----

# ANOVA using permutation tests
aovp4 = aovp(log(data$RUEzp_D65) ~ data$P * data$PL)
summary.aovp(aovp4)

# Prepare data for pairwise comparison
data$interaction <- interaction(data$P, data$PL)

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$RUEzp_D65, data$P, method = "BH")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$RUEzp_D65, data$PL, method = "BH")

# Pairwise comparisons (Conover-Iman post-hoc test)
conover_result<-conover.test(data$RUEzp_D65, data$interaction, method = "BH")

# Extract pairwise adjusted p-values and comparisons
p_adjusted <- conover_result$P.adjusted
comparisons <- conover_result$comparisons

# Create a data frame for easier manipulation
conover_result_df <- data.frame(
  comparisons = comparisons,
  P_adjusted = p_adjusted)

# Round relevant columns
conover_result_df$P_adjusted <- round(conover_result_df$P_adjusted, 3)

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
  mutate_at(c("Group"), as.factor)%>%
  mutate(Group = fct_relevel(Group, c("0.APL", "10.APL", "100.APL", "1000.APL", "0.NPL", "10.NPL", "100.NPL", "1000.NPL")))%>%
  arrange(Group)

# Preparing the data to make the plot
df4 <- data%>%summarySE(measurevar="RUEzp_D65", groupvars="interaction", na.rm = T)

# Adding the letters to the data.frame
df4$Comp <- cld_result$Letter

# Cleaning the data.frame
df4<-df4%>%
  mutate(P=c("0", "10", "100", "1000", "0", "10", "100", "1000"))%>%
  mutate(PL=c("APL", "APL", "APL", "APL", "NPL", "NPL", "NPL", "NPL"))%>%
  mutate_at(c("P", "PL"), as.factor)

# Figure 2d
D = ggplot(df4, aes(colour=P, y=log(RUEzp_D65), x=P)) + 
  geom_point(position=position_dodge(0.9), size=3) + 
  geom_errorbar(aes(ymin=log(RUEzp_D65-se), ymax=log(RUEzp_D65+se)), width=.2,
                position=position_dodge(0.9)) +
  geom_jitter(data = RUE, aes(y=log(RUEzp_D65)), alpha=0.3, position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  theme_few() +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        panel.border = element_blank(),
        axis.text.x = element_text(colour="black", size = 12, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 12, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = -1.5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = 0.5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none") +  
  geom_text(aes(label=Comp, y = log(RUEzp_D65 + se)+0.4), size = 4.5, color = "Gray25", show.legend = FALSE, position = position_dodge(0.9))+
  labs(x = expression("P-addition levels (µg/L)"),
       y = expression(paste("ln RUE"[ZP])))+
  guides(colour = guide_legend(expression(paste("P-addition levels"))))
D

## Figure 2----
Fig2<- (A | B) / (C | D)

Fig2<-Fig2+plot_annotation(tag_levels = c("a", "b", "c", "d"), tag_suffix = ')')& 
  theme(plot.tag = element_text(face="bold"))
Fig2

jpeg(filename = "Fig2.jpeg",
     width = 8,     
     height = 6.5,   
     units = "in",
     res = 600,        
     quality = 100)    
print(Fig2)
dev.off()

##=============================================
## Model info
##=============================================
# Extract model info----
# Function to extract SS and p-values
extract_ss_pvalue <- function(aovp_object) {
  summary_df <- as.data.frame(summary(aovp_object)[[1]])  # Extract the summary and convert to data frame
  summary_df$`Pr(Prob)` <- round(summary_df$`Pr(Prob)`, 3)  # Round p-values to 3 decimal places
  summary_df$`R Sum Sq` <- round(summary_df$`R Sum Sq`, 3)  # Round SS values to 3 decimal places
  return(summary_df[, c("R Sum Sq", "Pr(Prob)")])  # Return SS and p-value columns
}

# Extract SS and p-values from each summary
summary_aovp1 <- extract_ss_pvalue(aovp1)
summary_aovp2 <- extract_ss_pvalue(aovp2)
summary_aovp3 <- extract_ss_pvalue(aovp3)
summary_aovp4 <- extract_ss_pvalue(aovp4)

# Calculate the pseudo-F for ANOVA with permutation----
## CAS: P term----
# summary
summary_ovp1 <- summary(aovp1)
str(summary_ovp1)

# Extract necessary values
ss_terms <- summary_ovp1[[1]]$`R Sum Sq`
df_terms <- summary_ovp1[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[1]  
ssw <- ss_terms[4]  s

# Degrees of freedom
dfb <- df_terms[1]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## CAS: PL term----
# summary 
summary_ovp1 <- summary(aovp1)
str(summary_ovp1)

# Extract necessary values
ss_terms <- summary_ovp1[[1]]$`R Sum Sq`
df_terms <- summary_ovp1[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[2]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[2]  
dfw <- df_terms[4]
# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## CAS: PxPL term----
# summary 
summary_ovp1 <- summary(aovp1)
str(summary_ovp1)

# Extract necessary values
ss_terms <- summary_ovp1[[1]]$`R Sum Sq`
df_terms <- summary_ovp1[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[3]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[3]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## SD: P term----
# summary
summary_ovp2 <- summary(aovp2)
str(summary_ovp2)

# Extract necessary values
ss_terms <- summary_ovp2[[1]]$`R Sum Sq`
df_terms <- summary_ovp2[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[1] 
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[1]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## SD: PL term----
#summary
summary_ovp2 <- summary(aovp2)
str(summary_ovp2)

# Extract necessary values
ss_terms <- summary_ovp2[[1]]$`R Sum Sq`
df_terms <- summary_ovp2[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[2]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[2]  
dfw <- df_terms[4] 

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw

# Display the pseudo-F statistic
round(pseudo_F, 3) 

## SD: PxPL term----
# Display the summary of the aovp model
summary_ovp2 <- summary(aovp2)
str(summary_ovp2)

# Extract necessary values
ss_terms <- summary_ovp2[[1]]$`R Sum Sq`
df_terms <- summary_ovp2[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[3]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[3]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## S: P term----
#summary
summary_ovp3 <- summary(aovp3)
str(summary_ovp3)

# Extract necessary values
ss_terms <- summary_ovp3[[1]]$`R Sum Sq`
df_terms <- summary_ovp3[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[1]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[1]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## S: PL term----
# summary 
summary_ovp3 <- summary(aovp3)
str(summary_ovp3)

# Extract necessary values
ss_terms <- summary_ovp3[[1]]$`R Sum Sq`
df_terms <- summary_ovp3[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[2]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[2]  
dfw <- df_terms[4] 

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## S: PxPL term----
# summary
summary_ovp3 <- summary(aovp3)
str(summary_ovp3)

# Extract necessary values
ss_terms <- summary_ovp3[[1]]$`R Sum Sq`
df_terms <- summary_ovp3[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[3]  
ssw <- ss_terms[4] 

# Degrees of freedom
dfb <- df_terms[3]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## RUE: P term----
# summary
summary_ovp4 <- summary(aovp4)
str(summary_ovp4)

# Extract necessary values
ss_terms <- summary_ovp4[[1]]$`R Sum Sq`
df_terms <- summary_ovp4[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[1]  
ssw <- ss_terms[4]  

# Degrees of freedom
dfb <- df_terms[1]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## RUE: PL term----
# summary
summary_ovp4 <- summary(aovp4)
str(summary_ovp4)

# Extract necessary values
ss_terms <- summary_ovp4[[1]]$`R Sum Sq`
df_terms <- summary_ovp4[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[2]  
ssw <- ss_terms[4] 

# Degrees of freedom
dfb <- df_terms[2]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

## RUE: PxPL term----
# summary
summary_ovp4 <- summary(aovp4)
str(summary_ovp4)

# Extract necessary values
ss_terms <- summary_ovp4[[1]]$`R Sum Sq`
df_terms <- summary_ovp4[[1]]$Df

# Identify the sums of squares for between-groups and within-groups
ssb <- ss_terms[3]  
ssw <- ss_terms[4] 

# Degrees of freedom
dfb <- df_terms[3]  
dfw <- df_terms[4]  

# Calculate mean squares
msb <- ssb / dfb
msw <- ssw / dfw

# Calculate the pseudo-F statistic
pseudo_F <- msb / msw
round(pseudo_F, 3) 

##=============================================
## variation partitioning
##=============================================
#### mod1: RUE ~ Community structure (SD, CAS, and S), P-addition and HH (Figure 3a) ####

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

#### mod2: RUE ~ SD, CAS and S (Figure 3b) ####

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


#### mod3: RUE ~ SD, CAS and P-addition (Figure 3c) ####

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






#### mod4: RUE ~ SD, CAS and habitat heterogeneity (Figure 3d) ####

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
cols=c("palegreen", "limegreen","forestgreen","darkgreen") 

# Create a custom data frame with r and p values for each facet (PL)
annotation_data <- data.frame(
  PL = c("APL", "NPL"),
  label = c(
    "italic(r)==-0.51*','~italic(p)==0.045",  
    "italic(r)==-0.22*','~italic(p)==0.42"),
  x = c(1.5, 1.5),
  y = c(1.5, 1.5))

A = ggscatter(
  RUE, x = "SDshn_dens", y = "CAS_dens",
  color = "P",
  alpha = 0.7,
  size = 3,
  rug = T) +
  scale_x_continuous(limits = c(-1, 4.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
  scale_colour_manual("P", values = cols, labels = c("0", "10", "100", "1000")) +
  facet_wrap(~PL) +
  geom_smooth(method = "lm", se = T, color = "#CDC9A5", fill = "#EEEED1", size = 1, alpha = 0.4, linetype = 1) +
  #stat_cor(label.x = 0.1, method = "spearman") +
  theme(
    strip.background = element_rect(colour = "black", fill = "white"),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(colour = "black", size = 14, vjust = .5, angle = 0),
    axis.text.y = element_text(colour = "black", size = 14, angle = 0),
    axis.title.x = element_text(angle = 0, vjust = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.y = element_text(angle = 90, vjust = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.line = element_line(size = 0.5, colour = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.key = element_blank(),
    legend.background = element_blank()
  ) +
  labs(title = NULL, y = expression(paste("CWM"[BS], "(mm)")), x = "SD") +
  guides(color=guide_legend(title=expression(paste("P-addition levels (µg/L)")), order = 1))+
  # Add custom text for each facet
  geom_text(data = annotation_data, aes(x = x, y = y, label = label, group = PL), 
            parse=TRUE, size = 5, inherit.aes = FALSE)

A

annotation_data <- data.frame(
  PL = c("APL", "NPL"),
  label = c(
    "italic(r)==-0.05*','~italic(p)==0.85",  
    "italic(r)==0.02*','~italic(p)==0.95"),
  x = c(9.3, 9.3),
  y = c(3, 3))

B = ggscatter(
  RUE, x = "S", y = "CAS_dens",
  color = "P",
  alpha = 0.7,
  size=3,
  #add = "reg.line",
  #shape="P",
  rug=T) +
  scale_y_continuous(limits=c(0, 3.5), expand=c(0,0))+
  scale_x_continuous(limits=c(4, 16.5), expand=c(0,0))+
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  geom_smooth(method = "lm", se = T, color = "#CDC9A5", fill = "#EEEED1", size = 1, alpha = 0.4, linetype = 1) +
  #stat_cor(label.x=4.3, method = "spearman") +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(colour="black", size = 14, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 14, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank(), # Background underneath legend keys
        legend.background = element_blank()) +
  labs(title = NULL, x= "S", y = expression(paste("CWM"[BS], "(mm)")),) +
  # Add custom text for each facet
  geom_text(data = annotation_data, aes(x = x, y = y, label = label, group = PL), 
            parse=TRUE, size = 5, inherit.aes = FALSE)
B

annotation_data <- data.frame(
  PL = c("APL", "NPL"),
  label = c(
    "italic(r)==0.01*','~italic(p)==0.96",  
    "italic(r)==0.28*','~italic(p)==0.29"),
  x = c(9, 9),
  y = c(4, 4))

C = ggscatter(
  RUE, x = "S", y = "SDshn_dens",
  color = "P",
  alpha = 0.7,
  size=3,
  #add = "reg.line",
  #shape="P",
  rug=T) +
  scale_y_continuous(limits=c(-1, 4.9), expand=c(0,0))+
  scale_x_continuous(limits=c(4,16.5), expand=c(0,0))+
  scale_colour_manual("P",  values = cols, labels = c("0", "10", "100", "1000"))+
  facet_wrap(~PL) +
  geom_smooth(method = "lm", se = T, color = "#CDC9A5", fill = "#EEEED1", size = 1, alpha = 0.4, linetype = 1) +
  #stat_cor(label.x=4.3, method = "spearman") +
  theme(strip.background = element_rect(colour="black", fill="white"),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(colour="black", size = 14, vjus = .5,angle = 0),
        axis.text.y = element_text(colour="black", size = 14, angle = 0),
        axis.title.x = element_text(angle = 0, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjus = .5, size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(size=0.5, colour = "black"),
        legend.position = "none",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank(), # Background underneath legend keys
        legend.background = element_blank()) +
  labs(title = NULL, x= "S", y = "SD")+
  # Add custom text for each facet
  geom_text(data = annotation_data, aes(x = x, y = y, label = label, group = PL), 
            parse=TRUE, size = 5, inherit.aes = FALSE)
C

FigS1<- A / B / C

FigS1<-FigS1+plot_annotation(tag_levels = c("a", "b", "c"), tag_suffix = ')')& 
  theme(plot.tag = element_text(face="bold"))
FigS1

jpeg(filename = "FigS1.jpeg",
     width = 5.5,  
     height = 8.75,
     units = "in",
     res = 600,     
     quality = 100) 
plot(FigS1)
dev.off()