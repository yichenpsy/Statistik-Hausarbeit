#########HA7 Aufgabe1
library(openintro)
library(ggplot2)

## Daten laden
data(births, package="openintro") 
head(births)
help(births)

## b)
# Data auswählen
df <- births[c("weeks","sex_baby","smoke")]

# NA-Einträgen entfernen
df <- na.omit(df)

# 2x2 balanciertes Design
table(df$sex_baby, df$smoke)
Bedin_female_non<- df[sample (which( (df$sex_baby == "female") & (df$smoke == "nonsmoker") ), 19),]
Bedin_male_non<- df[sample (which( (df$sex_baby == "male") & (df$smoke == "nonsmoker") ), 19),]
Bedin_female_smo<- df[sample (which( (df$sex_baby == "female") & (df$smoke == "smoker") ), 19),]
Bedin_male_smo<- df[sample (which( (df$sex_baby == "male") & (df$smoke == "smoker") ), 19),]

df2<-rbind(Bedin_female_non, Bedin_male_non, Bedin_female_smo, Bedin_male_smo)
df2

## c)
range(df2$weeks)
interaction.plot(df2$sex_baby,df2$smoke, df2$weeks,
                 ylab = "Schwangerschaftswoche der Entbindung",
                 xlab = "Geschlecht des Babys",
                 main = "Schwangerschaftswoche ~ Geschlecht * Rauchen",
                 col=c("red", "blue"), 
                 trace.label="Rauchen", 
                 fixed=TRUE,
                 type="b",
                 lty=1, lwd = 1.2,
                 pch=19,
                 ylim = c(37,40),
                 cex.main=1.25, cex.lab=1.25, cex.axis=1.2)
                              

interaction.plot(df2$smoke, df2$sex_baby, df2$weeks,
                 ylab = "Schwangerschaftswoche der Entbindung",
                 xlab = " Rauchen der Mutter",
                 main = "Schwangerschaftswoche ~ Rauchen * Geschlecht",
                 col=c("red", "blue"), 
                 trace.label="Geschlecht", 
                 fixed=TRUE,
                 type="b",
                 lty=1, lwd = 1.2,
                 pch=19,
                 ylim = c(37,40),
                 cex.main=1.25, cex.lab=1.25, cex.axis=1.2)


## e)
contrasts(df2$sex_baby) <- contr.sum(2)
contrasts(df2$smoke) <- contr.sum(2)

reg <- lm(weeks ~ sex_baby*smoke, data=df2)
summary(reg)

##g)
anova(reg)
