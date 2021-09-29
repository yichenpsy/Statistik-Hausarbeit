#### HA6 Aufgabe 1####
library(dplyr)
# Daten vorbereiten
setwd ("/Users/yichenzhong/Desktop/HF_SE2_Statistik2/HA/HA6")
df <- read.table("reaction.dat.txt", header = TRUE)

##b) ALM
x1 <- c(ifelse(df$control=="C" & df$gender == "F", 1, 0))
x2 <- c(ifelse(df$control=="C" & df$gender == "M", 1, 0))
x3 <- c(ifelse(df$control=="T" & df$gender == "F", 1, 0))
x4 <- c(ifelse(df$control=="T" & df$gender == "M", 1, 0))
X <- cbind(x1, x2, x3, x4)

## c)
beta.dach <- solve(t(X) %*% X) %*% t(X) %*% df$time
c <- c(-1, -1, 1, 1)
(Psi.dach <- c%*%beta.dach) # 0.101892

## d)
df1 <- df
df1$gender <- ifelse(df1$gender == "F", "weiblich", ifelse( df1$gender == "M", "männlich", 0))
df1$control <- ifelse(df1$control == "C", "ohne Telefonat", ifelse( df1$control == "T", "mit Telefonat", 0))

SEM_fun <- function(x) {
  n <- length(x)
  SEM <- sd(x)/sqrt(n) 
  return(SEM)}
SEM <- tapply(df1$time, list(df1$gender,df1$control ), SEM_fun)
t(SEM)

Mittelwerte <- tapply(df1$time, list(df1$gender,df1$control ), mean)
t(Mittelwerte)

library(Hmisc)

interaction.plot(df1$control,df1$gender, df1$time,
                 ylab = "Reaktionszeit (Sekunden): Mittelwert +/- 1 SEM",
                 xlab = "Telefonieren",
                 main = "Mittelwerte der Reaktionszeit für die beiden Bedingungen des Telefonierens \n   im Vergleich der beiden Geschlechtergruppen",
                 col=c("blue","red"), 
                 trace.label="Geschlecht", 
                 fixed=TRUE,
                 type="b",
                 lty=1, lwd = 2,
                 pch=19,
                 ylim = c(1.3, 1.5),
                 cex.main=1, cex.lab=1.1, cex.axis=1)

errbar(x=rep(c(1,2),each= 2),
       Mittelwerte,
       Mittelwerte+SEM,
       Mittelwerte-SEM,
       add = TRUE,
       type = "n")

## e) Statistische Hypothesentest
# Damit wir die Funktion später nutzen können, definieren wir sie zunächst:
t.bruch <- function(y,X,c,Psi.0=0){            
  beta.dach      <-  solve(t(X) %*% X) %*% t(X) %*% y
  Psi.dach       <-  c%*%beta.dach
  n              <-  length(y) 
  p              <-  length(beta.dach)
  y.dach         <-  X %*% beta.dach 
  s2             <- (1/(n-p)) * sum((y-y.dach)^2)
  est.V.Psi.dach <-  s2 * t(c)%*% solve(t(X) %*% X)%*%c
  (Psi.dach - Psi.0) / (sqrt(est.V.Psi.dach))           
}

(t.emp <-t.bruch(df$time,X,c)) # 2.172129

# Entscheidung treffen
alpha <- 0.05
n <- nrow(df)
p <- length(beta.dach)
(t.krit <- qt(1-alpha,n-p)) # 1.67
abs(t.emp) > t.krit
# => H0 verwerfen

(p.value <- pt(t.emp,n-p, lower.tail = FALSE)) # 0.017
p.value < alpha
# => H0 verwerfen

##f)
mean.T <- mean(df[df$control=="T","time"]) #1.445571
mean.C <- mean(df[df$control=="C","time"]) #1.389613

ST <- sd(df$time)
d.dach <- (mean.T - mean.C)/ST #0.640132

## g)
n <- nrow(X)
p <- ncol(X)
deg <- n-p

ncp <- seq(0, 10, by=.1) # eine Reihe möglicher Nonzentralitätsparameter

alpha <- 0.05
t_krit <- qt(1-alpha,deg)

power <- 1-pt(t_krit,deg,ncp=ncp)
d.psi <- ncp/as.vector((t(c)%*% solve(t(X) %*% X)%*%c)^(-0.5))

tab <- data.frame(ncp,power,d.psi)

(ncp.1 <- tab$ncp[tab$power>=0.9][1]) # 3
(d.1 <- ncp.1/as.vector((t(c)%*% solve(t(X) %*% X)%*%c)^(-0.5))) #1.66856

# Plot
plot(d.psi,power,
     type='l',
     xlab='Effektstärke d.psi',ylab='Power',
     main='Power im t-Test der parametrischen Funktion, einseitig')
abline(h = tab$power[tab$power>0.90][1], lty = 2)
abline(v = tab$d.psi[tab$power>0.90][1], lty = 2)

## h)
power.apriori <- function(n1){ 
  d.psi <- 1
  x1 <- c(            rep(1,n1), rep(0,5*n1))
  x2 <- c(rep(0,n1),   rep(1,n1), rep(0,4*n1))
  x3 <- c(rep(0,2*n1), rep(1,n1), rep(0,3*n1))
  x4 <- c(rep(0,3*n1), rep(1,n1), rep(0,2*n1))
  x5 <- c(rep(0,4*n1), rep(1,n1), rep(0,n1))
  x6 <- c(rep(0,5*n1), rep(1,n1))
X <- cbind(x1,x2,x3,x4,x5,x6)
c <- c(1,-0.5,-0.5, 0, 0, 0) # Koeffizientenvektor
ncp <- as.vector((t(c)%*% solve(t(X) %*% X)%*%c)^(-0.5) %*% d.psi)
n <- nrow(X)
p <- ncol(X)
deg <- n-p
alpha <- 0.05
t.krit <- qt(1-alpha,deg)
(power <- 1-pt(t.krit, deg, ncp))
}

n1 <- 2
power.apriori(n1) 
while(power.apriori(n1) < 0.8){n1 <- n1 + 1}
n1  #10 
n1*6 #60

