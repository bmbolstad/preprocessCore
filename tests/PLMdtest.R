

library(preprocessCore)


values <- rnorm(100)
group.labels <- sample(0:4,replace=TRUE, 100)

results <- double(10000)
ngroups <- 2


for (i in 1:10000){
       values <- rnorm(100,sd=1)
       values <- values/sd(values)
       group.labels <- sample(0:(ngroups-1),replace=TRUE, 100)
       blah <- .C("R_split_test",as.double(values), as.integer(100), as.integer(ngroups), as.integer(group.labels),double(1))
       results[i] <- blah[[5]]
}

plot(sort(results),qchisq(0:9999/10000,ngroups-1))
lm(qchisq(0:9999/10000,ngroups-1) ~ sort(results))



boxplot(values ~ group.labels,ylim=c(-2,2))



sc <- median(abs(resid(lm(values ~ 1))))/0.6745
sum((resid(lm(values ~ 1))/sc)^2)/2
sum((resid(lm(values ~ as.factor(group.labels)))/sc)^2)/2


values <- rnorm(100)
group.labels <- sample(0:4,replace=TRUE, 100)
values[group.labels == 1] <- values[group.labels == 1] + 0.4


blah <- .C("R_split_test",as.double(values), as.integer(100), as.integer(5), as.integer(group.labels),double(1))

boxplot(values ~ group.labels,ylim=c(-2,2))



library(preprocessCore)

.C("R_test_get_design_matrix",as.integer(4),as.integer(5))



chips <- as.factor(rep(c(1,2,3,4,5,6),c(5,5,5,5,5,5)))
probes <- rep(c(1,3,4,5,6),6)
       
probes[c(1,6,11)] <- 2
##probes[24 + c(8,16,24)] <- 10
probes <- as.factor(probes)


model.matrix(~ -1 + probes)%*%contr.sum(6)


probes <- rep(c(1,3,4,5,6),6)
       
probes[c(1,6,11)] <- 2
probes[c(20,25,30)] <- 7
probes <- as.factor(probes)
model.matrix(~ -1 + probes)%*%contr.sum(7)




probes <- rep(c(1,3,4,5,6),6)
       
probes[c(1,6,11)] <- 2
probes[c(5,10,15)] <- 7
probes <- as.factor(probes)
model.matrix(~ -1 + probes)%*%contr.sum(7)



probes <- rep(c(1,3,4,5,6),6)
       
probes[c(1,6,11)] <- 2
probes[1+c(1,6,11)] <- 8
probes[2+c(1,6,11)] <- 9
probes[3+c(1,6,11)] <- 10
probes[c(5,10,15)] <- 7
probes <- as.factor(probes)
model.matrix(~ -1 + probes)%*%contr.sum(10)









true.probes <- c(4,3,2,1,-1,-2,-3,-4)

true.chips  <- c(8,9,10,11,12,13)


y <- outer(true.probes,true.chips,"+")



estimate.coefficients <- function(y){


colmean <- apply(y,2,mean)

y <- sweep(y,2,FUN="-",colmean)

rowmean <- apply(y,1,mean)
y <- sweep(y,1,FUN="-",rowmean)


list(y,colmean,rowmean)
}
estimate.coefficients(y)



y <- outer(true.probes,true.chips,"+")


estimate.coefficients(y)




y2 <- sweep(y,2,FUN="-",apply(y,2,mean))



c(3.875, 2.875,  1.875,  0.875,
 -1.125, -2.125, -3.125, -4, -2.25)




cp <- rep(c(1,2,3,4,5,6),rep(8,6))
pr <- rep(c(1,2,3,4,5,6,7,8),6)


pr[c(32,40,48)] <- 9




true.probes <- c(4,3,2,1,-1,-2,-3,-4)

true.chips  <- c(8,9,10,11,12,10)


y <- outer(true.probes,true.chips,"+") + rnorm(48,0,0.1)

y[8,4:6] <- c(11,12,10)+2 + rnorm(3,0,0.1)


lm(as.vector(y) ~  -1 + as.factor(cp) + C(as.factor(pr),"contr.sum"))


matplot(y,type="l")
matplot(matrix(fitted( lm(as.vector(y) ~  -1 + as.factor(cp) +
C(as.factor(pr),"contr.sum"))),ncol=6),type="l")


library(preprocessCore)
true.probes <- c(4,3,2,1,-1,-2,-3,-4)

true.chips  <- c(8,9,10,11,12,10)

y <- outer(true.probes,true.chips,"+") + rnorm(48,0,0.25)

y[8,4:6] <- c(11,12,10)+ 2.5 + rnorm(3,0,0.25)
y[5,4:6] <- c(11,12,10)+-2.5 + rnorm(3,0,0.25)



###.C("plmd_fit_R", as.double(y), as.integer(8), as.integer(6),
###		as.integer(2), as.integer(c(1,1,1,2,2,2) - 1),
###		double(6 +2*8),
###		double(48),
###		double(48))

###matplot(matrix(.C("plmd_fit_R", as.double(y), as.integer(8), as.integer(6),
###		as.integer(2), as.integer(c(1,1,1,2,2,2) - 1),
###		double(6 +2*8),
###		double(48),
###		double(48))[[7]],ncol=6))
###		


##.Call("R_plmd_model",y,0,1.3345,as.integer(c(1,1,1,2,2,2) - 1),as.integer(2))
rcModelPLM(y)
rcModelPLMd(y,c(1,1,1,2,2,2))

###R_plmd_model(SEXP Y, SEXP PsiCode, SEXP PsiK, SEXP Groups, SEXP Ngroups)





pr[seq(3,48,8)][1:3] <- 10

y[seq(3,48,8)][1:3] <- c(8,9,10) -3 + rnorm(3,0,0.1)
lm(as.vector(y) ~  -1 + as.factor(cp) + C(as.factor(pr),"contr.sum"))

