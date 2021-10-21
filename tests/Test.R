## -----------------------
## Test
## -----------------------
x = rbinom(100, 1, 0.3)

error = rnorm(100, 0, 1)

y = 1 + x + error

model.full <- lm(y ~ x)
model.reduced <- NULL

r = 1000
method = "Known"; sigma2 = 1; multi = 'normal'; boot.type = 1; alpha = 0.05; correct = TRUE; num.cores = 1


boot.ci(model.full, multi = "normal", method = "known", sigma2 = 1, boot.type = 1)

boot.ci(model.full, multi = "normal", boot.type = 1, num.cores = 5)

boot.ci(model.full, multi = "normal", boot.type = 2, num.cores = 5)

boot.ci(model.full, multi = "normal", boot.type = 3, num.cores = 5)

boot.ci(model.full, multi = "normal", boot.type = 4, num.cores = 5)

boot.ci(model.full, multi = "none", boot.type = 4, num.cores = 5)

boot.ci(model.full, multi = "rad", boot.type = 4, num.cores = 5)

# 50 0 none 1 Chisq TRUE TRUE 0.1

boot.ci(model.full = lm(y~x), multi = 'none', boot.type = 1)


# gamma
shape = 100
r = rgamma(1000, shape = shape, rate = sqrt(shape))
mean(r)
var(r)

hist(r)


 # -----------------

shape = 100
x = matrix(1:9, nrow = 3)
set.seed(111)
rgamma(9, shape = shape, rate = sqrt(shape/ x[,TRUE]^2) )
set.seed(111)
rgamma(9, shape = shape, rate = sqrt(shape/ x[,TRUE]^2)[1:9] )

solve(matrix(rep(1, 9), nrow = 3) )


# --------------
# TEST
# --------------
n = 100
m1 = m0 = 3
rhosq = 0.6
S = 0.6
m = m1 + m0

method = "Known"
# the Cov matrix of X
I_3 = diag(3)
A = matrix(rhosq/(m1*m0), nrow = m0, ncol = m1)
Sigma_x = rbind(cbind(I_3, A),
                cbind(A, I_3))


x = MASS::mvrnorm(n = n, mu = rep(0, 6), Sigma = Sigma_x)

shape = 100
errors = rgamma(n, shape =  shape, rate = sqrt(shape) ) - sqrt(shape)

beta = rep(sqrt(S^2/m1/(1-rhosq)), m)
y = 1 + x %*% beta + errors

x = as.data.frame(x)

model.full <- lm(y ~ ., data = x)
model.reduced <- lm(y ~ V1 + V2 + V3, data = x)

boot.ci(model.full, model.reduced, multi = "none", boot.type = 1, num.cores = 5)


x1 = rnorm(100, 4, 1)
x2 = rnorm(100, 6, 1)
e = rnorm(100)
outcome = 1 + 0.5*x1 + x2 + e

model.full = glm(outcome ~ x1 + x2)
model.reduced  = NULL



lm(y ~ x)


library(magrittr)

boot.model.full <- lm(form.full, data = boot.data)

boot.model.reduced <- lm(as.formula(paste0(all.vars(form.full)[1], " ~ 1")), data = boot.data)

car::Anova(boot.model.full, test.statistics = "Wald", white.adjust = FALSE)

mod1 = glm(form.full, data = boot.data)
mod2 = lm(form.full, data = boot.data)

car::Anova(mod1)
car::Anova(mod2)

mod2 = glm(as.formula(paste0(all.vars(form.full)[1], " ~ 1")), data = boot.data)
car::Anova(mod1, test.statistics = "Wald", white.adjust = TRUE)

mdscore::wald.test(boot.model.full)


mod1 = glm(outcome ~ x1)
mod2 = glm(outcome ~ x1 + x2)
mod3 = glm(outcome ~ x2)
# effect of x2
lmtest::waldtest(mod1, mod2, vcov = sandwich::vcovHC, test = "Chisq")
# effect of x1
lmtest::waldtest(mod3, mod2, vcov = sandwich::vcovHC, test = "Chisq")

# manual: x2
beta = coef(mod2)
var = sandwich::vcovHC(mod2)
naive = vcov(mod2)
wald = beta[3]^2/var[3, 3]

# manual x1
wald = beta[2]^2/var[2, 2]



# structure
title = "ANOES based on the RESI"

structure()

match.arg(method)

# Test: boot.ci
n = 100
x1 = rnorm(n, 4, 1)
x2 = rnorm(n, 6, 1)
e = rnorm(n)
outcome = 1 + 0.5*x1 + 3*x2 + e

## linear models
model1 = glm(outcome ~ x1  +  x2)
model0  = glm(outcome ~ x1)

boot.ci(model.full = model1, model.reduced = model0)

boot.ci(model.full = model1)

model0 = glm(outcome ~ 1)

boot.ci(model.full = model1, model.reduced = model0)


# logistic regression model

n = 100
x1 = rbinom(n, 1, 0.3)
x2 = rnorm(n, 0, 2)

p = 1/(1+exp(-c(1 + 2*x1 - 1*x2)))

y = rbinom(n, 1, p)

mod1 = glm(y ~ x1 + x2, family = "binomial")

data = mod1$model
data$resid = residuals(mod1, type = "response")

data$hat = predict(mod1, data, type = "response")
data$new_hat = data$hat + data$resid

boot.ci(mod1)



