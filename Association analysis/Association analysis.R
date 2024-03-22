
##### Logistic Model
model <- glm( pheno ~ protein + covariates, data=data, family="binomial")

##### Linear Model
model <- lm( pheno ~ protein + covariates, data=data)