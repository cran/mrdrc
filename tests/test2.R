#library(mrdrc)
#
#for (i in unique(exp.az$exp))
#{
#    exp.azSub <- subset(exp.az, exp == i)
#
#    exp.azSub.m <- pnsdrm(exp.azSub$conc, exp.azSub$response, type = "continuous",
#    model = "semi-parametric", fct = list(LL.4()), respLev = runif(5, 1, 50),
#    reference = NULL, level = 0.95, robust = FALSE, logex = TRUE)
#
#     exp.azSub.m
#
#    plot(exp.azSub.m) 
#}
