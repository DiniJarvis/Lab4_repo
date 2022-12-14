#' A RC class to handle the Linear Regression using Ordinary Linear Algebra
#'
#' This class contains various methods to handle special functions.
#' 
#' Package Description 
#' 
#' 
#' 
#' @param formula A Formula.
#' 
#' @param data A Data frame.
#'
#' @field reg_Coef to find regression coefficients
#' @field fit_Val for the fitted values
#' @field res for residuals
#' @field dof for degrees of freedom
#' @field Sigma_square for residual variance
#' @field Var_Beta for Variance of regression coefficients
#' @field t_Beta for t-values for each coefficient
#' @field pvalue for p-values for each coefficient
#' @field parse to parse the input data
#' @field stand_res for standardised residuals for plot2
#' @field variance for variance values 
#' @return nothing
#' @import methods
#' @exportClass linreg
#' @export linreg 
#' @importFrom ggplot2 theme_linedraw theme element_blank element_text stat_summary ggtitle xlab scale_x_continuous


linreg <- setRefClass( "linreg",
                       
                       fields = list(
                         formula = "formula",
                         data = "data.frame",
                         reg_Coef = "matrix",
                         fit_Val = "matrix",
                         res = "matrix" ,
                         dof = "numeric",
                         res_Var = "numeric",
                         var_Beta = "matrix",
                         t_Beta = "matrix",
                         pvalue = "matrix",
                         parse = "character",
                         stand_res = "matrix",
                         variance = "numeric"
                       ),
                       
                       methods = list(
                         initialize = function (formula, data)
                         {
                           stopifnot(all.vars(formula) %in% colnames(data))
                           stopifnot (is.data.frame(data))
                           formula <<- formula
                           data <<- data
                           X <- model.matrix(formula, data)
                           dep_y <- all.vars(formula)[1]
                           y <- as.matrix(data[dep_y])
                           parse <<- deparse(substitute(data))
                           #Regressions coefficients
                           reg_Coef <<- solve((t(X) %*% X)) %*% t(X) %*% y
                           #X <- QR
                           #Beta <- solve(R)%*%t(Q)%*%y
                           #Fitted values
                           fit_Val <<- X %*% reg_Coef
                           #Residuals
                           res <<- y - fit_Val
                           #Degrees of freedom
                           dof <<- nrow(X) - ncol(X)
                           #Residual variance
                           res_Var <<- as.numeric((t(res) %*% res) / dof)
                           #Variance of regression coefficients
                           var_Beta <<-
                             res_Var * solve((t(X) %*% X))
                           #t-values for each coefficient
                           t_Beta <<- reg_Coef / sqrt(diag(var_Beta))
                           #p values for reg coefficients
                           pvalue <<- 2 * pt(abs(t_Beta), dof,lower.tail = FALSE)
                           #variance value
                           variance <<- round(sqrt(res_Var), 2)
                           #standardised residual for plot2
                           stand_res <<-
                             sqrt(abs((res - mean(res)) / sqrt(res_Var)))
                         },
                         
                         # Prints out the coefficients and coefficient names, similar as done by the lm class.
                         print = function() {
                           cat(paste("Call: \n"))
                           cat(paste("linreg(formula = ",format(formula), ", data = ", parse, ")\n\n", sep = ""))
                           cat(paste("Coefficients:\n"))
                           coef <- structure(as.vector(reg_Coef), names= row.names(reg_Coef))
                           own_print(coef)
                         },
                         #vector of residuals e
                         resid = function(){
                           cat("Returning vector of residuals:", "\n")
                           return(as.vector(round(res,2)))
                         },
                         pred = function(){
                           cat("Returning predicted values :", "\n")
                           return(as.vector(round(fit_Val,2)))
                         },
                         coef = function(){
                           cat("Returning coefficients as a vector:", "\n")
                           return(as.vector(round(reg_Coef,2)))
                         },
                         plot = function(){
                           library(ggplot2)
                           library(ggThemeAssist)
                           theme <-  theme(
                             plot.background = element_rect(color = "black"),
                             panel.background = element_rect(fill = "white", color = NA),
                             panel.grid.major = element_line(color = "#1c1c19", size = 0.5),
                             panel.grid.major.x = element_blank(),
                             panel.grid.minor.x = element_blank(),
                             panel.grid.major.y = element_blank(),
                             panel.grid.minor.y = element_blank(),
                             axis.line = element_line(color = "#1c1c19", size = 0.5),
                             axis.text = element_text(color = "#1c1c19", size = 6),
                             axis.ticks = element_line(color = "#38ccd6", size = 0.5),
                             axis.title.x = element_text(color = "#38ccd6", size = 14,
                                                         face = "bold"),
                             axis.title.y = element_text(color = "#38ccd6", size = 14,
                                                         face = "bold"),
                             panel.grid.minor = element_line(color = "#1c1c19", size = 5),
                             plot.caption = element_text(size = 10, hjust =
                                                           0.5),
                             plot.margin = unit(c(1.2, 1.2, 1.2, 1.2), "cm"),
                             axis.text.x = element_text(size = 8),
                             axis.text.y = element_text(size = 8)
                           )
                           
                           
                           
                           plot1 <- ggplot(data.frame(fit_Val, res), aes(y=res, x=fit_Val)) + geom_point(shape=21, size=3, colour="black", fill="white")
                           plot1 <- plot1 + theme
                           plot1 <- plot1 + stat_summary(fun=median, colour="red", geom="line", aes(group = 1))
                           
                           plot1 <- plot1 + ggtitle("Residuals vs fitted") + xlab(paste("Fitted Values \n", "linreg(", format(formula), ")"))
                           plot1 <- plot1 + ylab("Residuals")
                           plot2 <- ggplot(data.frame(fit_Val, stand_res), aes(y=stand_res, x=fit_Val)) + geom_point(alpha = 0.6, shape=21, size=3, colour="black", fill="white")
                           plot2 <- plot2 + theme
                           plot2 <- plot2 + stat_summary(fun=median, colour="red", geom="line", aes(group = 1))
                           plot2 <- plot2 + ggtitle("Scale-Location") + xlab(paste("Fitted Values \n", "linreg(", format(formula), ")"))
                           plot2 <- plot2 + ylab(expression(bold(sqrt("Standardized Residuals"))))
                           plot2 <- plot2 + scale_x_continuous(breaks = seq(0.0, 1.5, by= 0.5))
                           plotlist <- list(plot1, plot2)
                           return(plotlist)
                           
                         },
                         
                         #summary()
                         summary = function(){
                           
                           cat(paste("linreg(formula = ", format(formula), ", data = ", parse, ") :\n\n ", sep = ""))
                           a<- setNames(as.data.frame(cbind(reg_Coef,as.matrix(sqrt(diag(var_Beta))),t_Beta, formatC(pvalue, format = "e", digits = 2), p_calc(pvalue))), c("Coefficients","Standard error","t-values", "p-values", ""))
                           own_print(a)
                           cat(paste("\n\n Residual standard error: ", sqrt(res_Var), " on ", dof, " degrees of freedom: ", sep = ""))
                         }
                         
                       ))

own_print<- function(a){
  print(a)
}

p_calc = function(p_val) {
  x <- ifelse(p_val > 0.1, " ",
              (ifelse(p_val > 0.05, " . ",
                      (ifelse(p_val > 0.01, "*",
                              (ifelse(p_val > 0.001, "**","***")))))))
  return(x)
}   

