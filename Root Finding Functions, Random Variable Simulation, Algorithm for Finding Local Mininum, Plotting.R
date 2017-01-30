```{r}

#----------------- finds root of a function using the Newton Raphson method ------

newtonraphson_show <- function(ftn, x0, xmin = x0-1, xmax = x0+1) {
  # applies Newton-Raphson to find x such that ftn(x)[1] == 0
  # x0 is the starting point
  # subsequent iterations are plotted in the range [xmin, xmax]

  # plot the function
  x <- seq(xmin, xmax, (xmax - xmin)/200)
  fx <- c()
  for (i in 1:length(x)) {
    fx[i] <- ftn(x[i])[1]
  }
  plot(x, fx, type = "l", xlab = "x", ylab = "f(x)",
    main = "zero f(x) = 0", col = "blue", lwd = 2)
  lines(c(xmin, xmax), c(0, 0), col = "blue")

  # do first iteration
  xold <- x0
  f.xold <- ftn(xold)
  xnew <- xold - f.xold[1]/f.xold[2]
  lines(c(xold, xold, xnew), c(0, f.xold[1], 0), col = "red")

  cat("last x value", xnew, "\n")
  for(i in 1:3){
    xold <- xnew;
    f.xold <- ftn(xold)
    xnew <- xold - f.xold[1]/f.xold[2]
    lines(c(xold, xold, xnew), c(0, f.xold[1], 0), col = "red")
    cat("last x value", xnew, "\n")
  }

  return(xnew)
}


a <- function(x){
  fx <- cos(x) - x
  dfx <- -sin(x) - 1
  return(c(fx,dfx))
} 

newtonraphson_show(a, 3, -1,3.1)

b <- function(x){
  fx <- log(x) - exp(-x)
  dfx <- 1/x + exp(-x)
  return(c(fx,dfx))
}

newtonraphson_show(b, 2, 1, 2.2)
```

```{r}

## ------------------------Finds root of a function using the secant method--------

par(mfrow = c(1,1))
secant_show <- function(ftn, x0, x1, xmin = min(x0, x1) - 1, xmax = max(x0, x1) + 1){
  # applies the secant algorithm to find x such that ftn(x) == 0
  # x0 and x1 are the starting points

  # plot the function
  x <- seq(xmin, xmax, (xmax - xmin)/200)
  fx <- sapply(x, ftn)
  plot(x, fx, type = "l", xlab = "x", ylab = "f(x)",
       main = "zero of f(x)", col = "yellow", lwd = 2)
  lines(c(xmin, xmax), c(0, 0), col = "black")
  
  # calculate your initial f0 and f1
  f0 <- ftn(x0)
  f1 <- ftn(x1)
  
  xolder <- x0
  xold <- x1
  f.xolder <- ftn(xolder)
  f.xold <- ftn(xold)
  xnew <- xold - f.xold*(xold - xolder)/(f.xold - f.xolder)
  lines(c(x0, x1, xnew), c(f0,f1,0), col = "red")

    # do first iteration
  
  cat("last x value", xnew, "\n")
  
  for (i in 1:3) {
    xolder <- xold
    xold <- xnew
    f.xolder <- ftn(xolder)
    f.xold <- ftn(xold)
    xnew <- xold - f.xold*(xold - xolder)/(f.xold - f.xolder)
    lines(c(xolder, xold, xnew), c(f.xolder,f.xold,0), col = "red")
    cat("last x value", xnew, "\n")
  }
  
  
  cat("last x value", xnew, "\n")

  return(xnew)
}

f1 <- function(x) cos(x) - x
secant_show(f1, 3, 2.8, -1, 3.5)

```

```{r}
## ------------ Simulating and plotting a Random Variable ---------------

cmf <- function(x){
  Fx <- 0
  if (x >= 1) {
    Fx <- (0.1)
    }
  
  if (x >= 2){
    Fx <- Fx + (0.3)}
  
  if (x >=5){
    Fx <- Fx + (0.6)
  }
  return(Fx)
}

sim <- function(F){
  U = runif(1)
  X <- 0 
  while(F(X)<U){
    X <- X+1
  }
  return(X)
}
x <- c(1,2,3,4,5)
Fx1 <- cmf(1)
Fx2 <- cmf(2)
Fx3 <- cmf(3)
Fx4 <- cmf(4)
Fx5 <- cmf(5)

plot(x = x, y = c(Fx1, Fx2, Fx3, Fx4, Fx5), type = "s", col = "red")

set.seed(3)
sim(cmf) 
sim(cmf)
sim(cmf)

## ------------ Simulating a Negative Binomial Random Variable -------

cmf_nb <- function(x, r, p){
Fx <- 0
    for (i in 0:x){
  Fx <- Fx + choose(i+r-1, i) * p^i * (1-p)^r  
  }
  return(Fx)
}  

sim_nb <- function(F, ...){
  U <- runif(1)
  X <- 0
  while(F(X, ...) < U){
    X = X + 1
  }
  return(X)
}

set.seed(3)
sim_nb(cmf_nb, 3, .5) 
sim_nb(cmf_nb, 3, .5)
sim_nb(cmf_nb, 3, .5)
```

```{r}
# ----Coordinate Descent Algorithm for Optimization to Find Local Minimum of a Function---

##### golden function is a modification of code provided by Eric Cai
golden = function(f, lower, upper, tolerance = 1e-5)
{
   golden.ratio = 2/(sqrt(5) + 1)

   ## Use the golden ratio to find the initial test points
   x1 <- lower + golden.ratio * (upper - lower)
   x2 <- upper - golden.ratio * (upper - lower)
   
   ## the arrangement of points is:
   ## lower ----- x2 --- x1 ----- upper

   ### Evaluate the function at the test points
   f1 <- f(x1)
   f2 <- f(x2)

   while (abs(upper - lower) > tolerance) {
        if (f2 > f1) {
        # the minimum is to the right of x2
        lower <- x2  # x2 becomes the new lower bound
        x2 <- x1     # x1 becomes the new x2
        f2 <- f1     # f(x1) now becomes f(x2)
        x1 <- lower + golden.ratio * (upper - lower)  
        f1 <- f(x1)  # calculate new x1 and f(x1)
        } else {
        # then the minimum is to the left of x1
        upper <- x1  # x1 becomes the new upper bound
        x1 <- x2     # x2 becomes the new x1
        f1 <- f2
        x2 <- upper - golden.ratio * (upper - lower)
        f2 <- f(x2)  # calculate new x2 and f(x2)
        }
    }
    (lower + upper)/2 # the returned value is the midpoint of the bounds
}

f <- function(x){ (x - 3)^2 }
golden(f, 0, 10)

g <- function(x,y) { 
    5 * x ^ 2 - 6 * x * y + 5 * y ^ 2
    }
x <- seq(-1.5,1.5, len=100)
y <- seq(-1.5,1.5, len=100)
z <- outer(x,y,g)



coord_descent <- function(g, start_x, start_y, lower, upper, tol = 1e-5, iterations = 15){
  
  contour(x,y,z, levels = seq(.5,5,by=.9)) # code to plot contour lines
  x_i <- -1.5
  y_i <- -1.5
  
  x <- start_x
  y <- start_y

  g_x <- function(x){
    g(x, y = start_y)
  }
  next_x <- golden(g_x,lower,upper)
  x <- next_x
  
  lines(c(start_x, next_x), c(start_y, start_y), col = c(2))
  
  g_y <- function(y){
   g(x = next_x,y)
  }
  next_y <- golden(g_y,lower,upper)
  y <- next_y
  
  lines(c(next_x, next_x), c(start_y, next_y), col = c(2))
  
  dif <- g(start_x,start_y) - g(next_x,next_y)
  iter <- 1

  while (abs(dif) > tol & iter < iterations) {
    iter <- iter + 1
    current <- c(x,y)
    print(current)
    previous_x <- next_x
    previous_y <- next_y
    
    g_x <- function(x){
      g(x, y = next_y)
    }
    next_x <- golden(g_x, lower, upper)
    x <- next_x
    
    lines(c(previous_x, next_x), c(previous_y, previous_y), col = c((iter+1):iterations))
    
    g_y <- function(y){
      g(x = next_x,y)
    }
    next_y <- golden(g_y, lower, upper)
    y <- next_y
   
    lines(c(next_x, next_x), c(previous_y, next_y),col = c((iter+1):iterations))
    
   dif <- g(previous_x,previous_y) - g(next_x, next_y)
  }
  
  minimum <- c(next_x, next_y)
  return(minimum)
}


coord_descent(g, -1.5, -1.5, 1.5, -1.5)

```