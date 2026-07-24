#Get q,x,w,mu true values via iterations

#specify true values of variance components
a<-sqrt(0.2)
d<-sqrt(0.2)
f<-1
s<-sqrt(0.2)
t<-sqrt(0)
e<-sqrt(0.2)

VF<-0.2
mat2pat<-1 #maternal relative to the paternal environmental transmission effect
K <- mat2pat^2 + 1
par0<-a*a+d*d+s*s+t*t+VF+e*e #the variance of the latent parenting trait in gen0
n<-sqrt(VF/(K*par0)) #paternal path coefficient 
m <- n*mat2pat #maternal path coefficient 
  
numits <- 500

DAT <- data.frame(q=rep(NA,numits),x=rep(NA,numits),w=rep(NA,numits),mu=rep(NA,numits))

#Gen 0 values
q<-1
x<-.2
w<-0
mu<-.15/(a^2*q+f^2*x+2*a*w*f+d^2+s^2+t^2+e^2)

#iterations
for (IT in 1:numits){
q_new<-1+mu*(q*a+w*f)^2
x_new<-(n^2+m^2)*(a^2*q+f^2*x+2*a*w*f+d^2+s^2+t^2+e^2)+2*m*n*mu*(a*q*a+f*x*f+a*w*f+a*w*f+d*d+s*s+t*t+e*e)^2
w_new<-.5*mu*(m+n)*(q*a+w*f)*(a*q*a+f*x*f+a*w*f+a*w*f+d*d+s*s+t*t+e*e)+.5*(q*a+w*f)*(m+n)

q<-q_new
x<-x_new
w<-w_new

mu_new<-.15/(a^2*q+f^2*x+2*a*w*f+d^2+s^2+t^2+e^2) 
mu<-mu_new

DAT[IT,] <- c(q,x,w,mu)
}

#equilibrium values for q,x,w,mu
DAT[500,]
