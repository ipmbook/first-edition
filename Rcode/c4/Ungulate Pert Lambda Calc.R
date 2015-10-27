## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Use the Ungulate IBM to illustrate the construction of an IPM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls(all=TRUE))

## R's working directory must be set to where this file sits, so the source()'s below run

## run the utility functions
source("../utilities/Standard Graphical Pars.R")

## run the ungulate IBM, fit demographic models & build the IPM
# source("../c2/Ungulate Calculations.R")

source("../c2/Ungulate Demog Funs.R")
min.size <- 1.60
max.size <- 3.70

## set the working directory to our figures location
# setwd("~/Repos/ipm_book/c4/figures")

## copy model parameters 
m.par <- m.par.true

## build the IPM associated with the  
IPM.sys <- mk_K(200, m.par, min.size, max.size)

## extract the mesh points and the kernel
meshpts <- IPM.sys$meshpts
h <- diff(meshpts[1:2])
K <- IPM.sys$K
P <- IPM.sys$P
F <- IPM.sys$F

## compute the eigen vectors / values 
IPM.eig.sys <- eigen(K)
## lambda
lambda <- Re(IPM.eig.sys$values[1])
## unnormalised stable size distribution ('w')
w.z <- Re(IPM.eig.sys$vectors[,1])
## ... as before (nothing new here)

## unnormalised reproductive value distribution ('v')...
v.z1 <- Re(eigen(t(K))$vectors[,1])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## 1. kernel sensitivity and elasticity
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

## kernel sensitivity
K.sens <- outer(v.z1, w.z, "*")/sum(v.z1 * w.z * h)
## kernel elasticity (the divide-by-dz is needed to ensure we are working
## with the kernel, which is a function, and not the iteration matrix)
K.elas <- K.sens * (K/h) / lambda
## does the kernel elasticity integrate to 1?...
sum(K.elas) * h^2
## ...yes

## calculate the elasticities of the survival and reproduction components
P <- IPM.sys$P / h
F <- IPM.sys$F / h
P.elas <- P * K.sens / lambda
F.elas <- F * K.sens / lambda

## get the relative contributions of survival and reproduction to lambda
sum(P.elas) * h^2
sum(F.elas) * h^2

## plot these
ikeep <- which(meshpts>2.3 & meshpts<3.5) # use to extract a region to plot
postscript("KernSensElas.eps", 
           width=8, height=8, horizontal=FALSE, paper="special")
## plot the sensitivity kernel and the three elasticity surfaces
set_graph_pars("panel4")
image(meshpts[ikeep], meshpts[ikeep], t(K.sens[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(K.sens[ikeep,ikeep]), add=TRUE)
add_panel_label("a")
image(meshpts[ikeep], meshpts[ikeep], t(K.elas[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(K.elas[ikeep,ikeep]), 
        add=TRUE, levels=seq.int(1, 21, by=2))
add_panel_label("b")
image(meshpts[ikeep], meshpts[ikeep], t(P.elas[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(P.elas)[ikeep,ikeep], 
        add=TRUE, levels=seq.int(1, 21, by=2))
add_panel_label("c")
image(meshpts[ikeep], meshpts[ikeep], t(F.elas[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(F.elas)[ikeep,ikeep], add=TRUE)
add_panel_label("d")
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## 2. calculate survival & probability of reproduction function perturbations 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## survival
dK_by_ds_z1z <-
    outer(meshpts, meshpts,
          function (z1, z, m.par)
              {g_z1z(z1, z, m.par) + pb_z(z, m.par) * (1/2) * pr_z(m.par) * c_z1z(z1, z, m.par)},
          m.par)
s.sens.z <- apply(K.sens * dK_by_ds_z1z, 2, sum) * h
s.elas.z <- s.sens.z * s_z(meshpts, m.par) / lambda

## probability of reproduction
dK_by_dpb_z1z <- outer(meshpts, meshpts,
                      function(z1, z, m.par)
                        {s_z(z, m.par) * (1/2) * pr_z(m.par) * c_z1z(z1, z, m.par)},
                      m.par)
pb.sens.z <- apply(K.sens * dK_by_dpb_z1z, 2, sum) * h
pb.elas.z <- pb.sens.z * pb_z(meshpts, m.par) / lambda

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. plot the survival / reproduction perturbation functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ikeep <- which(meshpts>2.0 & meshpts<3.6) # use to extract a region to plot
postscript("SurvReprSensElas.eps", 
           width=8, height=4.4, horizontal=FALSE, paper="special")
set_graph_pars("panel2")
plot(meshpts[ikeep], s.sens.z[ikeep], type="l", xlab="Mass, z", 
     ylab=expression(paste("Sensitivity / Elasticity of ",italic(s),"(",italic(z),")")))
lines(meshpts[ikeep], s.elas.z[ikeep], lty=2)
legend("topleft", legend=c("Sensitivity","Elasticity"), lty=c(1,2), bty="n")
add_panel_label("a")
plot(meshpts[ikeep], pb.sens.z[ikeep], type="l", xlab="Mass, z", 
     ylab=expression(paste("Sensitivity / Elasticity of ",italic(p)[italic(b)],"(",italic(z),")")))
lines(meshpts[ikeep], pb.elas.z[ikeep], lty=2)
legend("topleft", legend=c("Sensitivity","Elasticity"), lty=c(1,2), bty="n")
add_panel_label("b")
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 4. calculate the growth & offspring size kernel perturbation functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## growth 
dK_by_dg_z1z <- outer(meshpts, meshpts,
                      function(z1, z, m.par) s_z(z, m.par), m.par)
g.sens.z1z <- K.sens * dK_by_dg_z1z 
g.elas.z1z <- g.sens.z1z * outer(meshpts, meshpts, g_z1z, m.par) / lambda

## offspring size
dK_by_dc_z1z <- outer(meshpts, meshpts,
                      function(z1, z, m.par)
                        s_z(z, m.par) * pb_z(z, m.par) * (1/2) * pr_z(m.par),
                      m.par)
c.sens.z1z <- K.sens * dK_by_dc_z1z
c.elas.z1z <- c.sens.z1z * outer(meshpts, meshpts, c_z1z, m.par) / lambda

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. plot growth & offspring size kernel perturbation functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## set up the plots
ikeep <- which(meshpts>2.3 & meshpts<3.5) # use to extract a region to plot
postscript("GrowOffSizeSensElas.eps", 
           width=8, height=8, horizontal=FALSE, paper="special")
set_graph_pars("panel4")
## plot the growth sensitivity and elasticity surfaces
image(meshpts[ikeep], meshpts[ikeep], t(g.sens.z1z[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(g.sens.z1z[ikeep,ikeep]), add=TRUE)
add_panel_label("a")
image(meshpts[ikeep], meshpts[ikeep], t(g.elas.z1z[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(g.elas.z1z[ikeep,ikeep]), 
        add=TRUE, levels=seq.int(1,21,by=2))
add_panel_label("b")
## plot the offspring size kernel sensitivity and elasticity surfaces
image(meshpts[ikeep], meshpts[ikeep], t(c.sens.z1z[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(c.sens.z1z[ikeep,ikeep]), add=TRUE)
add_panel_label("c")
image(meshpts[ikeep], meshpts[ikeep], t(c.elas.z1z[ikeep,ikeep]),
      col=grey(seq(0.6, 1, length=100)),
      xlab="Mass (t), z", ylab="Mass (t+1), z\'")
contour(meshpts[ikeep], meshpts[ikeep], t(c.elas.z1z[ikeep,ikeep]), add=TRUE)
add_panel_label("d")
dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 6. total elasticity/sensitivity (+ sanity checks)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## sensitivity / elasticity associated with perturbing the whole function 
sum(s.sens.z * h)
sum(s.elas.z * h)
## survival elasiticity == 1  
## ...because the whole kernel is scaled by 's_z'

## sensitivity / elasticity associated with perturbing the whole function 
sum(pb.sens.z * h)
sum(pb.elas.z * h)
## probability reproduction elasticity != 1  
## ...because only the F kernel is scaled by 'p_b'

## integrated growth kernel elasticity !=1
sum(g.sens.z1z) * h^2
sum(g.elas.z1z) * h^2
## integrated offspring size kernel elasticity !=1
sum(c.sens.z1z) * h^2
sum(c.elas.z1z) * h^2

## integrated elasticities associated with the fecundity component functions are equal 
round(sum(pb.elas.z)*h, 8) == round(sum(c.elas.z1z)*h^2, 8)
## integrated growth and mother-daughter function elasticities sum to 1
sum(g.elas.z1z)*h^2 + sum(c.elas.z1z)*h^2

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 7. parameter sensitivities / elasticities
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## survival intercept - 1) use of equ ???? 
ds_by_dbeta0_z1z <-
    outer(meshpts, meshpts,
          function(z1, z, m.par) {
              nu <- m.par["surv.int"] + m.par["surv.z"] * z 
              exp(nu)/(1+exp(nu))^2
          }, m.par)

s.int.sens <- sum(K.sens * dK_by_ds_z1z * ds_by_dbeta0_z1z) * h^2
s.int.elas <- s.int.sens * m.par["surv.int"] / lambda

## survival size slope
ds_by_dbetaz_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          nu <- m.par["surv.int"] + m.par["surv.z"] * z 
          z * exp(nu)/(1+exp(nu))^2
        }, m.par)

s.slp.sens <- sum(K.sens * dK_by_ds_z1z * ds_by_dbetaz_z1z) * h^2
s.slp.elas <- s.slp.sens * m.par["surv.z"] / lambda

## probability of reproduction intercept
dpb_by_dbeta0_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          nu <- m.par["repr.int"] + m.par["repr.z"] * z 
          exp(nu)/(1+exp(nu))^2
        }, m.par)

pb.int.sens <- sum(K.sens * dK_by_dpb_z1z * dpb_by_dbeta0_z1z) * h^2
pb.int.elas <- s.int.sens * m.par["repr.int"] / lambda

## probability of reproduction size slope
dpb_by_dbetaz_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          nu <- m.par["repr.int"] + m.par["repr.z"] * z 
          z * exp(nu)/(1+exp(nu))^2
        }, m.par)

pb.slp.sens <- sum(K.sens * dK_by_dpb_z1z * dpb_by_dbetaz_z1z) * h^2
pb.slp.elas <- pb.slp.sens * m.par["repr.z"] / lambda

## growth intercept
dg_by_dbeta0_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          g_z1z(z1, z, m.par) * 
            (z1 - m.par["grow.int"] - m.par["grow.z"] * z) / m.par["grow.sd"]^2
        }, m.par)

g.int.sens <- sum(K.sens * dK_by_dg_z1z * dg_by_dbeta0_z1z) * h^2
g.int.elas <- g.int.sens * m.par["grow.int"] / lambda

## growth size slope
dg_by_dbetaz_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          g_z1z(z1, z, m.par) * 
            z * (z1 - m.par["grow.int"] - m.par["grow.z"] * z) / m.par["grow.sd"]^2
        }, m.par)

g.slp.sens <- sum(K.sens * dK_by_dg_z1z * dg_by_dbetaz_z1z) * h^2
g.slp.elas <- g.slp.sens * m.par["grow.z"] / lambda

## growth standard deviation
dg_by_dsigma_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          nu <- m.par["grow.int"] + m.par["grow.z"] * z
          g_z1z(z1, z, m.par) * ((z1 - nu)^2 - m.par["grow.sd"]^2) / m.par["grow.sd"]^3
        }, m.par)

g.sig.sens <- sum(K.sens * dK_by_dg_z1z * dg_by_dsigma_z1z) * h^2
g.sig.elas <- g.sig.sens * m.par["grow.sd"] / lambda

## offspring size intercept
dc_by_dbeta0_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          c_z1z(z1, z, m.par) * 
            (z1 - m.par["rcsz.int"] - m.par["rcsz.z"] * z) / m.par["rcsz.sd"]^2
        }, m.par)

c.int.sens <- sum(K.sens * dK_by_dc_z1z * dc_by_dbeta0_z1z) * h^2
c.int.elas <- c.int.sens * m.par["rcsz.int"] / lambda

## offspring size slope
dc_by_dbetaz_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          c_z1z(z1, z, m.par) * 
            z * (z1 - m.par["rcsz.int"] - m.par["rcsz.z"] * z) / m.par["rcsz.sd"]^2
        }, m.par)

c.slp.sens <- sum(K.sens * dK_by_dc_z1z * dc_by_dbetaz_z1z) * h^2
c.slp.elas <- c.slp.sens * m.par["rcsz.z"] / lambda

## growth standard deviation
dc_by_dsigma_z1z <-
  outer(meshpts, meshpts,
        function(z1, z, m.par) {
          nu <- m.par["rcsz.int"] + m.par["rcsz.z"] * z
          c_z1z(z1, z, m.par) * ((z1 - nu)^2 - m.par["rcsz.sd"]^2) / m.par["rcsz.sd"]^3
        }, m.par)

c.sig.sens <- sum(K.sens * dK_by_dc_z1z * dc_by_dsigma_z1z) * h^2
c.sig.elas <- c.sig.sens * m.par["rcsz.sd"] / lambda

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 8. perturbation functions associated with the expected size functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## expected growth function
mu.g.sens.z <- apply(K.sens * dK_by_dg_z1z * dg_by_dbeta0_z1z, 2, sum) * h
mu.g <- m.par["grow.int"] + m.par["grow.z"] * meshpts
mu.g.elas.z <- mu.g.sens.z * mu.g / lambda

## expected size of offspring function
mu.c.sens.z <- apply(K.sens * dK_by_dc_z1z * dc_by_dbeta0_z1z, 2, sum) * h
mu.c <- m.par["rcsz.int"] + m.par["rcsz.z"] * meshpts
mu.c.elas.z <- mu.c.sens.z * mu.c / lambda

## plot 
ikeep <- which(meshpts>2.0 & meshpts<3.6) # use to extract a region to plot
postscript("MeanGrowOffSizeSensElas.eps", 
           width=8, height=4.4, horizontal=FALSE, paper="special")
set_graph_pars("panel2")
plot(meshpts[ikeep], mu.g.sens.z[ikeep], type="l", ylim=c(0,12), xlab="Mass, z", 
     ylab=expression(paste("Sensitivity / Elasticity of ",italic(mu[g]),"(",italic(z),")")))
lines(meshpts[ikeep], mu.g.elas.z[ikeep], lty=2)
legend("topleft", legend=c("Sensitivity","Elasticity"), lty=c(1,2), bty="n")
add_panel_label("a")
plot(meshpts[ikeep], mu.c.sens.z[ikeep], type="l", ylim=c(0,2.5), xlab="Mass, z", 
     ylab=expression(paste("Sensitivity / Elasticity of ",italic(mu)[italic(c[0])],"(",italic(z),")")))
lines(meshpts[ikeep], mu.c.elas.z[ikeep], lty=2)
legend("topleft", legend=c("Sensitivity","Elasticity"), lty=c(1,2), bty="n")
add_panel_label("b")
dev.off()



