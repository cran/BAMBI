contour_vmsin <- function(par.mat, pi.mix, levels=exp(seq(-20,2,by=0.5)), ...)
{
  seq_0_2pi <- seq(0, 2*pi, length.out = 100)
  log.c.vmsin <- log_const_vmsin_all(par.mat)
  e.points <- as.matrix(expand.grid(seq_0_2pi, seq_0_2pi))
  dens <- vmsinmix_manyx(e.points, par.mat, pi.mix, log.c.vmsin)
  # dens <- apply(e.points, 1, function(x) vmsinmix(x[1], x[2], par.mat, pi.mix, log.c.vmsin))
  contour(seq_0_2pi, seq_0_2pi, matrix(dens, nrow=100), levels=levels)
}

contour_vmcos <- function(par.mat, pi.mix, levels=exp(seq(-20,2,by=0.5)), ...)
{
  seq_0_2pi <- seq(0, 2*pi, length.out = 100)
  uni_rand <- matrix(runif(2e4), ncol=2)
  log.c.vmcos <- log_const_vmcos_all(par.mat, uni_rand)
  e.points <- as.matrix(expand.grid(seq_0_2pi, seq_0_2pi))
  dens <- vmcosmix_manyx(e.points, par.mat, pi.mix, log.c.vmcos)
  # dens <- apply(e.points, 1, function(x) vmcosmix(x[1], x[2], par.mat, pi.mix, log.c.vmcos))
  contour(seq_0_2pi, seq_0_2pi, matrix(dens, nrow=100), levels=levels, ...)
}

contour_wnorm2 <- function(par.mat, pi.mix, levels=exp(seq(-20,2,by=1)), omega.2pi.mat, ...)
{
  seq_0_2pi <- seq(0, 2*pi, length.out = 100)
  log.c.wnorm2 <- log_const_wnorm2_all(par.mat)
  e.points <- as.matrix(expand.grid(seq_0_2pi, seq_0_2pi))
  dens <- wnorm2mix_manyx(e.points, par.mat, pi.mix, log.c.wnorm2, omega.2pi.mat)
  # dens <- apply(e.points, 1, function(x) wnorm2mix(x, par.mat, pi.mix, log.c.wnorm2, omega.2pi.mat))
  contour(seq_0_2pi, seq_0_2pi, matrix(dens, nrow=100), levels=levels, ...)
}
