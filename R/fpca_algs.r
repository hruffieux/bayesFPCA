vmp_gauss_fpca <- function(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, n_mc=100, plot_elbo=FALSE
) {

	# Establish necessary parameters:

	Sigma_zeta <- sigma_zeta^2*diag(L)
	T_vec <- sapply(Y, length)
	K <- dim(C[[1]])[2] - 2
	d <- (K+2)*(L+1)

	# Initialise VMP simulation:

	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {

		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Sigma_q_zeta[[i]] <- diag(L)
	}

	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)

	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)

	eta_1_sum <- 0
	eta_2_sum <- 0
	for(i in 1:N) {

		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
		eta_1_sum <- eta_1_sum + sum_val

		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- 1*eta_1_sum
	eta_2 <- -1/2*eta_2_sum
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)

	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, mu_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)

	eta_1 <- -0.5*sum(T_vec)
	eta_2 <- -0.5*sum(T_vec)
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"

	eta_vec$"p(zeta)->zeta" <- replicate(
		N, gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
	)

	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"

	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
	G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"

	eta_1 <- rep(0, d)
	eta_2 <- -0.5*as.vector(diag(d))
	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)

	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"

	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)

	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"

	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)

	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
	G$"p(sigsq_m|a_m)->a_m" <- "diag"

	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)

	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

	eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
	G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]

	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]

	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)

	elbo_res <- NULL
	converged <- FALSE
	iter <- 0

	while((!converged) & (iter < n_vmp)) {

		iter <- iter + 1

		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"

		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"

		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"

		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"

		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"

		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"

		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"

		# Update p(Y|nu,zeta,sigsq_eps) fragment:

		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)

		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)

		fpc_lik_fragment <- fpc_lik_frag(
			eta_in, G_in, C, Y, T_vec, L
		)

		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- fpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- fpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"eta"[[3]]

		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"G"[[1]]

		# Update p(nu|Sigma_nu) fragment:

		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)

		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)

		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
			eta_in, G_in, L, mu_beta, Sigma_beta
		)

		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]

		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]

		# Update p(sigsq_eps|a_eps) fragment:

		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)

		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
			1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
		)

		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]

		G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]

		# Update p(sigsq_m|a_m) fragment:

		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)

		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_m->p(sigsq_m|a_m)",
			1, G$"sigsq_m->p(sigsq_m|a_m)"
		)

		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]

		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]

		# Update p(sigsq_p|a_p) fragment:

		for(l in 1:L) {

			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)

			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
			)

			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]

			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
		}

		# Compute the entropy:

		ent <- 0

		eta_nu <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(nu|Sigma_nu)->nu"
		)
		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)

		ent <- ent + ent_nu

		ent_zeta <- 0
		for(i in 1:N) {

			eta_zeta <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
			ent_zeta <- ent_zeta + sum_val
		}

		ent <- ent + ent_zeta

		eta_sigsq_eps <- list(
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		G_sigsq_eps <- c(
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
			G$"p(sigsq_eps|a_eps)->sigsq_eps"
		)
		ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)

		ent <- ent + ent_sigsq_eps

		eta_a_eps <- list(
			eta_vec$"p(sigsq_eps|a_eps)->a_eps",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_a_eps <- c(
			G$"p(sigsq_eps|a_eps)->a_eps",
			G$"p(a_eps)->a_eps"
		)
		ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)

		ent <- ent + ent_a_eps

		eta_sigsq_m <- list(
			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		)
		G_sigsq_m <- c(
			G$"p(nu|Sigma_nu)->sigsq_m",
			G$"p(sigsq_m|a_m)->sigsq_m"
		)
		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)

		ent <- ent + ent_sigsq_m

		eta_a_m <- list(
			eta_vec$"p(sigsq_m|a_m)->a_m",
			eta_vec$"p(a_m)->a_m"
		)
		G_a_m <- c(
			G$"p(sigsq_m|a_m)->a_m",
			G$"p(a_m)->a_m"
		)
		ent_a_m <- entropy_igw(eta_a_m, G_a_m)

		ent <- ent + ent_a_m

		ent_sigsq_p <- 0
		ent_a_p <- 0
		for(l in 1:L) {

			eta_sigsq_p <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
			)
			G_sigsq_p <- list(
				G$"p(nu|Sigma_nu)->sigsq_p"[l],
				G$"p(sigsq_p|a_p)->sigsq_p"[l]
			)
			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
			ent_sigsq_p <- ent_sigsq_p + sum_val

			eta_a_p <- list(
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_a_p <- c(
				G$"p(sigsq_p|a_p)->a_p"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- entropy_igw(eta_a_p, G_a_p)
			ent_a_p <- ent_a_p + sum_val
		}

		ent <- ent + ent_sigsq_p + ent_a_p

		# Compute the cross-entropy:

		c_ent <- 0

		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		c_ent_p_Y <- cross_entropy_fpc_lik_frag(eta_in, G_in, C, Y, T_vec, L)

		c_ent <- c_ent + c_ent_p_Y

		eta_in <- list(
			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		)
		G_in <- list(
			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
		)
		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)

		c_ent <- c_ent + c_ent_p_nu

		c_ent_p_zeta <- 0
		for(i in 1:N) {

			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[,l],
				eta_vec$"p(zeta)->zeta"[,l]
			)
			sum_val <- cross_entropy_gauss_prior(
				eta_in, rep(0, L),
				Sigma_zeta, use_vech=TRUE
			)
			c_ent_p_zeta <- c_ent_p_zeta + sum_val
		}

		c_ent <- c_ent + c_ent_p_zeta

		eta_in <- list(
			eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
			eta_vec$"a_eps->p(sigsq_eps|a_eps)",
			eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		)
		G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
		G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
		c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

		c_ent <- c_ent + c_ent_p_sigsq_eps

		eta_in <- list(
			eta_vec$"a_eps->p(a_eps)",
			eta_vec$"p(a_eps)->a_eps"
		)
		G_in <- c(
			G$"a_eps->p(a_eps)",
			G$"p(a_eps)->a_eps"
		)
		c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

		c_ent <- c_ent + c_ent_p_a_eps

		eta_in <- list(
			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
			eta_vec$"a_m->p(sigsq_m|a_m)",
			eta_vec$"p(sigsq_m|a_m)->a_m"
		)
		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

		c_ent <- c_ent + c_ent_p_sigsq_m

		eta_in <- list(
			eta_vec$"a_m->p(a_m)",
			eta_vec$"p(a_m)->a_m"
		)
		G_in <- c(
			G$"a_m->p(a_m)",
			G$"p(a_m)->a_m"
		)
		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

		c_ent <- c_ent + c_ent_p_a_m

		c_ent_p_sigsq_p <- 0
		c_ent_p_a_p <- 0
		for(l in 1:L) {

			eta_in <- list(
				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
			)
			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val

			eta_in <- list(
				eta_vec$"a_p->p(a_p)"[,l],
				eta_vec$"p(a_p)->a_p"[,l]
			)
			G_in <- c(
				G$"a_p->p(a_p)"[l],
				G$"p(a_p)->a_p"[l]
			)
			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
			c_ent_p_a_p <- c_ent_p_a_p + sum_val
		}

		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p

		# Compute the ELBO

		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)

		if(plot_elbo) {

			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
		}

		if(iter > 1) {

			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)

			if(iter > 2) {

				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {

				criterion_2_satisfied <- FALSE
			}

			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)

			if(criterion_satisfied) {

				converged <- TRUE
			}
		}
	}

	# Get the list of natural parameter vectors:

	return(eta_vec)
}

#' @export
vmp_gauss_mfpca <- function(
	n_vmp, N, p, L, C, Y, mu_beta,
	Sigma_beta, A, time_g, C_g,
	criterion, plot_elbo = FALSE
) {

	# Establish necessary parameters:

	n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))
	K <- ncol(C[[1]][[1]]) - 2
	d <- (K+2)*(L+1)

	# Initialise VMP simulation:

	E_q_zeta <- vector("list", length = N)
	Cov_q_zeta <- vector("list", length = N)
	for(i in 1:N) {

		E_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
		Cov_q_zeta[[i]] <- diag(L)
	}

	eta_vec <- vector("list", length=32)
	names(eta_vec) <- c(
		"nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
		"zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"zeta->p(zeta)", "p(zeta)->zeta",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)

	G <- vector("list", length=24)
	names(G) <- c(
		"sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
		"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
		"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
		"a_eps->p(a_eps)", "p(a_eps)->a_eps",
		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
		"a_m->p(a_m)", "p(a_m)->a_m",
		"a_p->p(a_p)", "p(a_p)->a_p"
	)

	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- vector("list", length = p)
	eta_vec$"p(nu|Sigma_nu)->nu" <- vector("list", length = p)
	for(j in 1:p) {

		eta_1_sum <- 0
		eta_2_sum <- 0
		for(i in 1:N) {

			E_q_zeta_tilde <- c(1, E_q_zeta[[i]])
			Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta[[i]])
			E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

			sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
			eta_1_sum <- eta_1_sum + sum_val

			sum_val <- as.vector(kronecker(E_q_tcross_zeta_tilde, crossprod(C[[i]][[j]])))
			eta_2_sum <- eta_2_sum + sum_val
		}

		eta_1 <- eta_1_sum
		eta_2 <- -0.5*eta_2_sum
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"[[j]] <- c(eta_1, eta_2)

		eta_1 <- rep(0, d)
		eta_2 <- -0.5*as.vector(diag(d))
		eta_vec$"p(nu|Sigma_nu)->nu"[[j]] <- c(eta_1, eta_2)
	}

	D_L <- duplication.matrix(L)
	eta_1 <- Reduce(cbind, E_q_zeta)
	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L))))
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)

	eta_vec$"p(zeta)->zeta" <- replicate(
		N, gauss_prior_frag(rep(0, L), diag(L), use_vech = TRUE)
	)

	eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- vector("list", length = p)
	G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- vector("list", length = p)
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- vector("list", length = p)
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- vector("list", length = p)
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- vector("list", length = p)
	G$"p(sigsq_eps|a_eps)->a_eps" <- vector("list", length = p)
	for(j in 1:p) {

		eta_1 <- -0.5*sum(n[, j])
		eta_2 <- -0.5*sum(n[, j])
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[[j]] <- c(eta_1, eta_2)
		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[[j]] <- "full"

		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]] <- c(-3/2, -1/2)
		G$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]] <- "full"

		eta_vec$"p(sigsq_eps|a_eps)->a_eps"[[j]] <- c(-1/2, -1/2)
		G$"p(sigsq_eps|a_eps)->a_eps"[[j]] <- "diag"
	}

	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- vector("list", length = p)
	G$"p(nu|Sigma_nu)->sigsq_m" <- vector("list", length = p)
	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- vector("list", length = p)
	G$"p(nu|Sigma_nu)->sigsq_p" <- vector("list", length = p)
	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- vector("list", length = p)
	G$"p(sigsq_m|a_m)->sigsq_m" <- vector("list", length = p)
	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- vector("list", length = p)
	G$"p(sigsq_p|a_p)->sigsq_p" <- vector("list", length = p)
	eta_vec$"p(sigsq_m|a_m)->a_m" <- vector("list", length = p)
	G$"p(sigsq_m|a_m)->a_m" <- vector("list", length = p)
	eta_vec$"p(sigsq_p|a_p)->a_p" <- vector("list", length = p)
	G$"p(sigsq_p|a_p)->a_p" <- vector("list", length = p)
	for(j in 1:p) {

		eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- c(-K/2, -K/2)
		G$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- "full"

		eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- replicate(L, c(-K/2, -K/2))
		G$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- rep("full", L)

		eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- c(-3/2, -1/2)
		G$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- "full"

		eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]] <- replicate(L, c(-3/2, -1/2))
		G$"p(sigsq_p|a_p)->sigsq_p"[[j]] <- rep("full", L)

		eta_vec$"p(sigsq_m|a_m)->a_m"[[j]] <- c(-1/2, -1/2)
		G$"p(sigsq_m|a_m)->a_m"[[j]] <- "diag"

		eta_vec$"p(sigsq_p|a_p)->a_p"[[j]] <- replicate(L, c(-1/2, -1/2))
		G$"p(sigsq_p|a_p)->a_p"[[j]] <- rep("diag", L)
	}

	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
	eta_vec$"p(a_eps)->a_eps" <- vector("list", length = p)
	G$"p(a_eps)->a_eps" <- vector("list", length = p)
	eta_vec$"p(a_m)->a_m" <- vector("list", length = p)
	G$"p(a_m)->a_m" <- vector("list", length = p)
	eta_vec$"p(a_p)->a_p" <- vector("list", length = p)
	G$"p(a_p)->a_p" <- vector("list", length = p)
	for(j in 1:p) {

		eta_vec$"p(a_eps)->a_eps"[[j]] <- igw_prior_updates[[2]]
		G$"p(a_eps)->a_eps"[[j]] <- igw_prior_updates[[1]]

		eta_vec$"p(a_m)->a_m"[[j]] <- igw_prior_updates[[2]]
		G$"p(a_m)->a_m"[[j]] <- igw_prior_updates[[1]]

		eta_vec$"p(a_p)->a_p"[[j]] <- replicate(L, igw_prior_updates[[2]])
		G$"p(a_p)->a_p"[[j]] <- rep(igw_prior_updates[[1]], L)
	}

	elbo_res <- NULL
	converged <- FALSE
	iter <- 0
	while((!converged) & (iter < n_vmp)) {

		iter <- iter + 1

		cat("starting iteration", iter, "of", n_vmp, "\n")

		eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"

		eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"

		eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
		G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

		eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
		G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
		eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
		G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"

		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"

		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"

		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"

		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"

		# Update p(Y|nu,zeta,sigsq_eps) fragment:

		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)

		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)

		mfpc_lik_fragment <- mfpc_lik_frag(eta_in, G_in, C, Y, L)

		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- mfpc_lik_fragment$"eta"[[1]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mfpc_lik_fragment$"eta"[[2]]
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mfpc_lik_fragment$"eta"[[3]]

		G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mfpc_lik_fragment$"G"[[1]]

		# For j = 1, ..., p, update p(nu[j]|Sigma_nu[j]) fragment:

		for(j in 1:p) {

			eta_in <- list(
				eta_vec$"nu->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->nu"[[j]],
				eta_vec$"sigsq_m->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]],
				eta_vec$"sigsq_p->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]]
			)

			G_in <- list(
				G$"sigsq_m->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_m"[[j]],
				G$"sigsq_p->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_p"[[j]]
			)

			fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
				eta_in, G_in, L, mu_beta, Sigma_beta
			)

			eta_vec$"p(nu|Sigma_nu)->nu"[[j]] <- fpc_gauss_pen_fragment$"eta"[[1]]
			eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- fpc_gauss_pen_fragment$"eta"[[2]]
			eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- fpc_gauss_pen_fragment$"eta"[[3]]

			G$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- fpc_gauss_pen_fragment$"G"[[1]]
			G$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- fpc_gauss_pen_fragment$"G"[[2]]
		}

		# For j = 1, ..., p, update p(sigsq_m[j]|a_m[j]) fragment:

		for(j in 1:p) {

			eta_in <- list(
				eta_vec$"sigsq_m->p(sigsq_m|a_m)"[[j]],
				eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]],
				eta_vec$"a_m->p(sigsq_m|a_m)"[[j]],
				eta_vec$"p(sigsq_m|a_m)->a_m"[[j]]
			)

			iter_igw_fragment <- iter_igw_frag(
				eta_in, G$"a_m->p(sigsq_m|a_m)"[[j]],
				1, G$"sigsq_m->p(sigsq_m|a_m)"[[j]]
			)

			eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- iter_igw_fragment$"eta"[[1]]
			eta_vec$"p(sigsq_m|a_m)->a_m"[[j]] <- iter_igw_fragment$"eta"[[2]]

			G$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- iter_igw_fragment$"G"[[1]]
			G$"p(sigsq_m|a_m)->a_m"[[j]] <- iter_igw_fragment$"G"[[2]]
		}

		# For j = 1, ..., p, for l = 1, ..., L, update p(sigsq_p[j][l]|a_p[j][l]) fragment:

		for(j in 1:p) {

			for(l in 1:L) {

				eta_in <- list(
					eta_vec$"sigsq_p->p(sigsq_p|a_p)"[[j]][, l],
					eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l],
					eta_vec$"a_p->p(sigsq_p|a_p)"[[j]][, l],
					eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l]
				)

				iter_igw_fragment <- iter_igw_frag(
					eta_in, G$"a_p->p(sigsq_p|a_p)"[[j]][l],
					1, G$"sigsq_p->p(sigsq_p|a_p)"[[j]][l]
				)

				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l] <- iter_igw_fragment$"eta"[[1]]
				eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l] <- iter_igw_fragment$"eta"[[2]]

				G$"p(sigsq_p|a_p)->sigsq_p"[[j]][l] <- iter_igw_fragment$"G"[[1]]
				G$"p(sigsq_p|a_p)->a_p"[[j]][l] <- iter_igw_fragment$"G"[[2]]
			}
		}

		# Compute the entropy:

		ent <- 0

		for(j in 1:p) {

			eta_in <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"[[j]],
				eta_vec$"p(nu|Sigma_nu)->nu"[[j]]
			)
			ent_nu <- entropy_gauss(eta_in, use_vech = FALSE)

			ent <- ent + ent_nu
		}

		for(i in 1:N) {

			eta_in <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[, i],
				eta_vec$"p(zeta)->zeta"[, i]
			)
			ent_zeta <- entropy_gauss(eta_in, use_vech = TRUE)

			ent <- ent + ent_zeta
		}

		for(j in 1:p) {

			eta_sigsq_eps <- list(
				eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[[j]],
				eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]]
			)
			G_sigsq_eps <- c(
				G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[j],
				G$"p(sigsq_eps|a_eps)->sigsq_eps"[j]
			)
			ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)

			ent <- ent + ent_sigsq_eps

			eta_a_eps <- list(
				eta_vec$"p(sigsq_eps|a_eps)->a_eps"[[j]],
				eta_vec$"p(a_eps)->a_eps"[[j]]
			)
			G_a_eps <- c(
				G$"p(sigsq_eps|a_eps)->a_eps"[[j]],
				G$"p(a_eps)->a_eps"[[j]]
			)
			ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)

			ent <- ent + ent_a_eps

			eta_sigsq_m <- list(
				eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]],
				eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]]
			)
			G_sigsq_m <- c(
				G$"p(nu|Sigma_nu)->sigsq_m"[[j]],
				G$"p(sigsq_m|a_m)->sigsq_m"[[j]]
			)
			ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)

			ent <- ent + ent_sigsq_m

			eta_a_m <- list(
				eta_vec$"p(sigsq_m|a_m)->a_m"[[j]],
				eta_vec$"p(a_m)->a_m"[[j]]
			)
			G_a_m <- c(
				G$"p(sigsq_m|a_m)->a_m"[[j]],
				G$"p(a_m)->a_m"[[j]]
			)
			ent_a_m <- entropy_igw(eta_a_m, G_a_m)

			ent <- ent + ent_a_m

			eta_a_m <- list(
				eta_vec$"p(sigsq_m|a_m)->a_m"[[j]],
				eta_vec$"p(a_m)->a_m"[[j]]
			)
			G_a_m <- c(
				G$"p(sigsq_m|a_m)->a_m"[[j]],
				G$"p(a_m)->a_m"[[j]]
			)
			ent_a_m <- entropy_igw(eta_a_m, G_a_m)

			ent <- ent + ent_a_m

			for(l in 1:L) {

				eta_sigsq_p <- list(
					eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]][, l],
					eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l]
				)
				G_sigsq_p <- list(
					G$"p(nu|Sigma_nu)->sigsq_p"[[j]][l],
					G$"p(sigsq_p|a_p)->sigsq_p"[[j]][l]
				)
				ent_sigsq_p <- entropy_igw(eta_sigsq_p, G_sigsq_p)

				ent <- ent + ent_sigsq_p

				eta_a_p <- list(
					eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l],
					eta_vec$"p(a_p)->a_p"[[j]][,l]
				)
				G_a_p <- c(
					G$"p(sigsq_p|a_p)->a_p"[[j]][l],
					G$"p(a_p)->a_p"[[j]][l]
				)
				ent_a_p <- entropy_igw(eta_a_p, G_a_p)
				ent <- ent + ent_a_p
			}
		}

		# Compute the cross-entropy:

		c_ent <- 0

		eta_in <- list(
			eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
			eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)
		G_in <- list(
			G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
			G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
		)

		c_ent_p_Y <- cross_entropy_mfpc_lik_frag(eta_in, G_in, C, Y, L)

		c_ent <- c_ent + c_ent_p_Y

		for(j in 1:p) {

			eta_in <- list(
				eta_vec$"nu->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->nu"[[j]],
				eta_vec$"sigsq_m->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]],
				eta_vec$"sigsq_p->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]]
			)
			G_in <- list(
				G$"sigsq_m->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_m"[[j]],
				G$"sigsq_p->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_p"[[j]]
			)
			c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)

			c_ent <- c_ent + c_ent_p_nu
		}

		for(i in 1:N) {

			eta_in <- list(
				eta_vec$"zeta->p(zeta)"[, i],
				eta_vec$"p(zeta)->zeta"[, i]
			)
			c_ent_p_zeta <- cross_entropy_gauss_prior(eta_in, rep(0, L), diag(L), use_vech = TRUE)

			c_ent <- c_ent + c_ent_p_zeta
		}

		for(j in 1:p) {

			eta_in <- list(
				eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)"[[j]],
				eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]],
				eta_vec$"a_eps->p(sigsq_eps|a_eps)"[[j]],
				eta_vec$"p(sigsq_eps|a_eps)->a_eps"[[j]]
			)
			G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"[[j]]
			G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"[[j]]
			c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

			c_ent <- c_ent + c_ent_p_sigsq_eps

			eta_in <- list(
				eta_vec$"a_eps->p(a_eps)"[[j]],
				eta_vec$"p(a_eps)->a_eps"[[j]]
			)
			G_in <- c(
				G$"a_eps->p(a_eps)"[[j]],
				G$"p(a_eps)->a_eps"[[j]]
			)
			c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

			c_ent <- c_ent + c_ent_p_a_eps

			eta_in <- list(
				eta_vec$"sigsq_m->p(sigsq_m|a_m)"[[j]],
				eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]],
				eta_vec$"a_m->p(sigsq_m|a_m)"[[j]],
				eta_vec$"p(sigsq_m|a_m)->a_m"[[j]]
			)
			G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"[[j]]
			G_hyper <- G$"a_m->p(sigsq_m|a_m)"[[j]]
			c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

			c_ent <- c_ent + c_ent_p_sigsq_m

			eta_in <- list(
				eta_vec$"a_m->p(a_m)"[[j]],
				eta_vec$"p(a_m)->a_m"[[j]]
			)
			G_in <- c(
				G$"a_m->p(a_m)"[[j]],
				G$"p(a_m)->a_m"[[j]]
			)
			c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

			c_ent <- c_ent + c_ent_p_a_m

			for(l in 1:L) {

				eta_in <- list(
					eta_vec$"sigsq_p->p(sigsq_p|a_p)"[[j]][, l],
					eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l],
					eta_vec$"a_p->p(sigsq_p|a_p)"[[j]][, l],
					eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l]
				)
				G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[[j]][l]
				G_hyper <- G$"a_p->p(sigsq_p|a_p)"[[j]][l]
				c_ent_p_sigsq_p <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

				c_ent <- c_ent + c_ent_p_sigsq_p

				eta_in <- list(
					eta_vec$"a_p->p(a_p)"[[j]][, l],
					eta_vec$"p(a_p)->a_p"[[j]][, l]
				)
				G_in <- c(
					G$"a_p->p(a_p)"[[j]][l],
					G$"p(a_p)->a_p"[[j]][l]
				)
				c_ent_p_a_p <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

				c_ent <- c_ent + c_ent_p_a_p
			}
		}

		# Compute the ELBO

		elbo_new <- ent - c_ent
		elbo_res <- c(elbo_res, elbo_new)

		if(plot_elbo) {

			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab = "iterations", ylab = "ELBO")
		}

		if(iter > 1) {

			elbo_old <- elbo_res[iter - 1]
			criterion_1_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)

			if(iter > 2) {

				elbo_old <- elbo_res[iter - 2]
				criterion_2_satisfied <- (abs(elbo_new/elbo_old - 1) < criterion)
			} else {

				criterion_2_satisfied <- FALSE
			}

			criterion_satisfied <- (criterion_1_satisfied || criterion_2_satisfied)

			if(criterion_satisfied) {

				converged <- TRUE
			}
		}
	}

	# Get the list of natural parameter vectors:

	return(eta_vec)
}

# logistic_fpca <- function(
# 	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
# 	Sigma_beta, A, time_g, C_g, Psi_g,
# 	criterion, n_mc=100, plot_elbo=FALSE
# ) {
#
# 	# Establish necessary parameters:
#
# 	Sigma_zeta <- sigma_zeta^2*diag(L)
# 	T_vec <- sapply(Y, length)
# 	K <- dim(C[[1]])[2] - 2
# 	d <- (K+2)*(L+1)
#
# 	# Initialise VMP simulation:
#
# 	mu_q_zeta <- vector("list", length=N)
# 	Sigma_q_zeta <- vector("list", length=N)
# 	for(i in 1:N) {
#
# 		mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
# 		Sigma_q_zeta[[i]] <- diag(L)
# 	}
#
# 	eta_vec <- vector("list", length=24)
# 	names(eta_vec) <- c(
# 		"nu->p(Y|nu,zeta)", "p(Y|nu,zeta)->nu",
# 		"zeta->p(Y|nu,zeta)", "p(Y|nu,zeta)->zeta",
# 		"zeta->p(zeta)", "p(zeta)->zeta",
# 		"nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
# 		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
# 		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
# 		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
# 		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
# 		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
# 		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
# 		"a_m->p(a_m)", "p(a_m)->a_m",
# 		"a_p->p(a_p)", "p(a_p)->a_p"
# 	)
#
# 	G <- vector("list", length=16)
# 	names(G) <- c(
# 		"sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
# 		"sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
# 		"sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
# 		"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
# 		"a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
# 		"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
# 		"a_m->p(a_m)", "p(a_m)->a_m",
# 		"a_p->p(a_p)", "p(a_p)->a_p"
# 	)
#
# 	eta_1_sum <- 0
# 	eta_2_sum <- 0
# 	xi <- vector("list", length=N)
# 	for(i in 1:N) {
#
# 		xi[[i]] <- rep(1, T_vec[i])
# 		A_xi <- -tanh(xi[[i]]/2)/(4*xi[[i]])
#
# 		mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
# 		Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
# 		M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
#
# 		sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]] - 0.5)
# 		eta_1_sum <- eta_1_sum + sum_val
#
# 		M <- crossprod(C[[i]], diag(A_xi)%*%C[[i]])
# 		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, M))
# 		eta_2_sum <- eta_2_sum + sum_val
# 	}
# 	eta_1 <- eta_1_sum
# 	eta_2 <- eta_2_sum
# 	eta_vec$"p(Y|nu,zeta)->nu" <- c(eta_1, eta_2)
#
# 	D_L <- duplication.matrix(L)
# 	eta_1 <- Reduce(cbind, mu_q_zeta)
# 	eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
# 	eta_vec$"p(Y|nu,zeta)->zeta" <- rbind(eta_1, eta_2)
#
# 	eta_vec$"p(zeta)->zeta" <- replicate(
# 		N,
# 		gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
# 	)
#
# 	eta_1 <- rep(0, d)
# 	eta_2 <- -0.5*as.vector(diag(d))
# 	eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)
#
# 	eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
# 	G$"p(nu|Sigma_nu)->sigsq_m" <- "full"
#
# 	eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
# 	G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)
#
# 	eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
# 	G$"p(sigsq_m|a_m)->sigsq_m" <- "full"
#
# 	eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
# 	G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)
#
# 	eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
# 	G$"p(sigsq_m|a_m)->a_m" <- "diag"
#
# 	eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
# 	G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)
#
# 	igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
#
# 	eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
# 	G$"p(a_m)->a_m" <- igw_prior_updates[[1]]
#
# 	eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
# 	G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)
#
# 	elbo_res <- NULL
# 	elbo_new <- -Inf
# 	converged <- FALSE
# 	iter <- 0
#
# 	while((!converged) & (iter < n_vmp)) {
#
# 		elbo_old <- elbo_new
# 		iter <- iter + 1
#
# 		if(plot_elbo) {
#
# 			cat("starting iteration", iter, "of", n_vmp, "\n")
# 		}
#
# 		eta_vec$"nu->p(Y|nu,zeta)" <- eta_vec$"p(nu|Sigma_nu)->nu"
# 		eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta)->nu"
#
# 		eta_vec$"zeta->p(Y|nu,zeta)" <- eta_vec$"p(zeta)->zeta"
# 		eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta)->zeta"
#
# 		eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
# 		G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
# 		eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
# 		G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"
#
# 		eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
# 		G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
# 		eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
# 		G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"
#
# 		eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
# 		G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
# 		eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
# 		G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"
#
# 		eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
# 		G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
# 		eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
# 		G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
#
# 		# Update p(Y|nu,zeta) fragment:
#
# 		eta_in <- list(
# 			eta_vec$"nu->p(Y|nu,zeta)",
# 			eta_vec$"p(Y|nu,zeta)->nu",
# 			eta_vec$"zeta->p(Y|nu,zeta)",
# 			eta_vec$"p(Y|nu,zeta)->zeta"
# 		)
# 		logistic_fpc_lik_fragment <- logistic_fpc_lik_frag(eta_in, C, Y, L)
#
# 		eta_vec$"p(Y|nu,zeta)->nu" <- logistic_fpc_lik_fragment[[1]]
# 		eta_vec$"p(Y|nu,zeta)->zeta" <- logistic_fpc_lik_fragment[[2]]
#
# 		# Update p(nu|Sigma_nu) fragment:
#
# 		eta_in <- list(
# 			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
# 			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
# 			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
# 		)
#
# 		G_in <- list(
# 			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
# 			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
# 		)
#
# 		fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
# 			eta_in, G_in, L, mu_beta, Sigma_beta
# 		)
#
# 		eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
# 		eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
# 		eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]
#
# 		G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
# 		G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]
#
# 		# Update p(sigsq_m|a_m) fragment:
#
# 		eta_in <- list(
# 			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
# 			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
# 			eta_vec$"a_m->p(sigsq_m|a_m)",
# 			eta_vec$"p(sigsq_m|a_m)->a_m"
# 		)
#
# 		iter_igw_fragment <- iter_igw_frag(
# 			eta_in, G$"a_m->p(sigsq_m|a_m)",
# 			1, G$"sigsq_m->p(sigsq_m|a_m)"
# 		)
#
# 		eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
# 		eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]
#
# 		G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
# 		G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]
#
# 		# Update p(sigsq_p|a_p) fragment:
#
# 		for(l in 1:L) {
#
# 			eta_in <- list(
# 				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
# 				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
# 				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
# 				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
# 			)
#
# 			iter_igw_fragment <- iter_igw_frag(
# 				eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
# 				1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
# 			)
#
# 			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
# 			eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
#
# 			G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
# 			G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
# 		}
#
# 		# Compute the entropy:
#
# 		ent <- 0
#
# 		eta_nu <- list(
# 			eta_vec$"p(Y|nu,zeta)->nu",
# 			eta_vec$"p(nu|Sigma_nu)->nu"
# 		)
# 		ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)
#
# 		ent <- ent + ent_nu
#
# 		ent_zeta <- 0
# 		for(i in 1:N) {
#
# 			eta_zeta <- list(
# 				eta_vec$"p(Y|nu,zeta)->zeta"[,l],
# 				eta_vec$"p(zeta)->zeta"[,l]
# 			)
# 			sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
# 			ent_zeta <- ent_zeta + sum_val
# 		}
#
# 		ent <- ent + ent_zeta
#
# 		eta_sigsq_m <- list(
# 			eta_vec$"p(nu|Sigma_nu)->sigsq_m",
# 			eta_vec$"p(sigsq_m|a_m)->sigsq_m"
# 		)
# 		G_sigsq_m <- c(
# 			G$"p(nu|Sigma_nu)->sigsq_m",
# 			G$"p(sigsq_m|a_m)->sigsq_m"
# 		)
# 		ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)
#
# 		ent <- ent + ent_sigsq_m
#
# 		eta_a_m <- list(
# 			eta_vec$"p(sigsq_m|a_m)->a_m",
# 			eta_vec$"p(a_m)->a_m"
# 		)
# 		G_a_m <- c(
# 			G$"p(sigsq_m|a_m)->a_m",
# 			G$"p(a_m)->a_m"
# 		)
# 		ent_a_m <- entropy_igw(eta_a_m, G_a_m)
#
# 		ent <- ent + ent_a_m
#
# 		ent_sigsq_p <- 0
# 		ent_a_p <- 0
# 		for(l in 1:L) {
#
# 			eta_sigsq_p <- list(
# 				eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
# 				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
# 			)
# 			G_sigsq_p <- list(
# 				G$"p(nu|Sigma_nu)->sigsq_p"[l],
# 				G$"p(sigsq_p|a_p)->sigsq_p"[l]
# 			)
# 			sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
# 			ent_sigsq_p <- ent_sigsq_p + sum_val
#
# 			eta_a_p <- list(
# 				eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
# 				eta_vec$"p(a_p)->a_p"[,l]
# 			)
# 			G_a_p <- c(
# 				G$"p(sigsq_p|a_p)->a_p"[l],
# 				G$"p(a_p)->a_p"[l]
# 			)
# 			sum_val <- entropy_igw(eta_a_p, G_a_p)
# 			ent_a_p <- ent_a_p + sum_val
# 		}
#
# 		ent <- ent + ent_sigsq_p + ent_a_p
#
# 		# Compute the cross-entropy:
#
# 		c_ent <- 0
#
# 		eta_in <- list(
# 			eta_vec$"nu->p(Y|nu,zeta)",
# 			eta_vec$"p(Y|nu,zeta)->nu",
# 			eta_vec$"zeta->p(Y|nu,zeta)",
# 			eta_vec$"p(Y|nu,zeta)->zeta"
# 		)
#
# 		c_ent_p_Y <- cross_entropy_logistic_fpc_lik_frag(eta_in, C, Y, L)
#
# 		c_ent <- c_ent + c_ent_p_Y
#
# 		eta_in <- list(
# 			eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
# 			eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
# 			eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
# 		)
# 		G_in <- list(
# 			G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
# 			G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
# 		)
# 		c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)
#
# 		c_ent <- c_ent + c_ent_p_nu
#
# 		c_ent_p_zeta <- 0
# 		for(i in 1:N) {
#
# 			eta_in <- list(
# 				eta_vec$"zeta->p(zeta)"[,l],
# 				eta_vec$"p(zeta)->zeta"[,l]
# 			)
# 			sum_val <- cross_entropy_gauss_prior(
# 				eta_in, rep(0, L),
# 				Sigma_zeta, use_vech=TRUE
# 			)
# 			c_ent_p_zeta <- c_ent_p_zeta + sum_val
# 		}
#
# 		c_ent <- c_ent + c_ent_p_zeta
#
# 		eta_in <- list(
# 			eta_vec$"sigsq_m->p(sigsq_m|a_m)",
# 			eta_vec$"p(sigsq_m|a_m)->sigsq_m",
# 			eta_vec$"a_m->p(sigsq_m|a_m)",
# 			eta_vec$"p(sigsq_m|a_m)->a_m"
# 		)
# 		G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
# 		G_hyper <- G$"a_m->p(sigsq_m|a_m)"
# 		c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
#
# 		c_ent <- c_ent + c_ent_p_sigsq_m
#
# 		eta_in <- list(
# 			eta_vec$"a_m->p(a_m)",
# 			eta_vec$"p(a_m)->a_m"
# 		)
# 		G_in <- c(
# 			G$"a_m->p(a_m)",
# 			G$"p(a_m)->a_m"
# 		)
# 		c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
#
# 		c_ent <- c_ent + c_ent_p_a_m
#
# 		c_ent_p_sigsq_p <- 0
# 		c_ent_p_a_p <- 0
# 		for(l in 1:L) {
#
# 			eta_in <- list(
# 				eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
# 				eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
# 				eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
# 				eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
# 			)
# 			G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
# 			G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
# 			sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
# 			c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val
#
# 			eta_in <- list(
# 				eta_vec$"a_p->p(a_p)"[,l],
# 				eta_vec$"p(a_p)->a_p"[,l]
# 			)
# 			G_in <- c(
# 				G$"a_p->p(a_p)"[l],
# 				G$"p(a_p)->a_p"[l]
# 			)
# 			sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
# 			c_ent_p_a_p <- c_ent_p_a_p + sum_val
# 		}
#
# 		c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p
#
# 		# Compute the ELBO
#
# 		elbo_new <- ent - c_ent
# 		elbo_res <- c(elbo_res, elbo_new)
#
# 		if(plot_elbo) {
#
# 			plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
# 		}
#
#
# 		if(abs(elbo_new/elbo_old - 1) < criterion) {
#
# 			converged <- TRUE
# 		}
# 	}
#
# 	# Save the original q_nu:
#
# 	eta_nu <- list(eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta)->nu")
#
# 	q_nu <- gauss_q(eta_nu, use_vech = FALSE)
# 	mu_q_nu <- q_nu[[1]]
# 	Sigma_q_nu <- q_nu[[2]]
#
# 	# Set up the orthogonal decomposition:
#
# 	eta_in <- list(
# 		eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta)->nu",
# 		eta_vec$"p(zeta)->zeta", eta_vec$"p(Y|nu,zeta)->zeta"
# 	)
#
# 	fpc_rotns <- fpc_rotation(eta_in, time_g, C_g, Psi_g)
#
# 	mu_q_nu_mu <- fpc_rotns$"mu"[[1]]
# 	Sigma_q_nu_mu <- fpc_rotns$"mu"[[2]]
# 	mu_q_mu <- fpc_rotns$"mu"[[3]]
#
# 	mu_q_nu_psi <- fpc_rotns$"psi"[[1]]
# 	Sigma_q_nu_psi <- fpc_rotns$"psi"[[2]]
# 	M_q_Psi_star <- fpc_rotns$"psi"[[3]]
#
# 	mu_q_zeta <- fpc_rotns$"zeta"[[1]]
# 	Sigma_q_zeta <- fpc_rotns$"zeta"[[2]]
# 	M_q_Zeta_star <- fpc_rotns$"zeta"[[3]]
#
# 	res <- list(mu_q_mu, M_q_Psi_star)
# 	names(res) <- c("mu", "Psi")
# 	return(res)
# }
