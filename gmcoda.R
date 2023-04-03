#-------------------------------------------------------------------------------
# File  : gmcoda.R
# Aim   : Graphical model for multiple compositional vectors
# Author: Fang Huaying (Capital Normal University)
# Email : hyfang@cnu.edu.cn
# Date  : 10OCT2022
# Import: huge::huge.glasso
#-------------------------------------------------------------------------------
gmcoda <- function(data, isCnt = F, pseudoCnt = 0.5, lamb_min_ratio = 1E-2, 
	nlamb = 15) {
	n = nrow(data[[1]]);
	p0s = sapply(data, ncol);
	p = sum(p0s);
	#-Counts => Compositions
	if(isCnt) {
		data = lapply(data, function(xx) {
			(xx + pseudoCnt) / rowSums(xx + pseudoCnt);
		});
	};
	#-Central log ratio (CLR) transformation for compositional data
	S = var(do.call("cbind", lapply(data, function(xx) {
		log(xx) - rowMeans(log(xx));
	})));
	#-Generate lambda via lamb_min_ratio and nlamb
	lamb_max = max(abs(S[lower.tri(S)]));
	lamb_min = lamb_min_ratio * lamb_max;
	#-Store fit result
	lambs = exp(seq(log(lamb_max), log(lamb_min), length = nlamb));
	fit = list(lambda = lambs, nloglik = rep(NA, nlamb), icov = list(), 
		path = list(), df = rep(NA, nlamb));
	icov_ = diag(p);
	for(k in 1:nlamb) {
		fit_ = gmcoda_sub(A = S, p0s = p0s, iSig = icov_, 
			lambda = fit$lambda[k]);
		fit$nloglik[k] = fit_$nloglik;
		fit$icov[[k]] = fit_$iSig;
		fit$path[[k]] = (abs(fit_$iSig) >= 1E-6) + 0;
		diag(fit$path[[k]]) = 0;
		fit$df[k] = sum(fit$path[[k]])/2;
	};
	#-Compute BIC score for lambda selection
	fit$bic_score = fit$nloglik + log(n)/n * fit$df;
	fit$opt_index = which.min(fit$bic_score);
	fit$opt_path = fit$path[[fit$opt_index]];
	fit$opt_icov = fit$icov[[fit$opt_index]];
	fit$opt_lambda = fit$lambda[fit$opt_index];
	fit;
};
#----------------------------------------
#-Optimization with given lambda in gmcoda
gmcoda_sub <- function(A, p0s, iSig = NULL, lambda = 0.1, tol_err = 1E-4, 
	k_max = 100) {
	p = ncol(A);
	if(is.null(iSig)) iSig = diag(p);
	id0 = rep(1:length(p0s), times = p0s);
	mat0 = invRtOmegaR(iSig = iSig, p0s = p0s);
	err = 1;
	k = 0;
	fval_cur = Inf;
	while(err > tol_err && k < k_max) {
		mat1 = mat0[, id0] %*% iSig;
		mat2 = mat1 %*% A;
		mat3 = tcrossprod(mat1, mat2) + mat0;
		A2 = A - mat2[id0, ] - t(mat2)[, id0] + mat3[id0, id0];
		iSig2 = huge_glasso_mod(S = (A2 + t(A2))/2, lambda = lambda);
		mat0 = invRtOmegaR(iSig = iSig2, p0s = p0s);
		fval_new = obj_gmcoda(iSig = iSig2, A = A, mat0 = mat0, p0s = p0s, 
			lambda = lambda);
		err1 = max(abs(iSig2 - iSig) / (abs(iSig2) + 1));
		err2 = abs(fval_cur - fval_new) / (abs(fval_new) + 1);
		err = min(err1, err2);
		k = k + 1;
		iSig = iSig2;
		fval_cur = fval_new;
	};
	nloglik = fval_cur - lambda * (sum(abs(iSig)) - sum(diag(iSig)));
	list(iSig = iSig, nloglik = nloglik, n_iter = k, rel_err = err);
};
#----------------------------------------
#-(R^T Omega R)^{-1}
invRtOmegaR <- function(iSig, p0s) {
	d0 = length(p0s);
	id1 = cbind(c(0, cumsum(p0s)[-d0]) + 1, cumsum(p0s));
	mat0 = matrix(NA, nrow = d0, ncol = d0);
	for(i0 in 1:d0) {
		idx1 = id1[i0, 1]:id1[i0, 2];
		for(j0 in i0:d0) {
			idx2 = id1[j0, 1]:id1[j0, 2];
			mat0[i0, j0] = mat0[j0, i0] = sum(iSig[idx1, idx2]);
		};
	};
	mat0 = solve(mat0);
	(mat0 + t(mat0))/2;
};
#----------------------------------------
#-Objective function value of gmcoda (negative log likelihood + penalty)
obj_gmcoda <- function(iSig, A, mat0, p0s, lambda) {
	d0 = length(p0s);
	id1 = cbind(c(0, cumsum(p0s)[-d0]) + 1, cumsum(p0s));
	mat1 = sapply(1:d0, function(k0) colSums(iSig[id1[k0, 1]:id1[k0, 2], ]));
	nloglik = - determinant(iSig)$modulus - determinant(mat0)$modulus + 
		sum(iSig * A) - sum(crossprod(mat1, A %*% mat1) * mat0);
	pen = lambda * (sum(abs(iSig)) - sum(diag(iSig)));
	c(nloglik) + pen;
};
#----------------------------------------
# Modified huge::huge.glasso for quick preparation
# Input S must be covariance matrix
require(huge);
huge_glasso_mod <- function(S, lambda) {
	icov = diag(1/(diag(S) + lambda));
	z = which(rowSums(abs(S) > lambda) > 1);
	q = length(z);
	if (q > 0) {
		out.glasso = .Call("_huge_hugeglasso", S[z, z], lambda, F, F, F);
		icov[z, z] = out.glasso$icov[[1]]; 
   };
   icov;
};
#-------------------------------------------------------------------------------