#-------------------------------------------------------------------------------
#-File  : gmcoda.R
#-Aim   : Graphical model for multiple compositional vectors
#-Author: Huaying Fang (Capital Normal University)
#-Email : hyfang@cnu.edu.cn
#-Date  : 18OCT2023
#-Import: R packages glasso, huge and SpiecEasi
#-------------------------------------------------------------------------------
gmcoda_wrap = function(mat1, mat2, propFilt = 0.8, isCnt = T, only.gmcoda = F, nlambda = 40, lambda.min.ratio = 1E-2) {
	matOTU = list(mat1, mat2);
	#-Filter OTU count matrix: Keep OTUs/samples with more than 80% non-0s across samples/OTUs
	matOTUFilt = lapply(matOTU, function(xxx) {
		xxa = xxx[, colMeans(xxx > 0) >= propFilt];
		xxa[rowMeans(xxa > 0) >= propFilt, ];
	});
	xxa = intersect(rownames(matOTUFilt[[1]]), rownames(matOTUFilt[[2]]));
	data_list = lapply(1:2, function(xxk) matOTUFilt[[xxk]][xxa, ]);
	p0s = sapply(data_list, ncol);
	#-Partial correlation matrix from precision matrix
	omega2pcorr = function(mat) {
		aa = mat;
		bb = - aa / sqrt(diag(aa)) / rep(sqrt(diag(aa)), each = nrow(aa));
		diag(bb) = 1;
		return(bb);
	};
	#-gmcoda
	fit_gm = gmcoda(data_list, isCnt = isCnt, lamb_min_ratio = lambda.min.ratio, nlamb = nlambda);
	out = list(gmcoda = omega2pcorr(fit_gm$opt_icov));
	if(!only.gmcoda) {
		#-glasso
		set.seed(0);
		xxa = data_list;
		if(isCnt) xxa = lapply(data_list, function(xx_) return((xx_ + 0.5)/rowSums(xx_ + 0.5)));
		xxa = log(do.call("cbind", xxa));
		xxb = huge::huge(xxa, method = "glasso", nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, verbose = F);
		fit_gl = huge::huge.select(xxb, criterion = "stars", stars.thresh = 0.05, verbose = F);
		out$glasso = omega2pcorr(fit_gl$opt.icov);
		set.seed(NULL);
		#-SEGL (SpiecEasi)
		fit_segl = spiec.easi(data_list, method = "glasso", nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, pulsar.params = list(thresh = 0.05), verbose = F);
		out$segl = omega2pcorr(as.matrix(getOptiCov(fit_segl)));
		#-SEMB (SpiecEasi)
		fit_semb = spiec.easi(data_list, method = "mb", nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, pulsar.params = list(thresh = 0.05), verbose = F);
		out$semb = getRefit(fit_semb);
	};
	#-Output with edges as each row
	xxa = unlist(sapply(data_list, colnames));
	xxa2 = setNames(rep(1:length(data_list), sapply(data_list, ncol)), xxa);
	xxb = data.frame(OTU1 = xxa, OTU2 = rep(xxa, each = length(xxa)));
	xxa1 = lower.tri(out$gmcoda);
	xxb = xxb[xxa1, ];
	xxb$gmcoda = out$gmcoda[xxa1];
	if(!only.gmcoda) {
		xxb$glasso = out$glasso[xxa1];
		xxb$segl = out$segl[xxa1];
		xxb$semb = out$semb[xxa1];		
	};
	xxa1 = rowSums(abs(xxb[, -(1:2), drop = F])) > 0;
	xxb = xxb[xxa1, ];
	xxb1 = t(apply(cbind(xxa2[xxb$OTU1], xxa2[xxb$OTU2]), 1, sort));
	xxb = cbind(type = paste0(xxb1[,1], "o", xxb1[,2]), xxb);
	return(xxb);
};
#-------------------------------------------------------------------------------
gmcoda = function(data, isCnt = F, pseudoCnt = 0.5, lamb_min_ratio = 1E-2, nlamb = 15, wgt_lamb = 1) {
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
	#-Consider weights varying elements of precision matrix
	if(wgt_lamb == 1) {
		wgt_lamb = matrix(1, nrow = p, ncol = p);
	} else {
		wgt_lamb = wgt_lamb + t(wgt_lamb);
		wgt_lamb = wgt_lamb / max(wgt_lamb);
	};
	#-Generate lambda via lamb_min_ratio and nlamb
	lamb_max = max(abs(S[lower.tri(S)]));
	lamb_min = lamb_min_ratio * lamb_max;
	#-Store fit result
	lambs = exp(seq(log(lamb_max), log(lamb_min), length = nlamb));
	fit = list(lambda = lambs, nloglik = rep(NA, nlamb), icov = list(), path = list(), df = rep(NA, nlamb));
	icov_ = diag(p);
	for(k in 1:nlamb) {
		fit_ = gmcoda_sub(A = S, p0s = p0s, iSig = icov_, lambda = fit$lambda[k], wgt_lamb = wgt_lamb);
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
	return(fit);
};
#----------------------------------------
#-Optimization with given lambda in gmcoda
gmcoda_sub = function(A, p0s, iSig = NULL, lambda = 0.1, tol_err = 1E-4, k_max = 100, wgt_lamb) {
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
		iSig2 = huge_glasso_mod(S = (A2 + t(A2))/2, lambda = lambda, wgt_lamb = wgt_lamb);
		mat0 = invRtOmegaR(iSig = iSig2, p0s = p0s);
		fval_new = obj_gmcoda(iSig = iSig2, A = A, mat0 = mat0, p0s = p0s, lambda = lambda, wgt_lamb = wgt_lamb);
		err1 = max(abs(iSig2 - iSig) / (abs(iSig2) + 1));
		err2 = abs(fval_cur - fval_new) / (abs(fval_new) + 1);
		err = min(err1, err2);
		k = k + 1;
		iSig = iSig2;
		fval_cur = fval_new;
	};
	nloglik = fval_cur - lambda * (sum(abs(iSig * wgt_lamb)) - sum(diag(iSig * wgt_lamb)));
	return(list(iSig = iSig, nloglik = nloglik, n_iter = k, rel_err = err));
};
#----------------------------------------
#-(R^T Omega R)^{-1}
invRtOmegaR = function(iSig, p0s) {
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
	return((mat0 + t(mat0))/2);
};
#----------------------------------------
#-Objective function value of gmcoda (negative log likelihood + penalty)
obj_gmcoda = function(iSig, A, mat0, p0s, lambda, wgt_lamb) {
	d0 = length(p0s);
	id1 = cbind(c(0, cumsum(p0s)[-d0]) + 1, cumsum(p0s));
	mat1 = sapply(1:d0, function(k0) colSums(iSig[id1[k0, 1]:id1[k0, 2], ]));
	nloglik = - determinant(iSig)$modulus - determinant(mat0)$modulus + sum(iSig * A) - sum(crossprod(mat1, A %*% mat1) * mat0);
	pen = lambda * (sum(abs(iSig * wgt_lamb)) - sum(diag(iSig * wgt_lamb)));
	return(c(nloglik) + pen);
};
#----------------------------------------
#-Modified huge::huge.glasso for quick preparation
#-Input S must be covariance matrix
require(huge);
require(glasso);
huge_glasso_mod = function(S, lambda, wgt_lamb) {
	if(all(wgt_lamb == 1)) {
		icov = diag(1/(diag(S) + lambda));
		z = which(rowSums(abs(S) > lambda) > 1);
		q = length(z);
		if (q > 0) {
			out.glasso = .Call("_huge_hugeglasso", S[z, z], lambda, F, F, F);
			icov[z, z] = out.glasso$icov[[1]]; 
   		};
	} else {
		icov = glasso(s = S, rho = lambda * wgt_lamb)$wi;
	};
   return(icov);
};
#-------------------------------------------------------------------------------