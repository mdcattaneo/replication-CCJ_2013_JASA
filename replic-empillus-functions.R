##################################################################################
## R CODE FOR CATTANEO-CRUMP-JANSSON (2012, JASA): Generalized Jackknife for WAD
## DATE: 24-Jan-2013
## NOTE: This functions are used in the main file (replic-empillus.R)
## FUNCTIONS:
##     1) ROT Bandwidth Selectors
##     2) C Code Wrapper
##     3) GenJack Bias-Correction
##################################################################################

##############################################################
## FUNCTION 1: ROT Bandwidth Selectors
##############################################################
scale = function(x) {q = quantile(x, c(.25,.75)); return(min(sd(x), (q[2]-q[1])/1.349))}

ROT.h.1d = function(s, CB, CS, CH, P, n) {
    d = length(s)
    if (sign(CB) == sign(CS*CH)) {
        h.opt = (s[1]^P * prod(s) * d * abs(CB) / (P * abs(CS*CH) * n))^(1/(P+d))
    } else {
        h.opt = (s[1]^P * prod(s) * abs(CB) / (abs(CS*CH) * n))^(1/(P+d))
    }
    return(h.opt)
}
## TEST: ROT.h.1d(s=c(1,1,1), CB = -0.160354, CS = 1/8, CH = -3/(8*pi^(3/2)), P=4, n=700)

ROT.h.tr = function(s, b, CB, CS, CH, P, n) {
    d = length(s); b.P = crossprod(b^2,s^(-P)); b.2P = crossprod(b^2,s^(-2*P))

    h.opt = (prod(s)*abs(CB)*b.P*(sqrt((P-d)^2+4*P*d*crossprod(b)*b.2P/b.P^2) - sign(CS*CH/CB)*(P-d)) 
                / (abs(CS*CH)*2*P*b.2P*n))^(1/(P+d))

    return(h.opt)
}
## TEST: ROT.h.tr(s=c(1,1,1), b=c(1,2,3), CB = -0.160354, CS = 1/8, CH = -3/(8*pi^(3/2)), P=4, n=700)

ROT.hs.tr = function(s, b, CB, CS, CH, P, n) {
    d = 3
    AMSE = function(h) crossprod(b^2, (CB/(n*h[1]*h[2]*h[3]) + CS*CH*c(h[1]/s[1]^P,h[2]/s[2]^P,h[3]/s[3]^P)/prod(s))^2)
    DAMSE = function(h) c(-(1/(h[1]^3*h[2]^2*h[3]^2*n^2))*2*s[1]^(-2*(1+P))*s[2]^(-2-P)*s[3]^(-2-P)
                              *(b[3]^2*CB*s[1]^(1+2*P)*s[2]^(1+P)*s[3]*(CH*CS*h[1]*h[2]*h[3]^(1+P)*n
                                  +CB*s[1]*s[2]*s[3]^(1+P))+s[3]^P*(-b[1]^2*s[2]^P*(CH*CS*h[1]^(1+P)*h[2]*h[3]*n*P
                                      -CB*s[1]^(1+P)*s[2]*s[3])*(CH*CS*h[1]^(1+P)*h[2]*h[3]*n+CB*s[1]^(1+P)*s[2]*s[3])
                                      +b[2]^2*CB*s[1]^(1+2*P)*s[2]*s[3]*(CH*CS*h[1]*h[2]^(1+P)*h[3]*n+CB*s[1]*s[2]^(1+P)*s[3]))),

                          -(1/(h[1]^2*h[2]^3*h[3]^2*n^2))*2*s[1]^(-2-P)*s[2]^(-2*(1+P))*s[3]^(-2-P)
                              *(b[3]^2*CB*s[1]^(1+P)*s[2]^(1+2*P)*s[3]*(CH*CS*h[1]*h[2]*h[3]^(1+P)*n+CB*s[1]*s[2]*s[3]^(1+P))
                                  +s[3]^P*(b[1]^2*CB*s[1]*s[2]^(1+2*P)*s[3]*(CH*CS*h[1]^(1+P)*h[2]*h[3]*n
                                      +CB*s[1]^(1+P)*s[2]*s[3])-b[2]^2*s[1]^P*(CH*CS*h[1]*h[2]^(1+P)*h[3]*n*P
                                      -CB*s[1]*s[2]^(1+P)*s[3])*(CH*CS*h[1]*h[2]^(1+P)*h[3]*n+CB*s[1]*s[2]^(1+P)*s[3]))),

                          -(1/(h[1]^2*h[2]^2*h[3]^3*n^2))*2*s[1]^(-2-P)*s[2]^(-2-P)*s[3]^(-2*(1+P))
                              *(-b[3]^2*s[1]^P*s[2]^P*(CH*CS*h[1]*h[2]*h[3]^(1+P)*n*P
                                  -CB*s[1]*s[2]*s[3]^(1+P))*(CH*CS*h[1]*h[2]*h[3]^(1+P)*n+CB*s[1]*s[2]*s[3]^(1+P))
                                  +CB*s[1]*s[2]*s[3]^(1+2*P)*(b[1]^2*s[2]^P*(CH*CS*h[1]^(1+P)*h[2]*h[3]*n
                                      +CB*s[1]^(1+P)*s[1]*s[3])+b[2]^2*s[1]^P*(CH*CS*h[1]*h[2]^(1+P)*h[3]*n+CB*s[1]*s[2]^(1+P)*s[3])))
                          )

    h.opt = optim(par=rep(ROT.h.1d(s, CB, CS, CH, P, n),3), fn=AMSE, gr=DAMSE,
                  method="L-BFGS-B", lower=c(.3,.3,.3), upper=c(.9,.9,.9))

    return(h.opt$par)
}
## TEST: ROT.hs.tr(s=c(1,2/3,1), b=c(1,1,1), CB = -0.160354, CS = 1/8, CH = -3/(8*pi^(3/2)), P=4, n=700)

##############################################################
## FUNCTION 2: C Code Wrapper
##############################################################
C.wad = function(y,x,d,trim,P,n,h) {
    tmp11=rep.int(0,1);
    tmp21=rep.int(0,1);
    tmp51=rep.int(0,n*(n-1)/2); tmp52=rep.int(0,n*(n-1)/2); tmp61=rep.int(0,n); tmp62=rep.int(0,n); tmp71=rep.int(1,n)
    tmp81=rep.int(0,n); tmp82=rep.int(0,n);
    tmp91=rep.int(0,n); tmp92=rep.int(0,n);

    out.C = .C("wad",theta1.hat.LI=as.double(tmp11),
                     sigma1.hat.LI=as.double(tmp21),
                     K.ij=as.double(tmp51),DK.ij=as.double(tmp52),w.i=as.double(tmp61),Dw.i=as.double(tmp62),used=as.integer(tmp71),
                     f_LI_i=as.double(tmp81),Df_LI_i=as.double(tmp82),
                     e_LI_i=as.double(tmp91),De_LI_i=as.double(tmp92),
                     y=as.double(y),x=as.double(x),d=as.integer(d),trim=as.double(trim),
                     P=as.integer(P),n=as.integer(n),h=as.double(h))

    return(list(theta1.hat=out.C$theta1.hat.LI, sigma1.hat=out.C$sigma1.hat.LI, used=out.C$used))
}

##############################################################
## FUNCTION 3: GenJack Bias-Correction
##############################################################
genjack = function(data, d=3, nc.down, nc.up, thetalab) {
    data = cbind(data,NA); colnames(data)=c(colnames(data[,1:(ncol(data)-1)]),"theta1.hat-jack")

    hs = as.matrix(data[data[,"s"]==1,"c"]); h.len = nrow(hs)
    cs = NULL
    if (nc.down>0) for (i in 1:nc.down) cs = cbind(cs, hs[i:(h.len-nc.up-nc.down-1+i),]/hs[(1+nc.down):(h.len-nc.up),])
    if (nc.up>0) for (i in 1:nc.up) cs = cbind(cs, hs[(nc.down+1+i):(h.len-nc.up+i),]/hs[(1+nc.down):(h.len-nc.up),])
    
    cs.lambdas = cbind(1, cs, matrix(NA,nrow=nrow(cs),ncol=ncol(cs)+1))

    lambdas = function(cs,d) {
        cs = cs[!is.na(cs)]; J = length(cs)-1; M = rep(1,J+1)
        for (j in 0:(J-1)) M = rbind(M, cs^(2*j-d))
        M.inv = solve(M)
        return(M.inv[,1])
    }

    cs.lambdas[,(ncol(cs)+2):ncol(cs.lambdas)] = t(apply(cs.lambdas, 1, function(x) lambdas(x,d=d)))

    theta.hat.jack = function(x, nc.down, nc.up, thetalab, cs.lambdas) {
        theta.hats = x[(1+nc.down):(h.len-nc.up), thetalab]
        if (nc.down>0) for (i in 1:nc.down) theta.hats = cbind(theta.hats, x[i:(h.len-nc.up-nc.down-1+i),thetalab])
        if (nc.up>0) for (i in 1:nc.up) theta.hats = cbind(theta.hats, x[(nc.down+1+i):(h.len-nc.up+i),thetalab])

        out = rep(NA,nrow(x))
        out[(1+nc.down):(h.len-nc.up)] = rowSums(theta.hats*cs.lambdas[,(ncol(cs)+2):ncol(cs.lambdas)])
        return(out)
    }

    for (s in 1:max(data[,"s"])) {
        data[data[,"s"]==s, "theta1.hat-jack"] = theta.hat.jack(x=data[data[,"s"]==s,], nc.down=nc.down, nc.up=nc.up,
                                                                thetalab=thetalab, cs.lambdas=cs.lambdas)
    }

    return(data)
}




