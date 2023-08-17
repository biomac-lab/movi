# #' MVSE: A package for estimating (M)osquito-(b)orne (V)iral (S)uitability
# #'
# #' The MVSE package provides computational methods to estimate an index of
# #' 'transmission potential' for mosquito-borne viruses. MVSE allows for
# #' parameterization of particular viruses, host and mosquito-species. The methods
# #' offered are climate-driven and can therefore be location specific.
# #'
# #' For more information, please refer to the manual.
# #'
# #' @docType package
# #' @name MVSE
# NULL

# devtools::use_package("pbapply") # Defaults to imports
# devtools::use_package("genlasso") # Defaults to imports
# devtools::use_package("scales") # Defaults to imports

###########################################
#### CSV EXPORT FUNCTIONS
###########################################


#' Exports the estimated generation times.
#'
#' @param sep CSV field separator character.
#'
#' @return NULL
#'
exportGenerationTimes<- function(sep=','){
  filename<- paste0(.MVSEbuildOutPath(),'.genTimes.csv')
  out<- MVSE_results$Tg
  write.table(out,filename,col.names=TRUE,row.names=FALSE,sep=sep)
}

#' Exports the estimated entomological parameters of mosquito life-span, mosquito
#' incubation period and mosquito biting rate into one single CSV file.
#'
#' @param sep CSV field separator character.
#' @param Ns number of samples wanted.
#'
#' @return NULL
#'
exportEntoParameters<- function(sep=',', Ns=200){
  filename<- paste0(.MVSEbuildOutPath(),'.estimated_entomological.csv')

  ttt<-length(MVSE_data$date)
  res<- .calculatePostsVLS_EIP_IIP_HIP(Ns, ttt)
  VLS<- res[[1]]
  EIP<- res[[2]]

  A<- .calculatePostA(Ns, ttt)

  VLS.mean<- apply(VLS, MARG=2, FUN=mean)
  VLS.l<- apply(VLS, MARG=2, FUN=quantile, probs=c(0.025))
  VLS.u<- apply(VLS, MARG=2, FUN=quantile, probs=c(0.975))

  EIP.mean<- apply(EIP, MARG=2, FUN=mean)
  EIP.l<- apply(EIP, MARG=2, FUN=quantile, probs=c(0.025))
  EIP.u<- apply(EIP, MARG=2, FUN=quantile, probs=c(0.975))

  A.mean<- apply(A, MARG=2, FUN=mean)
  A.l<- apply(A, MARG=2, FUN=quantile, probs=c(0.025))
  A.u<- apply(A, MARG=2, FUN=quantile, probs=c(0.975))

  out<- data.frame(date=as.character(MVSE_data$date),
                   eip=EIP.mean,
                   eipu=EIP.u,
                   eipl=EIP.l,
                   lifespan=VLS.mean,
                   lifespanu=VLS.u,
                   lifespanl=VLS.l,
                   bitingrate=A.mean,
                   bitingrateu=A.u,
                   bitingratel=A.l
                 )
  write.table(out,filename,col.names=TRUE,row.names=FALSE,sep=sep)
}



#' Exports the estimated factors Alpha, Rho and Eta into csv files. Four files
#' are exported. One summarises the means and standard deviations of the estimated
#' distributions. The other 3, one for each factor, exports all estimated values
#' that compose the distributions.
#'
#' @param sep CSV field separator character.
#'
#' @return NULL
#'
exportEcoCoefficients<- function(sep=','){
  ALPHAS.mean<- round(mean(MVSE_results$alpha),2)
  ALPHAS.lower<- round(quantile(MVSE_results$alpha, probs=(0.025)),2)
  ALPHAS.upper<- round(quantile(MVSE_results$alpha, probs=(0.975)),2)
  #print(paste0('ALPHA: ', ALPHAS.mean, ' [', ALPHAS.lower, '-', ALPHAS.upper,']'))

  RHOS.mean<- round(mean(MVSE_results$rho),2)
  RHOS.lower<- round(quantile(MVSE_results$rho, probs=(0.025)),2)
  RHOS.upper<- round(quantile(MVSE_results$rho, probs=(0.975)),2)
  #print(paste0('RHO: ', RHOS.mean, ' [', RHOS.lower, '-', RHOS.upper,']'))

  ETAS.mean<- round(mean(MVSE_results$eta),2)
  ETAS.lower<- round(quantile(MVSE_results$eta, probs=(0.025)),2)
  ETAS.upper<- round(quantile(MVSE_results$eta, probs=(0.975)),2)
  #print(paste0('ETA: ', ETAS.mean, ' [', ETAS.lower, '-', ETAS.upper,']'))

  ## export summary of parameters
  eta<- cbind('alpha'=ETAS.mean,'alphal'=ETAS.lower,'alphau'=ETAS.upper)
  rho<- cbind('rho'=RHOS.mean,'rhol'=RHOS.lower,'rhou'=RHOS.upper)
  alpha<- cbind('eta'=ETAS.mean,'etal'=ETAS.lower,'etau'=ETAS.upper)
  out<- data.frame(cbind(eta,rho,alpha))
  write.table(out,paste0(.MVSEbuildOutPath(),'.factors.summary.csv'),col.names=TRUE,row.names=FALSE,sep=sep)

  ## export actual distributions of parameters
  out<- data.frame(alpha=MVSE_results$alpha)
  write.table(out,paste0(.MVSEbuildOutPath(),'.alpha.dist.csv'),col.names=TRUE,row.names=FALSE,sep=sep)
  out<- data.frame(eta=MVSE_results$eta)
  write.table(out,paste0(.MVSEbuildOutPath(),'.eta.dist.csv'),col.names=TRUE,row.names=FALSE,sep=sep)
  out<- data.frame(rho=MVSE_results$rho)
  write.table(out,paste0(.MVSEbuildOutPath(),'.rho.dist.csv'),col.names=TRUE,row.names=FALSE,sep=sep)
}

#' Exports the index P estimations into a csv file. Includes mean P and 95\% CI.
#'
#' @param sep CSV field separator character.
#'
#' @return NULL
#'
exportEmpiricalIndexP<- function(sep=','){
  filename<- paste0(.MVSEbuildOutPath(),'.estimated_indexP.csv')
  out<- data.frame(cbind(date=as.character(MVSE_data$date),indexP=MVSE_results$indexP,indexPlower=MVSE_results$indexPl,indexPupper=MVSE_results$indexPu))
  ##add all smoothing solutions to output
  if(length(MVSE_results$indexPSmooth)>1){
    smoothout<- c()
    for(ss in 1:length(MVSE_results$indexPSmooth)){
      smoothout<- cbind(smoothout, MVSE_results$indexPSmooth[[ss]])
    }
    colnames(smoothout)<- paste0('indexP',names(MVSE_results$indexPSmooth))
    out<- cbind(out, as.data.frame(smoothout))
  }
  write.table(out,filename,col.names=TRUE,row.names=FALSE,sep=sep)
}


#' Exports the factor Q estimations into a csv file. Includes mean Q and 95\% CI.
#'
#' @param sep CSV field separator character.
#'
#' @return NULL
#'
exportEmpiricalQ<- function(sep=','){
    filename<- paste0(.MVSEbuildOutPath(),'.estimated_Q.csv')
    out<- data.frame(cbind(date=as.character(MVSE_data$date),Q=MVSE_results$Q,Qlower=MVSE_results$Ql,Qupper=MVSE_results$Qu))
    ##add all smoothing solutions to output
    if(length(MVSE_results$QSmooth)>1){
      smoothout<- c()
      for(ss in 1:length(MVSE_results$QSmooth)){
        smoothout<- cbind(smoothout, MVSE_results$QSmooth[[ss]])
      }
      colnames(smoothout)<- paste0('Q',names(MVSE_results$QSmooth))
      out<- cbind(out, as.data.frame(smoothout))
    }
    write.table(out,filename,col.names=TRUE,row.names=FALSE,sep=sep)
}


#' Exports the V0 estimations into a csv file. Includes mean V0 and 95\% CI.
#'
#' @param sep CSV field separator character.
#'
#' @return NULL
#'
exportEmpiricalV0<- function(sep=','){
    filename<- paste0(.MVSEbuildOutPath(),'.estimated_V0.csv')
    out<- data.frame(cbind(date=as.character(MVSE_data$date),V0=MVSE_results$V0,V0lower=MVSE_results$V0l,Qupper=MVSE_results$V0u))
    ##add all smoothing solutions to output
    if(length(MVSE_results$V0Smooth)>1){
      smoothout<- c()
      for(ss in 1:length(MVSE_results$V0Smooth)){
        smoothout<- cbind(smoothout, MVSE_results$V0Smooth[[ss]])
      }
      colnames(smoothout)<- paste0('V0',names(MVSE_results$V0Smooth))
      out<- cbind(out, as.data.frame(smoothout))
    }
    write.table(out,filename,col.names=TRUE,row.names=FALSE,sep=sep)
}

###########################################
#### CSV IMPORT FUNCTIONS
###########################################


#' Reads the climatic series and applies some required transformations ('digests').
#' For required formats, please see information in \code{help_inputFormat()}.
#'
#' The effects of climateic variables on entomological parameters are calculated here.
#' These include the effects of temperature on mosquito lifespan and probability of
#' transmission per mosquito bite, and humidity on lifespan and biting rate.
#' Please see Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how
#' these effects are calculated.
#'
#' The function sets a series of global variables which can be
#' used by other functions. To see a list of globally accessible variables use
#' \code{help_globalVariables()}.
#'
#' @param filepath File path to file where the data to be used is.
#' @param D data.frame with the data to use.
#' @param sep Character separator for csv file.
#' @param NDTSmoothTH Number of time points to smooth the temperature and humidity with.
#'
#' @return NULL
#'
setEmpiricalClimateSeries<- function(filepath=NA, D=NA, sep=',', NDTSmoothTH=30){

  MVSE_output_tag<<- ""

  if(is.na(filepath) & !is.data.frame(D)) stop('Both filename and D cannot be undefined')
  if(!is.na(filepath) & is.data.frame(D)) stop('Both filename and D cannot be defined at the same time')

  if(!is.na(filepath) & !is.data.frame(D)){
    MVSE_filepath_input<<- filepath
    MVSE_data<<- read.table(MVSE_filepath_input, header=TRUE,sep=sep,stringsAsFactors=FALSE)
  }

  if(is.na(filepath) & is.data.frame(D)){
    MVSE_filepath_input<<- 'Data was from given matrix.'
    MVSE_data<<- D
  }

  ##temperature, humidity mandatory to calculate P
  # MVSE_data$H<<- .smoothUDSeries(MVSE_data$H,NDTSmoothTH)
  # MVSE_data$T<<- .smoothUDSeries(MVSE_data$T,NDTSmoothTH)
  # MVSE_data$R<<- .smoothUDSeries(MVSE_data$R,NDTSmoothTH)
  MVSE_data$H<<- MVSE_data$H
  MVSE_data$T<<- MVSE_data$T
  MVSE_data$R<<- MVSE_data$R
  MVSE_data$oH<<- MVSE_data$H
  MVSE_data$H<<-  MVSE_data$H/100 #normalize
  MVSE_data$oR<<- MVSE_data$R
  MVSE_data$R<<-  MVSE_data$R/max(MVSE_data$R, na.rm=TRUE) #normalize

  ##essential for P
  MVSE_data$A<<- .MVSEhum_effect_aV(MVSE_data$H,mean(MVSE_data$H))
  MVSE_data$Z<<- .MVSEtemp_effect_muV(MVSE_data$T)
  MVSE_data$Y<<- .MVSEhum_effect_muV(MVSE_data$H,mean(MVSE_data$H))
  MVSE_data$K<<- .MVSEtemp_epsVH(MVSE_data$T)
  MVSE_data$G<<- .MVSEtemp_effect_gammaV(MVSE_data$T)

  ##essential for Q, V0
  MVSE_data$B<<- .MVSEtemp_effect_muA(MVSE_data$T)
  MVSE_data$C<<- .MVSEtemp_effect_epsA(MVSE_data$T)
  MVSE_data$D<<- .MVSEtemp_effect_theta(MVSE_data$T)
  MVSE_data$E<<- .MVSEtemp_effect_C(MVSE_data$T)
  MVSE_data$F<<- .MVSErain_effect_C(MVSE_data$R,mean(MVSE_data$R))

}


#' Reads the climatic series and applies some required transformations ('digests').
#' For required formats, please see information in \code{help_inputFormat()}.
#'
#' The effects of climateic variables on entomological parameters are calculated here.
#' These include the effects of temperature on mosquito lifespan and probability of
#' transmission per mosquito bite, and humidity on lifespan and biting rate.
#' Please see Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how
#' these effects are calculated.
#'
#' The function sets a series of global variables which can be
#' used by other functions. To see a list of globally accessible variables use
#' \code{help_globalVariables()}.
#'
#' @param Trange Numerical range to be considered3 for temperature in Celsius.
#' @param Hrange Numerical range to be considered3 for humidity in the range 0 -100 \%
#' @param NtimePointsSmooth Number of time points to smooth the climatic series
#' with.
#'
#' @return NULL
#'
setTheoreticalClimateSeries<- function(Trange=seq(10,35,length.out=10), Hrange=seq(0,100,length.out=10), NtimePointsSmooth=3){

  expRange<<- expand.grid(T=Trange, H=Hrange)
  MVSE_data_theoretical<<- data.frame(T=as.numeric(expRange$T), H=as.numeric(expRange$H))

  MVSE_data_theoretical$H<<- MVSE_data_theoretical$H
  MVSE_data_theoretical$T<<- MVSE_data_theoretical$T

  MVSE_data_theoretical$oH<<- MVSE_data_theoretical$H
  MVSE_data_theoretical$H<<-   MVSE_data_theoretical$H/100

  MVSE_data_theoretical$A<<- .MVSEhum_effect_aV(MVSE_data_theoretical$H, mean(MVSE_data$H))
  MVSE_data_theoretical$Z<<- .MVSEtemp_effect_muV(MVSE_data_theoretical$T)
  MVSE_data_theoretical$Y<<- .MVSEhum_effect_muV(MVSE_data_theoretical$H, mean(MVSE_data$H))
  MVSE_data_theoretical$K<<- .MVSEtemp_epsVH(MVSE_data_theoretical$T)
  MVSE_data_theoretical$G<<- .MVSEtemp_effect_gammaV(MVSE_data_theoretical$T)

}

###########################################
#### MVSE ESTIMATION FUNCTIONS
###########################################

#' Estimates the factors Alpha, Eta and Rho.
#'
#' For the estimation of Eta and Rho
#' an MCMC routine is run (for which parameters are available). The function sets
#' the \code{MVSE_results} global variable which can be used by other functions.
#' To see a list of globally accessible variables use \code{help_globalVariables()}.
#'
#' @param nMCMC Number of MCMC steps to run the routine for.
#' @param bMCMC Number 0-1 representing the proportion of last steps to be considered3
#' not burn-in. E.g. 0.7 will use the last 0.3 (30\%) of MCMC steps for posteriors.
#' @param cRho Initial guess for factor Rho.
#' @param cEta Initial guess for factor Eta.
#' @param gauJump Gaussian standard deviation for Rho and Eta jumps in the MCMC.
#'
#' @return NULL
#'
estimateEcoCoefficients<- function(nMCMC=30000, bMCMC=0.5, cRho=0.5, cEta=2.0, gauJump=0.05){
  print("Estimating Alpha distribution...")
  ALPHAS<- .estimateFactorAlpha()
  print("Estimating Eta and Rho distributions...")
  listRes<- .estimateFactorsRhoEta(nMCMC, bMCMC, cRho, cEta, gauJump)
  RHOS<- listRes[[1]]
  ETAS<- listRes[[2]]
  MVSE_results<<- list(alpha=ALPHAS,rho=RHOS,eta=ETAS)
  print("Done.")
}


#' Estimates the (one-cycle) generation time; that is, the expected time for one
#' human case to transmit to a vector, and back to a human. This function uses the
#' estimated posteriors of vector lifespan and incubation period, together with the
#' priors for human infectious period and incubation period. The generation time is
#' separated into 2 time contributions: the time it takes for the human to transmit
#' to the vector, and the time it takes for the vector to transmit to a second human.
#' Both of these contributions are simulated and returned to the user. The actual
#' generation time is the sum of the two contributions.
#'
#' It is assumed that transmission may take place at any time point during the
#' infectious period of the human / vector. That is, while the infectious period
#' may be X, the 'effective waiting time for transmission' may be shorter. We take
#' the effective waiting time for transmission of an infectious period X as
#' \code{runif(1, min=0, max=X)}. The human infectious is taken from the prior;
#' the vector infectious period is taken from its lifespan minus the incubation
#' period.
#'
#' Generation time is calculated independently for every time step t. For a particupar
#' step, Ns samples are drawn from the human priors of infectious and incubation
#' periods. Ns samples are also taken for the posteriors of lifespan and incubation
#' period of the vector. The uniform sampling described above is used for every Ns sample.
#' The human and vector contributions to the generation times are just the sum of the
#' incubation period samples and 'effective waiting time for transmission'.
#'
#' The function adds the resulting estimations to \code{MVSE_results} global
#' variable under the name 'Tg' (\code{MVSE_results$Tg}). The output is a data.frame
#' with 22 columns (Tg, Tgu, Tgl, HTg, HTgu, HTgl, VTg, VTgu, VTgl, VLS, VLSu, VLSl,
#' EIP, EIPu, EIPl, HIP, HIPu, HIPl, IIP, IIPu, IIPl, FrNoVecCap). The first 21
#' are organized in sets of 3, with mean, 95 CI upper and lower bounds.
#' Tg stands for generation time, HTg human contribution, VTg vector contribution,
#' VLS vector lifespan, EIP extrinsic incubation period, HIP human infectious period,
#' IIP instrinsic incubation period. The last column FrNoVecCap is the proportion of
#' samples per time step which had vector lifespan longer than the incubation period.
#'
#'
#' @param Ns Number of samples to be drawn from priors and posteriors per time step.
#'
#' @return NULL
#'
simulateGenerationTime<- function(Ns= 200){
        print('Simulating generation times...')

        ttt<-length(MVSE_data$date)
        res<- .calculatePostsVLS_EIP_IIP_HIP(Ns, ttt)
        VLS<- res[[1]]
        EIP<- res[[2]]
        IIP<- res[[3]]
        HIP<- res[[4]]

        print('Organizing output...')
        sampleIP<- function(X){
          return(runif(length(X), min=0, max=X))
        }

        freqNoVecCapacity<- function(X){
          return(sum(X<=0)/length(X))
        }

        VLS.mean<- apply(VLS, MARG=2, FUN=mean)
        VLS.u<- apply(VLS, MARG=2, FUN=quantile, probs=c(0.975))
        VLS.l<- apply(VLS, MARG=2, FUN=quantile, probs=c(0.025))
        EIP.mean<- apply(EIP, MARG=2, FUN=mean)
        EIP.u<- apply(EIP, MARG=2, FUN=quantile, probs=c(0.975))
        EIP.l<- apply(EIP, MARG=2, FUN=quantile, probs=c(0.025))

        HIP.mean<- apply(HIP, MARG=2, FUN=mean)
        HIP.u<- apply(HIP, MARG=2, FUN=quantile, probs=c(0.975))
        HIP.l<- apply(HIP, MARG=2, FUN=quantile, probs=c(0.025))
        IIP.mean<- apply(IIP, MARG=2, FUN=mean)
        IIP.u<- apply(IIP, MARG=2, FUN=quantile, probs=c(0.975))
        IIP.l<- apply(IIP, MARG=2, FUN=quantile, probs=c(0.025))

            VEffIP<- VLS-EIP
            VFrNoCap<- apply(VEffIP, MARG=2, FUN=freqNoVecCapacity)
              ##use only the ones with vectorical capacity (LS-EIP>0)
              getVEffIpWithVecCap<- function(X, VEffIP, EIP){
                  ip<- VEffIP[,X]
                  eip<- EIP[,X]
                  withVecCap<- which(ip>0)
                  ip<- ip[withVecCap]
                  eip<- eip[withVecCap]
                  return(eip+ip)
              }
              VRealTg<- lapply(1:ttt, FUN=getVEffIpWithVecCap, VEffIP, EIP)
              VRealTg.mean<- sapply(VRealTg, FUN=mean)
              VRealTg.u<- sapply(VRealTg, FUN=quantile, probs=c(0.975))
              VRealTg.l<- sapply(VRealTg, FUN=quantile, probs=c(0.025))
              ##use only the ones with vectorical capacity (LS-EIP>0)
              getHEffIpWithVecCap<- function(X, VEffIP, IIP, HIP){
                  ip<- VEffIP[,X]
                  iip<- IIP[,X]
                  hip<- HIP[,X]
                  withVecCap<- which(ip>0)
                  iip<- iip[withVecCap]
                  hip<- hip[withVecCap]
                  return(iip+hip)
              }
              HRealTg<- lapply(1:ttt, FUN=getHEffIpWithVecCap, VEffIP, IIP, HIP)
              HRealTg.mean<- sapply(HRealTg, FUN=mean)
              HRealTg.u<- sapply(HRealTg, FUN=quantile, probs=c(0.975))
              HRealTg.l<- sapply(HRealTg, FUN=quantile, probs=c(0.025))

        addTgs<- function(X, HRealTg, VRealTg){
          return(HRealTg[[X]]+VRealTg[[X]])
        }
        RealTg<- sapply(1:ttt, FUN=addTgs, HRealTg, VRealTg)
        RealTg.mean<- sapply(RealTg, FUN=mean)
        RealTg.u<- sapply(RealTg, FUN=quantile, probs=c(0.975))
        RealTg.l<- sapply(RealTg, FUN=quantile, probs=c(0.025))

        TgDF<- data.frame(Tg=RealTg.mean, Tgu=RealTg.u, Tgl=RealTg.l,
                          HTg=HRealTg.mean, HTgu=HRealTg.u, HTgl=HRealTg.l,
                          VTg=VRealTg.mean, VTgu=VRealTg.u, VTgl=VRealTg.l,
                          VLS=VLS.mean, VLSu=VLS.u, VLSl=VLS.l,
                          EIP=EIP.mean, EIPu=EIP.u, EIPl=EIP.l,
                          HIP=HIP.mean, HIPu=HIP.u, HIPl=HIP.l,
                          IIP=IIP.mean, IIPu=IIP.u, IIPl=IIP.l,
                          FrNoVecCap= VFrNoCap
                        )

        MVSE_results$Tg<<- TgDF
}

#' Estimates the index P by simulation using a theoretical climate range.
#' That is, the estimated posteriors of factors Eta and Rho are sampled,
#' as well as the distribution of factor Alpha. These independent samples are then
#' used to calculate index P given its mathematical expression. Since this
#' expression is dependent on climatic variables in time, the result is a time
#' series of the index P. The user can choose the number of samples,
#' resulting in \code{nSample} estimations of P in time. These independent
#' estimations are then used to calculate the mean and 95\% CI of the index P.
#'
#' In contrast to the function \code{simulateEmpiricalIndexP} which uses the
#' observed climatic variables, this function uses a theoretical range for climate.
#' The range can be set in the function's parameters.
#'
#' The function adds the resulting index P estimations to the
#' \code{MVSE_results} global variable which can be used by other functions.
#' All independent estimations of index P in time are also stored3
#' in the \code{MVSE_indexP_theoretical} global variable (matrix, rows as independent estimations,
#' cols as time steps). To see a list of globally accessible variables use \code{help_globalVariables()}
#'
#' The index P follows the expression of R0 from Lourenco et al. doi.org/10.7554/eLife.29820.001
#' ignoring the number of mosquitoes per human (V/N).
#'
#' @param nSample Number of samples to draw from factors Alpha, Eta and Rho.
#'
#' @return NULL
#'
simulateTheoreticalIndexP<- function(nSample=1000){

  sN<- nSample
  ss<- sample(1:length(MVSE_results$alpha), sN)
  sMCMC_alpha<- MVSE_results$alpha[ss]
  sMCMC_rhos<- MVSE_results$rho[ss]
  sMCMC_etas<- MVSE_results$eta[ss]

  print("Simulating theoretical index P given distributions of Eta, Alpha and Rho...")

  Nsteps<- nrow(MVSE_data_theoretical)
  indexP<- rep(NA,sN*(Nsteps))
  indexP<- matrix(indexP, ncol=(Nsteps))
  pb <- txtProgressBar(min=0, max = sN-1, style = 3)
  for(ss in 1:sN){
    #for each sample of the factors' distributions
    eta<- sMCMC_etas[ss]
    rho<- sMCMC_rhos[ss]
    alpha<- sMCMC_alpha[ss]
    #sample other priors
    muH<- 1/MVSE_prior_muH_rdist(1)
    gammaH<- 1/MVSE_prior_gammaH_rdist(1)
    deltaH<- 1/MVSE_prior_deltaH_rdist(1)
    epsilonHV<- MVSE_prior_epsilonHV_rdist(1)
    ##calculate ento-epi parameters
    muV_t<- eta*getTempEffLifespan(tpe='the')*(1+getHumEffLifespan(tpe='the'))^rho
    gammaV_t <- alpha*getTempEffIncPer(tpe='the')
    a_t<- MVSE_prior_a_mean*(1+MVSE_data_theoretical$A)^rho
    epsilonVH_t<- getTempEffProbTransVH(tpe='the')
    betaHV<- a_t*epsilonVH_t
    betaVH<- a_t*epsilonHV
    ##correct negative to zero for bio meaning
    muV_t[which(muV_t<0)]<- 0
    gammaV_t[which(gammaV_t<0)]<- 0
    a_t[which(a_t<0)]<- 0
    epsilonVH_t[which(epsilonVH_t<0)]<- 0
    betaHV[which(betaHV<0)]<- 0
    betaVH[which(betaVH<0)]<- 0
    #calculate P
    indexP[ss,]<- (betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t)) #index P
    indexP[ss,which(indexP[ss,]<0)]<- 0
    setTxtProgressBar(pb, ss)
  }
  close(pb)

  MVSE_results$indexP_theoretical<<- colMeans(indexP)
  MVSE_results$indexPl_theoretical<<- apply(indexP, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.025), na.rm=TRUE)} )
  MVSE_results$indexPu_theoretical<<- apply(indexP, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.975), na.rm=TRUE)} )
  MVSE_indexP_theoretical<<- indexP

}

#' Estimates the index P by simulation. That is, the estimated posteriors of factors
#' Eta and Rho are sampled, as well as the distribution of factor Alpha. These
#' independent samples are then used to calculate index P given its mathematical
#' expression. Since this expression is dependent on climatic variables in time,
#' the result is a time series of the index P. The user can choose the number of
#' samples, resulting in \code{nSample} estimations of P in time. These independent
#' estimations are then used to calculate the mean and 95\% CI of the index P.
#'
#' The function adds the resulting index P estimations to the
#' \code{MVSE_results} global variable which can be used by other functions. Smoothed,
#' averaged estimations can also be added to \code{MVSE_results} using the parameter
#' \code{smoothing}. All independent estimations of index P in time are also stored3
#' in the \code{MVSE_indexP} global variable (matrix, rows as independent estimations,
#' cols as time steps). To see a list of globally accessible variables use
#' \code{help_globalVariables()}
#'
#' The index P follows the expression of R0 from Lourenco et al. doi.org/10.7554/eLife.29820.001
#' ignoring the number of mosquitoes per human (V/N).
#'
#' @param nSample Number of samples to draw from factors Alpha, Eta and Rho.
#' @param smoothing Set of integers that define how many averaging smoothers should
#' be calculated for the index P. E.g. c(7, 30) will calculate index P series using
#' +- 7 and +- 30 time points per existing point (up and down that point).
#'
#' @return NULL
#'
simulateEmpiricalIndexP<- function(nSample=1000, smoothing=c()){

  sN<- nSample
  ss<- sample(1:length(MVSE_results$alpha), size=sN)
  sMCMC_alpha<- MVSE_results$alpha[ss]
  sMCMC_rhos<- MVSE_results$rho[ss]
  sMCMC_etas<- MVSE_results$eta[ss]

  print("Simulating empirical index P, Q, V0 given distributions of Eta, Alpha and Rho...")

  Nsteps<- nrow(MVSE_data)
  repvals<- rep(NA,sN*(Nsteps))
  indexP<- matrix(repvals, ncol=(Nsteps))
  Q<- matrix(repvals, ncol=(Nsteps))
  V0<- matrix(repvals, ncol=(Nsteps))

  pb <- txtProgressBar(min=0, max = sN-1, style = 3)
  for(ss in 1:sN){
    #for each sample of the factors' distributions
    eta<- sMCMC_etas[ss]
    rho<- sMCMC_rhos[ss]
    alpha<- sMCMC_alpha[ss]
    #sample other priors
    muH<- 1/MVSE_prior_muH_rdist(1)
    gammaH<- 1/MVSE_prior_gammaH_rdist(1)
    deltaH<- 1/MVSE_prior_deltaH_rdist(1)
    epsilonHV<- MVSE_prior_epsilonHV_rdist(1)
    ##calculate ento-epi parameters
    muV_t<- eta*getTempEffLifespan()*(1+getHumEffLifespan())^rho
    gammaV_t <- alpha*getTempEffIncPer()
    a_t<- MVSE_prior_a_mean*(1+MVSE_data$A)^rho
    epsilonVH_t<- getTempEffProbTransVH()
    betaHV<- a_t*epsilonVH_t
    betaVH<- a_t*epsilonHV
    # print(paste(epsilonHV,muV_t[1],gammaV_t[1],a_t[1],epsilonVH_t[1]))
    ##correct negative to zero for bio meaning
    muV_t[which(muV_t<0)]<- 0
    gammaV_t[which(gammaV_t<0)]<- 0
    a_t[which(a_t<0)]<- 0
    epsilonVH_t[which(epsilonVH_t<0)]<- 0
    betaHV[which(betaHV<0)]<- 0
    betaVH[which(betaVH<0)]<- 0

    #calculate P
    indexP[ss,]<- (betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t)) #index P
    indexP[ss,which(indexP[ss,]<0)]<- 0

    ##calculate Q
    f<- 0.5 ##TODO: this should perhaps not be hardcoded?
    muA_t<- eta*getTempEffAquaLifeSpan()
    EA_t<- pmax(0, getTempEffAquaDev())
    THETA_t<- pmax(0, getTempEffOviPos())
    CS_t<- pmin(getTempEffEggEcc()*(1+getRainEffEggEcc())^rho,1)
    Q[ss,]<- pmax(0, EA_t/(EA_t+muA_t)*(f*THETA_t*CS_t)/muV_t)

    ##calculate V0
    V0[ss,]<- pmax(0, (MVSE_data$R+1)*(1-1/Q[ss,])*(EA_t/muV_t) )

    setTxtProgressBar(pb, ss)
  }
  close(pb)

  MVSE_results$indexP<<- colMeans(indexP)
  MVSE_results$indexPl<<- apply(indexP, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.025), na.rm=T)} )
  MVSE_results$indexPu<<- apply(indexP, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.975), na.rm=T)} )
  MVSE_indexP<<- indexP

  MVSE_results$Q<<- colMeans(Q)
  MVSE_results$Ql<<- apply(Q, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.025), na.rm=T)} )
  MVSE_results$Qu<<- apply(Q, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.975), na.rm=T)} )
  MVSE_Q<<- Q

  MVSE_results$V0<<- colMeans(V0)
  MVSE_results$V0l<<- apply(V0, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.025), na.rm=T)} )
  MVSE_results$V0u<<- apply(V0, MARGIN=2, FUN= function(X){ quantile(X, probs=c(0.975), na.rm=T)} )
  MVSE_V0<<- V0

  if( length(smoothing)>0 ){
    MVSE_results$indexPSmooth<<- list()
    MVSE_results$QSmooth<<- list() ###NEW
    MVSE_results$V0Smooth<<- list() ###NEW
    for(ss in smoothing){
      ii<- length(MVSE_results$indexPSmooth)+1
      # MVSE_results$indexPSmooth[[ii]]<<- .smoothUDSeries(MVSE_results$indexP, ss)
      # MVSE_results$QSmooth[[ii]]<<- .smoothUDSeries(MVSE_results$Q, ss)
      # MVSE_results$V0Smooth[[ii]]<<- .smoothUDSeries(MVSE_results$V0, ss)
      MVSE_results$indexPSmooth[[ii]]<<- MVSE_results$indexP
      MVSE_results$QSmooth[[ii]]<<- MVSE_results$Q
      MVSE_results$V0Smooth[[ii]]<<- MVSE_results$V0
    }
    names(MVSE_results$indexPSmooth)<<- paste0('smooth',smoothing)
    names(MVSE_results$QSmooth)<<- paste0('smooth',smoothing)
    names(MVSE_results$V0Smooth)<<- paste0('smooth',smoothing)
  }else{
    ##user does not want smoothing to be done and any previous smoothing should be deleted
    MVSE_results$indexPSmooth<<- list()
    MVSE_results$QSmooth<<- list()
    MVSE_results$V0Smooth<<- list()
  }

}


###########################################
#### SMOOTHERS AND FITTERS
###########################################

#' Calculates the mean behaviour of a year using the climatic data and the estimated P
#' for several years (if the input is for several years).
#'
#' @return data.frame with colums mean temperature (meanT), stdev temperature (sdT),
#' mean humidity (meanH), stdev humidity (sdH), mean index P (meanP), stdev index P (sdP)
#'
getTypicalYear<- function(){

  monthday<- (format(as.Date(MVSE_data$date, format="%Y-%m-%d"),"%m-%d"))
  uMD<- unique(monthday)
  uMD<- uMD[-(which(uMD=="02-29"))] #for simplicity

  typicalYear<- c()
  for(umd in uMD){
    indexes<- which(monthday==umd)

    meanT<- mean(MVSE_data$T[indexes])
    sdT<- sd(MVSE_data$T[indexes])
    meanH<- mean(MVSE_data$H[indexes])
    sdH<- sd(MVSE_data$H[indexes])

    meanP<- mean(MVSE_results$indexP[indexes])
    sdP<- sd(MVSE_results$indexP[indexes])

    typicalYear<- rbind(typicalYear, c(meanT, sdT, meanH, sdH, meanP, sdP))
  }
  class(typicalYear)<- 'numeric'
  typicalYear<- data.frame(typicalYear, stringsAsFactors=FALSE)
  colnames(typicalYear)<- c("meanT", "sdT", "meanH", "sdH", "meanP", "sdP")
  # typicalYear$date<- as.Date(uMD, format='%m-%d')
  typicalYear$date<- uMD
  return(typicalYear)
}

#' Calculates the mean monthly index P.
#'
#' @return data.frame with colums mean index P (meanP), stdev index P (sdP), date
#' in the form 'year-month' (YM).
#'
getMonthlyIndexP<- function(){

  yearmonth<- (format(as.Date(MVSE_data$date, format="%Y-%m-%d"),"%Y-%m"))
  uYM<- unique(yearmonth)

  result<- c()
  for(uym in uYM){
    indexes<- which(yearmonth==uym)

    meanP<- mean(MVSE_results$indexP[indexes])
    sdP<- sd(MVSE_results$indexP[indexes])

    result<- rbind(result, c(meanP, sdP))
  }
  class(result)<- 'numeric'
  result<- data.frame(result, stringsAsFactors=FALSE)
  result$YM<- uYM
  colnames(result)<- c("meanP", "sdP", "YM")
  return(result)
}

###########################################
#### GLOBAL MVSE ENVIRONMENT FUNCTIONS
###########################################


#' Set the prior for the mosquito lifespan.
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setMosqLifeExpPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_lev_mean<<- pmean
  MVSE_prior_lev_sd<<- psd
  MVSE_prior_lev_dist<<- pdist
  if(MVSE_prior_lev_dist=='gaussian'){
    MVSE_prior_lev_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_lev_mean, sd=MVSE_prior_lev_sd) }
    MVSE_prior_lev_ddist<<- function(X) { dnorm(X, mean=MVSE_prior_lev_mean, sd=MVSE_prior_lev_sd) }
  }else if(MVSE_prior_lev_dist=='gamma'){
    MVSE_prior_lev_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_lev_mean*MVSE_prior_lev_mean/(MVSE_prior_lev_sd*MVSE_prior_lev_sd), scale=(MVSE_prior_lev_sd*MVSE_prior_lev_sd)/MVSE_prior_lev_mean) }
    MVSE_prior_lev_ddist<<- function(X) { dgamma(X, shape=MVSE_prior_lev_mean*MVSE_prior_lev_mean/(MVSE_prior_lev_sd*MVSE_prior_lev_sd), scale=(MVSE_prior_lev_sd*MVSE_prior_lev_sd)/MVSE_prior_lev_mean) }
  }else if(MVSE_prior_lev_dist=='exponential'){
    MVSE_prior_lev_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_lev_mean) }
    MVSE_prior_lev_ddist<<- function(X) { dexp(X, rate=1/MVSE_prior_lev_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}

#' Get the prior for the mosquito lifespan.
#'
getMosqLEprior<- function(){
  print(paste('mean:',MVSE_prior_lev_mean,'sd:',MVSE_prior_lev_sd,'dist:',MVSE_prior_lev_dist))
}

#' Set the prior for the mosquito incubation period.
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setMosqIncPerPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_ic_mean<<- pmean
  MVSE_prior_ic_sd<<- psd
  MVSE_prior_ic_dist<<- pdist
  if(MVSE_prior_ic_dist=='gaussian'){
    MVSE_prior_ic_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_ic_mean, sd=MVSE_prior_ic_sd) }
    MVSE_prior_ic_ddist<<- function(X) { dnorm(X, mean=MVSE_prior_ic_mean, sd=MVSE_prior_ic_sd) }
  }else if(MVSE_prior_ic_dist=='gamma'){
    MVSE_prior_ic_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_ic_mean*MVSE_prior_ic_mean/(MVSE_prior_ic_sd*MVSE_prior_ic_sd), scale=(MVSE_prior_ic_sd*MVSE_prior_ic_sd)/MVSE_prior_ic_mean) }
    MVSE_prior_ic_ddist<<- function(X) { dgamma(X, shape=MVSE_prior_ic_mean*MVSE_prior_ic_mean/(MVSE_prior_ic_sd*MVSE_prior_ic_sd), scale=(MVSE_prior_ic_sd*MVSE_prior_ic_sd)/MVSE_prior_ic_mean) }
  }else if(MVSE_prior_ic_dist=='exponential'){
    MVSE_prior_ic_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_ic_mean) }
    MVSE_prior_ic_ddist<<- function(X) { dexp(X, rate=1/MVSE_prior_ic_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}

#' Get the prior for the mosquito incubation period.
#'
getMosqIPprior<- function(){
  print(paste('mean:',MVSE_prior_ic_mean,'sd:',MVSE_prior_ic_sd,'dist:',MVSE_prior_ic_dist))
}

#' Set the prior for the mosquito biting rate.
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setMosqBitingPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_a_mean<<- pmean
  MVSE_prior_a_sd<<- psd
  MVSE_prior_a_dist<<- pdist
  if(MVSE_prior_a_dist=='gaussian'){
    MVSE_prior_a_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_a_mean, sd=MVSE_prior_a_sd) }
    MVSE_prior_a_ddist<<- function(X) { dnorm(X, mean=MVSE_prior_a_mean, sd=MVSE_prior_a_sd) }
  }else if(MVSE_prior_a_dist=='gamma'){
    MVSE_prior_a_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_a_mean*MVSE_prior_a_mean/(MVSE_prior_a_sd*MVSE_prior_a_sd), scale=(MVSE_prior_a_sd*MVSE_prior_a_sd)/MVSE_prior_a_mean) }
    MVSE_prior_a_ddist<<- function(X) { dgamma(X, shape=MVSE_prior_a_mean*MVSE_prior_a_mean/(MVSE_prior_a_sd*MVSE_prior_a_sd), scale=(MVSE_prior_a_sd*MVSE_prior_a_sd)/MVSE_prior_a_mean) }
  }else if(MVSE_prior_a_dist=='exponential'){
    MVSE_prior_a_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_a_mean) }
    MVSE_prior_a_ddist<<- function(X) { dexp(X, rate=1/MVSE_prior_a_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}


#' Get the prior for the mosquito biting rate.
#'
#'
getMosqBRprior<- function(){
  print(paste('mean:',MVSE_prior_a_mean,'sd:',MVSE_prior_a_sd,'dist:',MVSE_prior_a_dist))
}

#' Set the file name and path of the climatic series input.
#'
#' @param tag A character string.
#'
#' @return NULL
#'
setOutputFilePathAndTag<- function(tag){
  MVSE_output_tag<<- tag
}


#' Get the file name and path of the climatic series input.
#'
getClimateInputFile<- function(){
  print(MVSE_filepath_input)
}

#' Set the prior for the human life-span.
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setHumanLifeExpPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_muH_mean<<- (pmean*365)
  MVSE_prior_muH_sd<<- (psd*365)
  MVSE_prior_muH_dist<<- pdist
  if(MVSE_prior_muH_dist=='gaussian'){
    MVSE_prior_muH_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_muH_mean, sd=MVSE_prior_muH_sd) }
  }else if(MVSE_prior_muH_dist=='gamma'){
    MVSE_prior_muH_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_muH_mean*MVSE_prior_muH_mean/(MVSE_prior_muH_sd*MVSE_prior_muH_sd), scale=(MVSE_prior_muH_sd*MVSE_prior_muH_sd)/MVSE_prior_muH_mean) }
  }else if(MVSE_prior_muH_dist=='exponential'){
    MVSE_prior_muH_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_muH_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}



#' Set the prior for the human incubation period.
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setHumanIncPerPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_gammaH_mean<<- pmean
  MVSE_prior_gammaH_sd<<- psd
  MVSE_prior_gammaH_dist<<- pdist
  if(MVSE_prior_gammaH_dist=='gaussian'){
    MVSE_prior_gammaH_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_gammaH_mean, sd=MVSE_prior_gammaH_sd) }
  }else if(MVSE_prior_gammaH_dist=='gamma'){
    MVSE_prior_gammaH_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_gammaH_mean*MVSE_prior_gammaH_mean/(MVSE_prior_gammaH_sd*MVSE_prior_gammaH_sd), scale=(MVSE_prior_gammaH_sd*MVSE_prior_gammaH_sd)/MVSE_prior_gammaH_mean) }
  }else if(MVSE_prior_gammaH_dist=='exponential'){
    MVSE_prior_gammaH_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_gammaH_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}

#' Set the prior for the human infectious period.
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setHumanInfPerPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_deltaH_mean<<- pmean
  MVSE_prior_deltaH_sd<<- psd
  MVSE_prior_deltaH_dist<<- pdist
  if(MVSE_prior_deltaH_dist=='gaussian'){
    MVSE_prior_deltaH_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_deltaH_mean, sd=MVSE_prior_deltaH_sd) }
  }else if(MVSE_prior_deltaH_dist=='gamma'){
    MVSE_prior_deltaH_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_deltaH_mean*MVSE_prior_deltaH_mean/(MVSE_prior_deltaH_sd*MVSE_prior_deltaH_sd), scale=(MVSE_prior_deltaH_sd*MVSE_prior_deltaH_sd)/MVSE_prior_deltaH_mean) }
  }else if(MVSE_prior_deltaH_dist=='exponential'){
    MVSE_prior_deltaH_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_deltaH_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}

#' Set the prior human to mosquito probability of transmission per infectious bite
#'
#' @param pmean A number.
#' @param psd A number.
#' @param pdist A string naming the desired3 distribution: implemented 'gaussian'
#' 'exponential' and 'gamma'.
#'
#' @return NULL
#'
setHumanMosqTransProbPrior<- function(pmean, psd, pdist){
  ##setting prior mean and sd into global variables
  MVSE_prior_epsilonHV_mean<<- pmean
  MVSE_prior_epsilonHV_sd<<- psd
  MVSE_prior_epsilonHV_dist<<- pdist
  if(MVSE_prior_epsilonHV_dist=='gaussian'){
    MVSE_prior_epsilonHV_rdist<<- function(N) { rnorm(N, mean=MVSE_prior_epsilonHV_mean, sd=MVSE_prior_epsilonHV_sd) }
  }else if(MVSE_prior_epsilonHV_dist=='gamma'){
    MVSE_prior_epsilonHV_rdist<<- function(N) { rgamma(N, shape=MVSE_prior_epsilonHV_mean*MVSE_prior_epsilonHV_mean/(MVSE_prior_epsilonHV_sd*MVSE_prior_epsilonHV_sd), scale=(MVSE_prior_epsilonHV_sd*MVSE_prior_epsilonHV_sd)/MVSE_prior_epsilonHV_mean) }
  }else if(MVSE_prior_epsilonHV_dist=='exponential'){
    MVSE_prior_epsilonHV_rdist<<- function(N) { rexp(N, rate=1/MVSE_prior_epsilonHV_mean) }
  }else{
    stop(paste0('Type of distribution "',pdist,'" not supported.'))
  }
}



#' Gets the time series for the effect of temperature on aquatic lifespan. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffAquaLifeSpan<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$B)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$B[t])
  if(tpe=='the') stop('Theoretical option is not implemented for Q and V0 calculations.')
}

#' Gets the time series for the effect of temperature on aquatic development rate. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffAquaDev<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$C)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$C[t])
  if(tpe=='the') stop('Theoretical option is not implemented for Q and V0 calculations.')
}

#' Gets the time series for the effect of temperature on oviposition. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffOviPos<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$D)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$D[t])
  if(tpe=='the') stop('Theoretical option is not implemented for Q and V0 calculations.')
}

#' Gets the time series for the effect of temperature on the egg success rate. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffEggEcc<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$E)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$E[t])
  if(tpe=='the') stop('Theoretical option is not implemented for Q and V0 calculations.')
}

#' Gets the time series for the effect of humidity on the egg success rate. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getRainEffEggEcc<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$F)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$F[t])
  if(tpe=='the') stop('Theoretical option is not implemented for Q and V0 calculations.')
}


#' Gets the time series for the effect of humidity on the biting rate. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getHumEffBitingRate<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$A)
  if(tpe=='the' &  is.na(t)) return(MVSE_data_theoretical$A)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$A[t])
  if(tpe=='the' &  !is.na(t)) return(MVSE_data_theoretical$A[t])
}


#' Gets the time series for the effect of humidity on mosquito lifespan. Please see
#' Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this effect
#' is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getHumEffLifespan<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$Y)
  if(tpe=='the' &  is.na(t)) return(MVSE_data_theoretical$Y)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$Y[t])
  if(tpe=='the' &  !is.na(t)) return(MVSE_data_theoretical$Y[t])
}

#' Gets the time series for the effect of temperature on mosquito lifespan. Please
#' see Lourenco et al. doi.org/10.7554/eLife.29820.001 for details on how this
#' effect is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffLifespan<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$Z)
  if(tpe=='the' &  is.na(t)) return(MVSE_data_theoretical$Z)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$Z[t])
  if(tpe=='the' &  !is.na(t)) return(MVSE_data_theoretical$Z[t])
}

#' Gets the time series for the effect of temperature on the probability of transmission
#' from vector to human. Please see Lourenco et al. doi.org/10.7554/eLife.29820.001
#' for details on how this effect is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffProbTransVH<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$K)
  if(tpe=='the' &  is.na(t)) return(MVSE_data_theoretical$K)
  if(tpe=='emp' &  !is.na(t)) return(MVSE_data$K[t])
  if(tpe=='the' &  !is.na(t)) return(MVSE_data_theoretical$K[t])
}

#' Gets the time series for the effect of temperature on the extrinsic incubation
#' period (mosquito). Please see Lourenco et al. doi.org/10.7554/eLife.29820.001
#' for details on how this effect is calculated.
#'
#' @param tpe defines the type of query: 'emp' (default) for empirical range; any other
#' value will return the theoretical range (used in the index P suitability maps).
#' @param t should be numeric and if not NA can be used to obtain the value at a
#' particular time step t.
#'
#' @return effect as array.
#'
getTempEffIncPer<- function(tpe='emp', t=NA){
  if(tpe=='emp' &  is.na(t)) return(MVSE_data$G)
  if(tpe=='the' &  is.na(t)) return(MVSE_data_theoretical$G)
  if(tpe=='emp' & !is.na(t)) return(MVSE_data$G[t])
  if(tpe=='the' & !is.na(t)) return(MVSE_data_theoretical$G[t])
}


#' Gets the posterior of the vector incubation period. User can choose particular
#' time step t and use mEtaRho samples to calculate the incubation period from the
#' posterior of the ecological coefficient alpha.
#' If the number of samples is not given, it is assumed the user wants a 'mean'
#' estimation, and the mean of the alpha ecological coefficient is used over all
#' time steps.
#'
#' @param mAlpha if given and numeric, sets the number of samples to be taken from
#' from the alpha ecological parameter that is used to calculate the incubation
#' period.
#' @param t if given and numeric, sets the desired time step to get a distribution
#' of the incubation period on that step.
#'
#' @return incubation period.
#'
getPosteriorAVForTimet<- function(mRho=NA, t=NA){
  if(is.na(mRho)) pA<- (MVSE_prior_a_mean*(1+getHumEffBitingRate())^mean(MVSE_results$rho))
  if(!is.na(mRho)){
     rr<- sample(MVSE_results$rho, mRho, replace=TRUE)
     if(!is.na(t)) pA<- (MVSE_prior_a_mean*(1+getHumEffBitingRate(t=t))^rr)
     if( is.na(t)) stop('Parameter t needs to be given.')
   }
  return(pA)
}


#' Gets the posterior of the vector incubation period. User can choose particular
#' time step t and use mEtaRho samples to calculate the incubation period from the
#' posterior of the ecological coefficient alpha.
#' If the number of samples is not given, it is assumed the user wants a 'mean'
#' estimation, and the mean of the alpha ecological coefficient is used over all
#' time steps.
#'
#' @param mAlpha if given and numeric, sets the number of samples to be taken from
#' from the alpha ecological parameter that is used to calculate the incubation
#' period.
#' @param t if given and numeric, sets the desired time step to get a distribution
#' of the incubation period on that step.
#'
#' @return incubation period.
#'
getPosteriorGammaVForTimet<- function(mAlpha=NA, t=NA){
  if(is.na(mAlpha)) pGam<- 1/(mean(MVSE_results$alpha)*getTempEffIncPer())
  if(!is.na(mAlpha)){
     aa<- sample(MVSE_results$alpha, mAlpha, replace=TRUE)
     if(!is.na(t)) pGam<- 1/(aa*getTempEffIncPer(t=t))
     if( is.na(t)) stop('Parameter t needs to be given.')
   }
  return(pGam)
}

#' Gets the posterior of the vector life span. User can choose particular
#' time step t and use mEtaRho samples to calculate the life span from the posteriors
#' of the ecological coefficients eta and rho.
#' If the number of samples is not given, it is assumed the user wants a 'mean'
#' estimation, and the mean of the alpha ecological coefficient is used over all
#' time steps.
#'
#' @param mEtaRho if given and numeric, sets the number of samples to be taken from
#' from the eta and rho ecological parameters that are used to calculate the
#' life span.
#' @param t if given and numeric, sets the desired time step to get a distribution
#' of the life span on that step.
#'
#' @return incubation period.
#'
getPosteriorLifespanVForTimet<- function(mEtaRho=NA, t=NA){
  if(is.na(mEtaRho)){
      pMu<- 1/(mean(MVSE_results$eta)*getTempEffLifespan()*(1+getHumEffLifespan())^mean(MVSE_results$rho))
  }else{
    if(!is.na(mEtaRho)){
      if(is.na(t)) stop('Parameter t needs to be given.')
      ss<- sample(1:length(MVSE_results$eta), mEtaRho, replace=TRUE)
      ee<- MVSE_results$eta[ss]
      rr<- MVSE_results$rho[ss]
      pMu<- 1/(ee*getTempEffLifespan(t=t)*(1+getHumEffLifespan(t=t))^rr)
    }else{
      stop('Both mEta and mRho have to be given if one of them already is.')
    }
   }
  return(pMu)
}



###########################################
#### PLOTTING FUNCTIONS
###########################################

#' Plots and exports generation time and variables related to its calculation.
#'
#' @param outfilename A name tag for the output file.
#' @param entoLim two numbers setting the y-axis limits for the entomological solutions
#' @param humLim two numbers setting the y-axis limits for the human solutions
#' @param genLim two numbers setting the y-axis limits for the generation solutions
#'
#' @return NULL
#'
plotGenerationTimes<- function(outfilename='debug_genTimes', entoLim= c(), humLim= c(), genLim= c()){

  if(length(entoLim)==0){
    entoLim<- c(0,max(MVSE_results$Tg$VLSu, na.rm=TRUE))
  }

  if(length(humLim)==0){
    humLim<- c(0,max(MVSE_results$Tg$HIPu, na.rm=TRUE))
  }

  if(length(genLim)==0){
    genLim<- c(0,max(MVSE_results$Tg$Tgu, na.rm=TRUE))
  }

  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=10,h=6,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)
    MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')
    plot(MVSE_data$date, MVSE_results$Tg$VLS, t='l',xaxt='n', col='black', ylim=entoLim, main='vector-related variables', ylab='values / ranges', xlab='time')
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$VLSl, MVSE_results$Tg$VLSu, colIn=alpha('black',0.4), colLines=NA)
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$EIPl, MVSE_results$Tg$EIPu, colIn=alpha('purple',0.4), colLines=NA)
    lines(MVSE_data$date, MVSE_results$Tg$VLS, t='l', col='black')
    lines(MVSE_data$date, MVSE_results$Tg$EIP, t='l', col='purple')
    labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
    axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
    legend('topleft',legend=c('lifespan mean & 95% CI','incubation period mean & 95% CI'),col=c('black','purple'), lwd=2, bg='white')
    box(lwd=1.5)

    plot(MVSE_data$date,MVSE_results$Tg$HIP, t='l',xaxt='n', col='black', ylim=humLim, main='human-related variables', ylab='values / ranges', xlab='time')
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$HIPl, MVSE_results$Tg$HIPu, colIn=alpha('red3',0.4), colLines=NA)
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$IIPl, MVSE_results$Tg$IIPu, colIn=alpha('blue3',0.4), colLines=NA)
    lines(MVSE_data$date,MVSE_results$Tg$HIP, t='l', col='red3')
    lines(MVSE_data$date,MVSE_results$Tg$IIP, t='l', col='blue3')
    labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
    axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
    legend('bottomleft',legend=c('incubation period mean & 95% CI','infectious period mean & 95% CI'),col=c('red3','blue3'), lwd=2, bg='white')
    box(lwd=1.5)

    plot(MVSE_data$date, MVSE_results$Tg$Tg, t='l',xaxt='n', col='white', ylim=genLim, main='generation times', ylab='values / ranges', xlab='time')
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$Tgl, MVSE_results$Tg$Tgu, colIn=alpha('green4',0.4), colLines=NA)
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$HTgl, MVSE_results$Tg$HTgu, colIn=alpha('blue4',0.4), colLines=NA)
    .plotbetweenY(MVSE_data$date, MVSE_results$Tg$VTgl, MVSE_results$Tg$VTgu, colIn=alpha('magenta',0.4), colLines=NA)
    lines(MVSE_data$date, MVSE_results$Tg$VTg, t='l', col='magenta')
    lines(MVSE_data$date, MVSE_results$Tg$HTg, t='l', col='blue4')
    lines(MVSE_data$date, MVSE_results$Tg$Tg, t='l', col='green4')
    labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
    axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
    legend('bottomleft',legend=c('total mean & 95% CI','human contribution mean & 95% CI','vector contribution mean & 95% CI'),col=c('green4','blue4','magenta'), lwd=2, bg='white')
    box(lwd=1.5)

    plot(MVSE_data$date, 1-MVSE_results$Tg$FrNoVecCap, t='l',xaxt='n', ylim=c(0,1), main='capacity (samples with vector lifespan > incubation period)', ylab='proportion', xlab='time')
    labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
    axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
    box(lwd=1.5)
  a<- dev.off()
}

#' Plots and exports a simple figure (mostly) for debug purposes. The distributions
#' of factors Alpha, Eta and Rho are plotted using R's default \code{density()}
#' function.
#'
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEcoCoefPosteriors<- function(outfilename='debug_dist', etaLim=c(), alphaLim=c(), rhoLim=c()){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=10*0.5,h=6*0.5,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)


  Ya<- density(MVSE_results$alpha, from=0, to=max(MVSE_results$alpha))
  Ye<- density(MVSE_results$eta, from=0, to=max(MVSE_results$eta))
  Yr<- density(MVSE_results$rho, from=0, to=max(MVSE_results$rho))
  Ya$y<- Ya$y/Ya$n
  Ye$y<- Ye$y/Ye$n
  Yr$y<- Yr$y/Yr$n

  if(length(etaLim)==0){
    etaLim<- c(min(Ya$x,Ye$x,Yr$x),max(Ya$x,Ye$x,Yr$x)*1.1)
  }
  if(length(alphaLim)==0){
    alphaLim<- c(min(Ya$x,Ye$x,Yr$x),max(Ya$x,Ye$x,Yr$x)*1.1)
  }
  if(length(rhoLim)==0){
    rhoLim<- c(min(Ya$x,Ye$x,Yr$x),max(Ya$x,Ye$x,Yr$x)*1.1)
  }

  plot(Ya$x, Ya$y, t='l', col='red3', lwd=2, main='posterior for eco-coefficient Alpha', xlab='value',ylab='frequency', xlim=alphaLim, ylim=c(min(Ya$y),max(Ya$y)*1.3) )
  .plotbellow(Ya$x, Ya$y, colIn='red3', colLines='red3')
  box(lwd=1.5)

  plot(Ye$x, Ye$y, t='l',col='blue', lwd=2, main='posterior for eco-coefficient Eta', xlab='value',ylab='frequency', xlim=etaLim, ylim=c(min(Ye$y),max(Ye$y)*1.3) )
  .plotbellow(Ye$x, Ye$y, colIn='blue', colLines='blue')
  box(lwd=1.5)

  plot(Yr$x, Yr$y, t='l',col='green3', lwd=2, main='posterior for eco-coefficient Rho', xlab='value',ylab='frequency', xlim=rhoLim, ylim=c(min(Yr$y),max(Yr$y)*1.3) )
  .plotbellow(Yr$x, Yr$y, colIn='green3', colLines='green3')
  box(lwd=1.5)

  Ya$y<- Ya$y/max(Ya$y)
  Ye$y<- Ye$y/max(Ye$y)
  Yr$y<- Yr$y/max(Yr$y)
  plot(Ya$x, Ya$y, t='l', col='white', lwd=2, main='posteriors for eco-coefficients', xlab='value',ylab='frequency', xlim=alphaLim, ylim=c(min(Ya$y),max(Ya$y)*1.3) )
  .plotbellow(Ya$x, Ya$y, colIn=alpha('red4',0.33), colLines='red2')
  .plotbellow(Ye$x, Ye$y, colIn=alpha('blue3',0.33), colLines='blue3')
  .plotbellow(Yr$x, Yr$y, colIn=alpha('green3',0.33), colLines='green3')
  lines(Ya$x, Ya$y, t='l', col='red2', lwd=3, lty=2)
  lines(Yr$x, Yr$y, t='l',col='green3', lwd=3, lty=2)
  lines(Ye$x, Ye$y, t='l',col='blue3', lwd=3, lty=2)
  legend('topright', col=c('red4','green3','blue'), legend=c('alpha','rho','eta'), lwd=3, cex=1.5)
  box(lwd=1.5)

  a<-dev.off()
}

#' Plots the response, in terms of theoretical index P, into a bubble chart with climated
#' as input variables and months identified in the range of the climatic vars.
#'
#' @param outfilename A name tag for the output file.
#' @param MinVarNorm scale colours to minimum desired value
#' @param MaxVarNorm scale colours to maximum desired value
#' @param TempLim Two numbers for the limit of temperature to be plotted.
#' @param HumLim Two numbers for the limit of humidity to be plotted.
#' @param cols Colours to be used in the colour key.
#'
#' @return NULL
#'
plotTheoreticalSuitRespMap<- function(outfilename='debug_theoretical_response', MinVarNorm=NA, MaxVarNorm=NA, TempLim=c(), HumLim=c(), cols= rev(c('red3','gold','skyblue','violet','purple'))){

  if(is.na(MaxVarNorm) | is.na(MinVarNorm) ){
    print("No MaxVarNorm and / or MinVarNorm given, adjusting to minimum / maximum in data set.")
    MaxVarNorm<- max(MVSE_results$indexP_theoretical)*1.001
    MinVarNorm<- min(MVSE_results$indexP_theoretical)*0.999
  }

  if(length(TempLim)!=2){
    print("No TempLim given, adjusting to maximum in data set.")
    TempLim<- range(MVSE_data_theoretical$T)
  }

  if(length(HumLim)!=2){
    print("No HumLim given, adjusting to maximum in data set.")
    HumLim<- range(MVSE_data_theoretical$H)*100
  }

  if(MaxVarNorm< max(MVSE_results$indexP_theoretical)){
    stop(paste("The parameter MaxVarNorm of this function has to be >= than the max in the data - which is ",max(MVSE_results$indexP_theoretical)))
  }

  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=7,bg='white')

  ff<- 1.3
  dta<- data.frame(T=MVSE_data_theoretical$T, H=MVSE_data_theoretical$oH, P=MVSE_results$indexP_theoretical)
  dta_emp<- data.frame(T=MVSE_data$T, H=MVSE_data$oH, P=MVSE_results$indexP)

  smored3P<- round(dta$P,5)
  smored3P_emp<- round(dta_emp$P,5)

  layout(matrix(1:2,ncol=2), width = c(0.82,0.18), height = c(1,1))
  par(mar=c(4.5, 4.5, 2, 0.5))

  minP<- min(smored3P)
  maxP<- max(smored3P)

  minP_emp<- min(smored3P_emp)
  meanP_emp<- mean(smored3P_emp)
  maxP_emp<- max(smored3P_emp)

  f <- colorRamp(cols, bias=1.1, interpolate = c("linear"))
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  rr<- range01(c(smored3P,smored3P_emp,MinVarNorm,MaxVarNorm))
  rr<- rr[1:(length(rr)-2)]
  colsPlotAll <- rgb(f(rr)/255)
  colsPlot<- head(colsPlotAll, length(smored3P))
  colsPlot_emp<- tail(colsPlotAll, length(smored3P_emp))

  plot(dta$H, dta$T, cex=ff, col=colsPlot, pch=20, xlim=HumLim, ylim=TempLim, xlab='humidity', ylab='temperature (celsius)')
  points(dta_emp$H, dta_emp$T, cex=ff, col='grey22', pch=21)
  points(dta_emp$H, dta_emp$T, cex=ff, col=colsPlot_emp, pch=20)
  box(lwd=1.5)

  par(mar=c(4.5, 1, 2, 2))
  cbr<- 50
  tbr<- 10
  legend_image <- as.raster(matrix( rev( rgb(f(seq(0,1,l=cbr))/255) ) , ncol=1))
  plot(c(0,2),c(0,MaxVarNorm),type = 'n', axes = F,xlab = '', ylab = '', main = 'index P')
  text(x=1.6, y=seq(0,MaxVarNorm,l=tbr), labels=round(seq(MinVarNorm,MaxVarNorm,l=tbr),2), font=2, cex=0.8)
  rasterImage(legend_image, 0, 0, 1, MaxVarNorm, interpolate=FALSE)

  a<-dev.off()
}


#' Plots the response, in terms of empirical generation time, into a bubble chart with climated
#' as input variables and months identified in the range of the climatic vars.
#'
#' @param outfilename A name tag for the output file.
#' @param MinVarNorm scale colours to minimum desired value
#' @param MaxVarNorm scale colours to maximum desired value
#' @param TempLim Two numbers for the limit of temperature to be plotted.
#' @param HumLim Two numbers for the limit of humidity to be plotted.
#' @param cols Colours to be used in the colour key.
#'
#' @return NULL
#'
plotEmpiricalGenTimeMap<- function(outfilename='debug_empirical_genTime', MinVarNorm=NA, MaxVarNorm=NA, TempLim=c(), HumLim=c(), cols= rev(c('red3','gold','skyblue','violet','purple'))){

  variable<- MVSE_results$Tg$Tg

  if(is.na(MaxVarNorm) | is.na(MinVarNorm) ){
    print("No MaxVarNorm and / or MinVarNorm given, adjusting to minimum / maximum in data set.")
    MaxVarNorm<- max(variable)*1.001
    MinVarNorm<- min(variable)*0.999
  }

  if(length(TempLim)!=2){
    print("No TempLim given, adjusting to maximum in data set.")
    TempLim<- range(MVSE_data$T)
  }

  if(length(HumLim)!=2){
    print("No HumLim given, adjusting to maximum in data set.")
    HumLim<- range(MVSE_data$H)*100
  }


  if(MaxVarNorm< max(variable)){
    stop(paste("The parameter MaxVarNorm of this function has to be >= than the max in the data - which is ",max(variable)))
  }

  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=7,bg='white')

  ff<- 1.3
  dta<- data.frame(T=MVSE_data$T, H=MVSE_data$oH, P=variable, M=as.numeric(format(as.Date(MVSE_data$date),"%m")))

  dta.m.T.m<- c()
  dta.m.H.m<- c()
  dta.m.T.r<- c()
  dta.m.H.r<- c()
  M<- 1:12
  for(mm in M){
    dta.m<- dta[which(dta$M==mm),]
    dta.m.T.m<- c(dta.m.T.m, median(dta.m$T))
    dta.m.H.m<- c(dta.m.H.m, median(dta.m$H))
  }

  smored3P<- round(variable,5)

  layout(matrix(1:2,ncol=2), width = c(0.82,0.18), height = c(1,1))
  par(mar=c(4.5, 4.5, 2, 0.5))

  minP<- min(smored3P)
  meanP<- mean(smored3P)
  maxP<- max(smored3P)

  f <- colorRamp(cols)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  rr<- range01(c(smored3P,MaxVarNorm,MinVarNorm))
  rr<- rr[1:(length(rr)-2)]
  colsPlot <- rgb(f(rr)/255)
  plot(dta$H, dta$T, cex=ff, col='grey22', pch=21, xlim=HumLim, ylim=TempLim, xlab='humidity', ylab='temperature (celsius)')
  points(dta$H, dta$T, cex=ff, col=colsPlot, pch=20)
  lines(c(dta.m.H.m,dta.m.H.m[1]), c(dta.m.T.m,dta.m.T.m[1]), lwd=5 , col='black')
  lines(c(dta.m.H.m,dta.m.H.m[1]), c(dta.m.T.m,dta.m.T.m[1]), lwd=1 , col='white')
  M<- 1:12

  points <- length(M)
    xdif<- diff(range(dta.m.H.m))
    ydif<- diff(range(dta.m.T.m))
    radius <- mean(xdif,ydif)
    center_x <- mean(dta.m.H.m)
    center_y <- mean(dta.m.T.m)
    xf<- (xdif/ydif)/ (max(xdif)/max(ydif))
    yf<- (ydif/xdif)/ (max(ydif)/max(ydif))
    drawCirclePoints <- function(points, radius, center_x, center_y, xf, yf) {
      slice <- 2 * pi / points
      ang<- atan2(dta.m.T.m[1] - center_y, dta.m.H.m[1] - center_x)
      angle <- slice * rev(seq(0, points, by = 1)) -ang
      angle <- seq(ang,ang-slice*points,-slice)
      newX <- center_x + radius * cos(angle)* xf
      newY <- center_y + radius * sin(angle)* yf
      for(mm in 1:points){
        lines(c(dta.m.H.m[mm],newX[mm]),c(dta.m.T.m[mm],newY[mm]))
        points(c(dta.m.H.m[mm],newX[mm]),c(dta.m.T.m[mm],newY[mm]), col='white', pch=20)
        points(newX[mm], newY[mm], pch=21, fg='black', bg='white',cex=3.5)
        text(newX[mm], newY[mm], mm, font=2, cex=0.9, col='black')
      }
    }
  drawCirclePoints(points, radius, center_x, center_y, xf, yf)
  box(lwd=1.5)
  par(mar=c(4.5, 1, 2, 2))
  cbr<- 50
  tbr<- 10
  legend_image <- as.raster(matrix( rev( rgb(f(seq(0,1,l=cbr))/255) ) , ncol=1))
  plot(c(0,2),c(0,MaxVarNorm),type = 'n', axes = F,xlab = '', ylab = '', main = 'gen. time')
  text(x=1.6, y=seq(0,MaxVarNorm,l=tbr), labels=round(seq(MinVarNorm,MaxVarNorm,l=tbr),2), font=2, cex=0.8)
  rasterImage(legend_image, 0, 0, 1, MaxVarNorm, interpolate=FALSE)

  a<-dev.off()
}

#' Plots the response, in terms of empirical capacity, into a bubble chart with climated
#' as input variables and months identified in the range of the climatic vars.
#'
#' Capacity is the expected proportion of human to mosquito to human transmission events
#' that would have mosquito lifespan longer than the mosquito incubation period. Capacity
#' towards zero implies unlikeliness that the mosquitoes will live longer enough to transmit.
#'
#' @param outfilename A name tag for the output file.
#' @param MinVarNorm scale colours to minimum desired value
#' @param MaxVarNorm scale colours to maximum desired value
#' @param TempLim Two numbers for the limit of temperature to be plotted.
#' @param HumLim Two numbers for the limit of humidity to be plotted.
#' @param cols Colours to be used in the colour key.
#'
#' @return NULL
#'
plotEmpiricalVecCapMap<- function(outfilename='debug_empirical_vecCap', MinVarNorm=NA, MaxVarNorm=NA, TempLim=c(), HumLim=c(), cols= rev(c('red3','gold','skyblue','violet','purple'))){

  variable<- 1- MVSE_results$Tg$FrNoVecCap

  if(is.na(MaxVarNorm) | is.na(MinVarNorm) ){
    print("No MaxVarNorm and / or MinVarNorm given, adjusting to minimum / maximum in data set.")
    MaxVarNorm<- max(variable)*1.001
    MinVarNorm<- min(variable)*0.999
  }

  if(length(TempLim)!=2){
    print("No TempLim given, adjusting to maximum in data set.")
    TempLim<- range(MVSE_data$T)
  }

  if(length(HumLim)!=2){
    print("No HumLim given, adjusting to maximum in data set.")
    HumLim<- range(MVSE_data$H)*100
  }

  if(MaxVarNorm< max(variable)){
    stop(paste("The parameter MaxVarNorm of this function has to be >= than the max in the data - which is ",max(variable)))
  }

  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=7,bg='white')

  ff<- 1.3
  dta<- data.frame(T=MVSE_data$T, H=MVSE_data$oH, P=variable, M=as.numeric(format(as.Date(MVSE_data$date),"%m")))

  dta.m.T.m<- c()
  dta.m.H.m<- c()
  dta.m.T.r<- c()
  dta.m.H.r<- c()
  M<- 1:12
  for(mm in M){
    dta.m<- dta[which(dta$M==mm),]
    dta.m.T.m<- c(dta.m.T.m, median(dta.m$T))
    dta.m.H.m<- c(dta.m.H.m, median(dta.m$H))
  }

  smored3P<- round(variable,5)

  layout(matrix(1:2,ncol=2), width = c(0.82,0.18), height = c(1,1))
  par(mar=c(4.5, 4.5, 2, 0.5))

  minP<- min(smored3P)
  meanP<- mean(smored3P)
  maxP<- max(smored3P)

  f <- colorRamp(cols)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  rr<- range01(c(smored3P,MaxVarNorm,MinVarNorm))
  rr<- rr[1:(length(rr)-2)]
  colsPlot <- rgb(f(rr)/255)
  plot(dta$H, dta$T, cex=ff, col='grey22', pch=21, xlim=HumLim, ylim=TempLim, xlab='humidity', ylab='temperature (celsius)')
  points(dta$H, dta$T, cex=ff, col=colsPlot, pch=20)
  lines(c(dta.m.H.m,dta.m.H.m[1]), c(dta.m.T.m,dta.m.T.m[1]), lwd=5 , col='black')
  lines(c(dta.m.H.m,dta.m.H.m[1]), c(dta.m.T.m,dta.m.T.m[1]), lwd=1 , col='white')
  M<- 1:12
  points <- length(M)
    xdif<- diff(range(dta.m.H.m))
    ydif<- diff(range(dta.m.T.m))
    radius <- mean(xdif,ydif)
    center_x <- mean(dta.m.H.m)
    center_y <- mean(dta.m.T.m)
    xf<- (xdif/ydif)/ (max(xdif)/max(ydif))
    yf<- (ydif/xdif)/ (max(ydif)/max(ydif))
    drawCirclePoints <- function(points, radius, center_x, center_y, xf, yf) {
      slice <- 2 * pi / points
      ang<- atan2(dta.m.T.m[1] - center_y, dta.m.H.m[1] - center_x)
      angle <- slice * rev(seq(0, points, by = 1)) -ang
      angle <- seq(ang,ang-slice*points,-slice)
      newX <- center_x + radius * cos(angle)* xf
      newY <- center_y + radius * sin(angle)* yf
      for(mm in 1:points){
        lines(c(dta.m.H.m[mm],newX[mm]),c(dta.m.T.m[mm],newY[mm]))
        points(c(dta.m.H.m[mm],newX[mm]),c(dta.m.T.m[mm],newY[mm]), col='white', pch=20)
        points(newX[mm], newY[mm], pch=21, fg='black', bg='white',cex=3.5)
        text(newX[mm], newY[mm], mm, font=2, cex=0.9, col='black')
      }
    }
  drawCirclePoints(points, radius, center_x, center_y, xf, yf)
  box(lwd=1.5)
  par(mar=c(4.5, 1, 2, 2))
  cbr<- 50
  tbr<- 10
  legend_image <- as.raster(matrix( rev( rgb(f(seq(0,1,l=cbr))/255) ) , ncol=1))
  plot(c(0,2),c(0,MaxVarNorm),type = 'n', axes = F,xlab = '', ylab = '', main = 'capacity')
  text(x=1.6, y=seq(0,MaxVarNorm,l=tbr), labels=round(seq(MinVarNorm,MaxVarNorm,l=tbr),2), font=2, cex=0.8)
  rasterImage(legend_image, 0, 0, 1, MaxVarNorm, interpolate=FALSE)

  a<-dev.off()
}




#' Plots the response, in terms of empirical index P, into a bubble chart with climated
#' as input variables and months identified in the range of the climatic vars.
#'
#' @param outfilename A name tag for the output file.
#' @param MinVarNorm scale colours to minimum desired value
#' @param MaxVarNorm scale colours to maximum desired value
#' @param TempLim Two numbers for the limit of temperature to be plotted.
#' @param HumLim Two numbers for the limit of humidity to be plotted.
#' @param cols Colours to be used in the colour key.
#'
#' @return NULL
#'
plotEmpiricalSuitRespMap<- function(outfilename='debug_empirical_response',  MinVarNorm=NA, MaxVarNorm=NA,  TempLim=c(), HumLim=c(), cols= rev(c('red3','gold','skyblue','violet','purple'))){

  if(is.na(MaxVarNorm) | is.na(MinVarNorm) ){
    print("No MaxVarNorm and / or MinVarNorm given, adjusting to minimum / maximum in data set.")
    MaxVarNorm<- max(MVSE_results$indexP)*1.001
    MinVarNorm<- min(MVSE_results$indexP)*0.999
  }

  if(length(TempLim)!=2){
    print("No TempLim given, adjusting to maximum in data set.")
    TempLim<- range(MVSE_data$T)
  }

  if(length(HumLim)!=2){
    print("No HumLim given, adjusting to maximum in data set.")
    HumLim<- range(MVSE_data$H)*100
  }

  if(MaxVarNorm< max(MVSE_results$indexP)){
    stop(paste("The parameter MaxVarNorm of this function has to be >= than the max in the data - which is ",max(MVSE_results$indexP)))
  }

  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=7,bg='white')

  ff<- 1.3
  dta<- data.frame(T=MVSE_data$T, H=MVSE_data$oH, P=MVSE_results$indexP, M=as.numeric(format(as.Date(MVSE_data$date),"%m")))

  dta.m.T.m<- c()
  dta.m.H.m<- c()
  dta.m.T.r<- c()
  dta.m.H.r<- c()
  M<- 1:12
  for(mm in M){
    dta.m<- dta[which(dta$M==mm),]
    dta.m.T.m<- c(dta.m.T.m, median(dta.m$T))
    dta.m.H.m<- c(dta.m.H.m, median(dta.m$H))
  }

  smored3P<- round(MVSE_results$indexP,5)

  layout(matrix(1:2,ncol=2), width = c(0.82,0.18), height = c(1,1))
  par(mar=c(4.5, 4.5, 2, 0.5))
  minP<- min(smored3P)
  meanP<- mean(smored3P)
  maxP<- max(smored3P)
  f <- colorRamp(cols)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  rr<- range01(c(smored3P, MinVarNorm, MaxVarNorm))
  rr<- rr[1:(length(rr)-2)]
  colsPlot <- rgb(f(rr)/255) #adding MaxVarNorm makes sure it respects maximum wanted
  plot(dta$H, dta$T, cex=ff, col='grey22', pch=21, xlim=HumLim, ylim=TempLim, xlab='humidity', ylab='temperature (celsius)')
  points(dta$H, dta$T, cex=ff, col=colsPlot, pch=20)
  lines(c(dta.m.H.m,dta.m.H.m[1]), c(dta.m.T.m,dta.m.T.m[1]), lwd=5 , col='black')
  lines(c(dta.m.H.m,dta.m.H.m[1]), c(dta.m.T.m,dta.m.T.m[1]), lwd=1 , col='white')
  M<- 1:12

  points <- length(M)
    xdif<- diff(range(dta.m.H.m))
    ydif<- diff(range(dta.m.T.m))
    radius <- mean(xdif,ydif)
    center_x <- mean(dta.m.H.m)
    center_y <- mean(dta.m.T.m)
    xf<- (xdif/ydif)/ (max(xdif)/max(ydif))
    yf<- (ydif/xdif)/ (max(ydif)/max(ydif))
    drawCirclePoints <- function(points, radius, center_x, center_y, xf, yf) {
      slice <- 2 * pi / points
      ang<- atan2(dta.m.T.m[1] - center_y, dta.m.H.m[1] - center_x)
      angle <- slice * rev(seq(0, points, by = 1)) -ang
      angle <- seq(ang,ang-slice*points,-slice)
      newX <- center_x + radius * cos(angle)* xf
      newY <- center_y + radius * sin(angle)* yf
      for(mm in 1:points){
        lines(c(dta.m.H.m[mm],newX[mm]),c(dta.m.T.m[mm],newY[mm]))
        points(c(dta.m.H.m[mm],newX[mm]),c(dta.m.T.m[mm],newY[mm]), col='white', pch=20)
        points(newX[mm], newY[mm], pch=21, fg='black', bg='white',cex=3.5)
        text(newX[mm], newY[mm], mm, font=2, cex=0.9, col='black')
      }
    }
  drawCirclePoints(points, radius, center_x, center_y, xf, yf)
  box(lwd=1.5)
  par(mar=c(4.5, 1, 2, 2))
  cbr<- 50
  tbr<- 10
  legend_image <- as.raster(matrix( rev( rgb(f(seq(0,1,l=cbr))/255) ) , ncol=1))
  plot(c(0,2),c(0,MaxVarNorm),type = 'n', axes = F,xlab = '', ylab = '', main = 'index P')
  text(x=1.6, y=seq(0,MaxVarNorm,l=tbr), labels=round(seq(MinVarNorm,MaxVarNorm,l=tbr),2), font=2, cex=0.8)
  rasterImage(legend_image, 0, 0, 1, MaxVarNorm, interpolate=FALSE)
  # abline(h=MaxVarNorm)
  # abline(h=max(MVSE_results$indexP))
  # abline(h=min(MVSE_results$indexP))
  a<-dev.off()
}


#' Plots and exports a simple figure (mostly) for debug purposes. The mean index
#' P and 95\% CI are plotted.
#'
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotClimate<- function(outfilename='debug_climate'){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=10*0.5,h=6*0.5,bg='white')

  par(mar=c(4, 4, 1.5, 5),cex=0.9)
  MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

  plot(MVSE_data$date,MVSE_data$oH,t='l',xaxt='n', ylab='humidity', xlab='date', main='climate', col='magenta', lwd=2, ylim=c(20,100))
  par(new=TRUE)
  plot(MVSE_data$date,MVSE_data$T,t='l',col="cadetblue4",xaxt="n",yaxt="n",xlab="",lwd=2,ylab="", ylim=c(10,40))
  axis(4)
  mtext("temperature",side=4,line=3)
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  legend('topleft',legend=c('humidity','temperature'),col=c('magenta','cadetblue4'), lwd=2, cex=0.8)
  box(lwd=1.5)

  plot(MVSE_data$date,MVSE_data$R,t='l',xaxt='n', ylab='rainfall', xlab='date', main='climate', col='limegreen', lwd=2, ylim=c(0,1))
  par(new=TRUE)
  plot(MVSE_data$date,MVSE_data$T,t='l',col="cadetblue4",xaxt="n",yaxt="n",xlab="",lwd=2,ylab="", ylim=c(10,40))
  axis(4)
  mtext("temperature",side=4,line=3)
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  legend('topleft',legend=c('rainfall','temperature'),col=c('limegreen','cadetblue4'), lwd=2, cex=0.8)
  box(lwd=1.5)

  a<-dev.off()
}

#' Plots and exports the mean Q
#' and 95\% CI under smoothing and no smoothing.
#'
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEmpiricalQ<- function(outfilename='debug_Q'){
      outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
      pdf(outfilename,w=8,h=6,bg='white')
      par(mar=c(4, 4, 1.5, 4),cex=0.9)
      MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

      plot(MVSE_data$date,MVSE_results$Q,t='l',xaxt='n', ylab='values', xlab='date', main='estimated Q', col='skyblue', lwd=3, ylim=c(0,max(MVSE_results$Qu)*1.1))
      lines(MVSE_data$date,MVSE_results$Ql,col='grey77')
      lines(MVSE_data$date,MVSE_results$Qu,col='grey77')
      lines(MVSE_data$date,MVSE_results$Q,col='skyblue')
      labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
      axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
      abline(h=1, col='black', lwd=3, lty=2)
      legend('topleft',legend=c('95% CI','mean'),col=c('grey77','skyblue'), lwd=2)
      box(lwd=1.5)

      if(length(MVSE_results$QSmooth)){
      for(ss in 1:length(MVSE_results$QSmooth)){
        mmain<- paste('estimated Q - ',names(MVSE_results$QSmooth)[ss])
        plot(MVSE_data$date,MVSE_results$Q,t='l',xaxt='n', ylab='values', xlab='date', main=mmain, col='skyblue', lwd=3, ylim=c(0,max(MVSE_results$Qu)*1.1))
        lines(MVSE_data$date,MVSE_results$Ql,col='grey77')
        lines(MVSE_data$date,MVSE_results$Qu,col='grey77')
        lines(MVSE_data$date,MVSE_results$Q,col='skyblue')
        lines(MVSE_data$date,MVSE_results$QSmooth[[ss]],col='black',lwd=3)
        labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
        axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
        abline(h=1, col='black', lwd=3, lty=2)
        legend('topleft',legend=c('95% CI','mean'),col=c('grey77','skyblue'), lwd=2)
        box(lwd=1.5)
      }}

      a<-dev.off()
}

#' Plots and exports the mean V0
#' and 95\% CI under smoothing and no smoothing.
#'
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEmpiricalV0<- function(outfilename='debug_V0'){
      outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
      pdf(outfilename,w=8,h=6,bg='white')
      par(mar=c(4, 4, 1.5, 4),cex=0.9)
      MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

      plot(MVSE_data$date,MVSE_results$V0,t='l',xaxt='n', ylab='values', xlab='date', main='estimated V0', col='limegreen', lwd=3, ylim=c(0,max(MVSE_results$V0u)*1.1))
      lines(MVSE_data$date,MVSE_results$V0l,col='grey77')
      lines(MVSE_data$date,MVSE_results$V0u,col='grey77')
      lines(MVSE_data$date,MVSE_results$V0,col='skyblue')
      labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
      axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
      abline(h=1, col='black', lwd=3, lty=2)
      legend('topleft',legend=c('95% CI','mean'),col=c('grey77','limegreen'), lwd=2)
      box(lwd=1.5)

      if(length(MVSE_results$V0Smooth)){
      for(ss in 1:length(MVSE_results$V0Smooth)){
        mmain<- paste('estimated Q - ',names(MVSE_results$V0Smooth)[ss])
        plot(MVSE_data$date,MVSE_results$V0,t='l',xaxt='n', ylab='values', xlab='date', main=mmain, col='limegreen', lwd=3, ylim=c(0,max(MVSE_results$V0u)*1.1))
        lines(MVSE_data$date,MVSE_results$V0l,col='grey77')
        lines(MVSE_data$date,MVSE_results$V0u,col='grey77')
        lines(MVSE_data$date,MVSE_results$V0,col='limegreen')
        lines(MVSE_data$date,MVSE_results$V0Smooth[[ss]],col='black',lwd=3)
        labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
        axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
        abline(h=1, col='black', lwd=3, lty=2)
        legend('topleft',legend=c('95% CI','mean'),col=c('grey77','limegreen'), lwd=2)
        box(lwd=1.5)
      }}
      a<-dev.off()
}


#' Plots and exports the mean index P, Q*P and V0*P
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEmpiricalIndexPPV0PQ<- function(outfilename='debug_mean_P_PV0_PQ'){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=6,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)
  MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

  plot(MVSE_data$date,MVSE_results$indexP/max(MVSE_results$indexP),t='l',xaxt='n', ylab='values (norm)', xlab='date', main='estimated index P', col='tomato3', lwd=3, ylim=c(0,1))
  lines(MVSE_data$date,(MVSE_results$indexP*MVSE_results$Q)/max(MVSE_results$indexP*MVSE_results$Q),col='skyblue')
  lines(MVSE_data$date,(MVSE_results$indexP*MVSE_results$V0)/max(MVSE_results$indexP*MVSE_results$V0),col='limegreen')
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  legend('bottomleft',legend=c('P','QxP','V0xP'),col=c('tomato3','skyblue','limegreen'), lwd=2)
  box(lwd=1.5)

  a<-dev.off()
}


#' Plots and exports the mean index P, Q and V0
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEmpiricalIndexPV0Q<- function(outfilename='debug_mean_P_Q_V0'){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=6,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)
  MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

  plot(MVSE_data$date,MVSE_results$indexP/max(MVSE_results$indexP),t='l',xaxt='n', ylab='values (norm)', xlab='date', main='estimated index P', col='tomato3', lwd=3, ylim=c(0,1))
  lines(MVSE_data$date,MVSE_results$V0/max(MVSE_results$V0),lwd=3,col='limegreen')
  lines(MVSE_data$date,MVSE_results$Q/max(MVSE_results$Q),lwd=3,col='skyblue')
  lines(MVSE_data$date,(MVSE_results$indexP*MVSE_results$Q)/max(MVSE_results$indexP*MVSE_results$Q),lwd=3,col='black')
  lines(MVSE_data$date,(MVSE_results$indexP*MVSE_results$V0)/max(MVSE_results$indexP*MVSE_results$V0),lwd=3,col='gold')
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  legend('bottomleft',legend=c('P','Q','V0','QxP','V0xP'),col=c('tomato3','skyblue','limegreen','black','gold'), lwd=2)
  box(lwd=1.5)

  a<-dev.off()
}


#' Plots and exports the mean index
#' P and 95\% CI under smoothing and no smoothing.
#'
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEmpiricalIndexP<- function(outfilename='debug_indexP'){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=8,h=6,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)
  MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

  plot(MVSE_data$date,MVSE_results$indexP,t='l',xaxt='n', ylab='values', xlab='date', main='estimated index P', col='tomato3', lwd=3, ylim=c(0,3.5))
  lines(MVSE_data$date,MVSE_results$indexPl,col='grey77')
  lines(MVSE_data$date,MVSE_results$indexPu,col='grey77')
  lines(MVSE_data$date,MVSE_results$indexP,col='tomato3')
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  abline(h=1, col='black', lwd=3, lty=2)
  legend('topleft',legend=c('95% CI','mean'),col=c('grey77','tomato3'), lwd=2)
  box(lwd=1.5)

  if(length(MVSE_results$indexPSmooth)){
  for(ss in 1:length(MVSE_results$indexPSmooth)){
    mmain<- paste('estimated index P - ',names(MVSE_results$indexPSmooth)[ss])
    plot(MVSE_data$date,MVSE_results$indexP,t='l',xaxt='n', ylab='values', xlab='date', main=mmain, col='tomato3', lwd=3, ylim=c(0,3.5))
    lines(MVSE_data$date,MVSE_results$indexPl,col='grey77')
    lines(MVSE_data$date,MVSE_results$indexPu,col='grey77')
    lines(MVSE_data$date,MVSE_results$indexP,col='tomato3')
    lines(MVSE_data$date,MVSE_results$indexPSmooth[[ss]],col='black',lwd=3)
    labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
    axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
    abline(h=1, col='black', lwd=3, lty=2)
    legend('topleft',legend=c('95% CI','mean'),col=c('grey77','tomato3'), lwd=2)
    box(lwd=1.5)
  }}

  a<-dev.off()
}


#' Plots and exports a simple figure (mostly) for debug purposes. The entomological
#' parameters that depend on the factors Alpha, Eta and Rho are plotted. These
#' parameters are the mosquito life-span (1/death rate), the mosquito incubation
#' period (1/rate incubation to infection), mosquito biting rate. Both the
#' resulting time series and distributions are presented. The distributions should
#' mimic the priors used before estimating the factors Alpha, Eta and Rho. The
#' priors are simulated when plotted and also presented.
#'
#' @param outfilename A name tag for the output file.
#' @param Ns Number of samples wanted for posteriors.
#' @param entoPostLim Two numbers for the y limits to be used in the posterior plots.
#' @param bitPostLim Two numbers for the y limits to be used in the posterior plots.
#'
#' @return NULL
#'
plotEntoParameters<- function(outfilename='debug_ento', Ns=200, entoPostLim=c(), bitPostLim=c()){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=10*0.5,h=6*0.5,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9,cex.main=0.8)
  MVSE_data$date<- as.Date(MVSE_data$date, format='%Y-%m-%d')

  ttt<-length(MVSE_data$date)
  res<- .calculatePostsVLS_EIP_IIP_HIP(Ns, ttt)
  VLS<- res[[1]]
  EIP<- res[[2]]
  VLS.mean<- apply(VLS, MARG=2, FUN=mean)
  EIP.mean<- apply(EIP, MARG=2, FUN=mean)

  plot(MVSE_data$date,VLS.mean, main="(mean) mosquito life-span and incubation period",xaxt='n', t='l', lwd=2, ylab="value in days", xlab="time", ylim=c(0,max(EIP.mean,VLS.mean)))
  lines(MVSE_data$date,EIP.mean, col='purple', lwd=2); box(lwd=1.5)
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  legend('bottomright',legend=c('mosquito lifespan','incubation period'),col=c('black','purple'),lwd=2,lty=c(1,1),cex=0.8)
  box(lwd=1.5)

  Y1<- density(EIP.mean, from=0, to=30)
  Y2<- density(VLS.mean, from=0, to=30)
  priorGam<- density(MVSE_prior_ic_rdist(N=20000))
  priorMu<- density(MVSE_prior_lev_rdist(N=20000))
  priorA<- density(MVSE_prior_a_rdist(N=20000))

  if(length(entoPostLim)==0){
      entoPostLim<- c(min(Y1$x,Y2$x,priorGam$x,priorMu$x),max(Y1$x,Y2$x,priorGam$x,priorMu$x))
  }

  if(length(bitPostLim)==0){
      bitPostLim<- c(min(Y1$x,priorA$x),max(Y1$x,priorA$x))
  }

  Y1$y<- Y1$y/max(Y1$y)
  Y2$y<- Y2$y/max(Y2$y)
  priorGam$y<- priorGam$y/max(priorGam$y)
  priorMu$y<- priorMu$y/max(priorMu$y)
  plot(Y1$x, Y1$y, t='l', lwd=2, main="distributions for mosquito life-span and incubation period", col='purple', ylab="density", xlab="value in days", xlim=entoPostLim, ylim=c(min(Y1$y,Y2$y,priorGam$y,priorMu$y),max(Y1$y,Y2$y,priorGam$y,priorMu$y)))
    .plotbellow(Y1$x, Y1$y, colIn='purple', colLines='purple')
    .plotbellow(Y2$x, Y2$y, colIn='black', colLines='black')
        lines(Y1$x, Y1$y, t='l', lwd=2, col='beige')
        lines(Y1$x, Y1$y, t='l', lwd=2, col='purple')
        lines(Y2$x, Y2$y, t='l', lwd=2, col='beige')
        lines(Y2$x, Y2$y, t='l', lwd=2, col='black')
        lines(priorGam$x,priorGam$y, lwd=2, lty=1, col='beige')
        lines(priorGam$x,priorGam$y, lwd=2, lty=2, col='purple')
        lines(priorMu$x,priorMu$y, lwd=2, lty=1, col='beige')
        lines(priorMu$x,priorMu$y, lwd=2, lty=2, col='black')
  box(lwd=1.5)
  legend('topright',legend=c('prior','estimated life-span','estimated inc. period'),col=c('grey','black','purple'),lwd=2,lty=c(2,1,1),cex=0.8)

  A<- .calculatePostA(Ns, ttt)
  A.mean<- apply(A, MARG=2, FUN=mean)
  plot(MVSE_data$date, A.mean, main="mosquito biting rate", xaxt='n', t='l', ylim=c(min(A.mean)*0.9,max(A.mean)*1.1), lwd=2, col="orange", ylab="value", xlab="time")
  labDates<- seq(MVSE_data$date[1], tail(MVSE_data$date, 1), by = "months")
  axis.Date(side = 1, MVSE_data$date, at = labDates, format = "%b %y", las = 1)
  box(lwd=1.5)

  Y1<- density(A.mean, from=0, to=1)
  Y1$y<- Y1$y/max(Y1$y)
  priorA$y<- priorA$y/max(priorA$y)
  plot(Y1$x, Y1$y, t='l', lwd=2, main="distribution for biting rate", col="orange", ylab="density", xlab="value", xlim=bitPostLim, ylim=c(min(Y1$y,priorA$y),max(Y1$y,priorA$y)))
  .plotbellow(Y1$x, Y1$y, colIn='orange', colLines='orange')
  lines(priorA$x, priorA$y, lwd=2, lty=1, col='beige')
  lines(priorA$x, priorA$y, lwd=2, lty=2, col='orange')
  box(lwd=1.5)
  legend('topright',legend=c('estimated biting rate','prior'),col=c('orange','grey'),lwd=2,lty=c(1,2),cex=0.8)

  a<-dev.off()
}


#' Plots and exports a simple figure (mostly) for debug purposes. The MCMC chains
#' are presented for the factors Eta and Rho.
#'
#' @param outfilename A name tag for the output file.
#'
#' @return NULL
#'
plotEcoCoefMCMCChains<- function(outfilename='debug_chains'){
  outfilename<- paste0(.MVSEbuildOutPath(),'_',outfilename,'.pdf')
  pdf(outfilename,w=7,h=4,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)

  plot((MVSE_results$rho), ylim=c(0,max(MVSE_results$rho)*1.3), col='green3', lwd=2, main='MCMC chain for eco-coefficient Rho after burnin', xlab='MCMC step',ylab='value', t='l')
  abline(h=mean(MVSE_results$rho), lty=1, lwd=2)
  abline(h=quantile(MVSE_results$rho, probs=c(0.025,0.975)), lty=2, lwd=2)
  box(lwd=1.5)
  legend('topright',legend=c('mean and 95% CI','chain'),col=c('black','green3'),lwd=2)

  plot((MVSE_results$eta), ylim=c(0,max(MVSE_results$eta)*1.3), col='blue', lwd=2, main='MCMC chain for eco-coefficient Eta after burnin', xlab='MCMC step',ylab='value', t='l')
  abline(h=mean(MVSE_results$eta), lty=1, lwd=2)
  abline(h=quantile(MVSE_results$eta, probs=c(0.025,0.975)), lty=2, lwd=2)
  box(lwd=1.5)
  legend('topright',legend=c('mean and 95% CI','chain'),col=c('black','blue'),lwd=2)
  a<-dev.off()
}

################################


#' Choose samples from the posterior and fit filter.
#'
#' @param nSamples An integer defining the number of samples to apply filtering to. Default is all the samples in the input.
#' @param breakyears A logical variable indicating weather the data should be broken to years and analyzed separately.
#' @param CIlevel A numeric value between 0 and 1, indicating the level of confidence for generated confidence intervals (CI).
#' @return Returns a list with the average filter prediction, CI bounds of the average filter values,
#'         CI bounds for all filter values, and average index P values. This list is also saved as csv files.
#'         Additionally, ir produces a pdf plot of the filter with overlayed average index P.
#'
#'
expectedPosterior<-function(nSamples=NA,breakyears=FALSE,CIlevel=0.95){

  thisIndexP<- MVSE_indexP
  thisDates<- as.Date(MVSE_data$date)

  nSamples<- dim(thisIndexP)[1]

  ## taking a certain column from a dataframe/matrix
  takecol<- function(x,column){
    return(x[,column])
  }

  # t-based CI for the average fit
  cierr<- function(dat,conf){
    return(qt(1-conf/2, df=length(dat)-1)*sd(dat)/sqrt(length(dat)))
  }

  bestfit=list()
  averageP=list()
  if(breakyears==TRUE){
    breakdate<- thisDates[1]
    breakdatem<- strftime(breakdate,'%m')
    breakdated<- strftime(breakdate,'%d')
    Yearinds<- which(strftime(thisDates,'%m')==breakdatem & strftime(thisDates,'%d')==breakdated )
    numYears<- length(Yearinds)
    rowsout<- max(diff(Yearinds))
    splitvec<- rep(1:length(Yearinds),times=diff(c(Yearinds,length(thisDates)+1)))
  }

  for(i in 1:nSamples){
    if(length(dim(thisIndexP))>0){
      values<- thisIndexP[i,]
    } else {
      values<- thisIndexP }
    if (breakyears==TRUE){
      valueslist<- split(values,splitvec)
      ## determining the number of years using the timescalein paramter
      bestfit[[i]]<- matrix(ncol=numYears,nrow=rowsout)
      averageP[[i]]<- matrix(ncol=numYears,nrow=rowsout)

      for(j in 1:numYears){
        lvaluescurr<- valueslist[[j]]
        predfiltse<- .getFit(lvaluescurr)
        if(length(predfiltse$fit)!=length(bestfit[[i]][,j])){
          length(predfiltse$fit)<- length(bestfit[[i]][,j])
          length(lvaluescurr)<- length(bestfit[[i]][,j])
        }
        bestfit[[i]][,j]<- predfiltse$fit
        averageP[[i]][,j]<- lvaluescurr
      }
    } else{
      lvaluescurr<- values
      predfiltse<- .getFit(lvaluescurr)
      bestfit[[i]]<- predfiltse$fit
      averageP[[i]]<- lvaluescurr
    }
  }

  if (breakyears==TRUE){
    pdf('ExpectedPosteriorBroken.pdf',w=10,h=6,bg='white')

    if(numYears>2){
      par(mfrow=c(ceiling(sqrt(numYears)),ceiling(sqrt(numYears))))
    } else if(numYears==2){
      par(mfrow=c(1,2))
    }
    maxy<- max(colMeans(thisIndexP[1:nSamples,],na.rm=T))
    miny<- min(colMeans(thisIndexP[1:nSamples,],na.rm=T))

    for(i in 1:numYears){
      result<- matrix(unlist(lapply(bestfit,takecol,i)),ncol=nSamples,nrow=rowsout)
      result<- na.omit(result)
      resultAverage<- matrix(unlist(lapply(averageP,takecol,i)),ncol=nSamples,nrow=rowsout)
      resultAverage<- na.omit(resultAverage)

      currdates<- thisDates[Yearinds[i]:(c(Yearinds,length(thisDates)+1)[i+1]-1)]

      meanFit<- rowMeans(result,na.rm = T)
      meanFitCILow<- rowMeans(result,na.rm = T)-apply(result,1,cierr,CIlevel)
      meanFitCIHigh<- rowMeans(result,na.rm = T)+apply(result,1,cierr,CIlevel)
      AllFitCILow<- apply(result,1,quantile,(1-CIlevel)/2,na.rm=TRUE)
      AllFitCIHigh<- apply(result,1,quantile,1-(1-CIlevel)/2,na.rm=TRUE)
      Average<- rowMeans(resultAverage,na.rm = T)

      par(mar=c(4, 4, 1.5, 4),cex=0.9)
      epiyear<- paste(unique(strftime(currdates,'%Y')),collapse = '-')
      plot(currdates,meanFit, type='l',ylab='Index P',xlab='',xaxt="n",ylim=c(miny,maxy), main = paste('Epi year ',epiyear))
      axis.Date(1,currdates,format= "%y-%m", las = 2,xpd=TRUE,srt=45)
      box(lwd=2.5)
      lines(currdates,Average,col='tomato3')
      box(lwd=1.5)
      outdat<- data.frame(cbind(meanFit,meanFitCILow,meanFitCIHigh,AllFitCILow,AllFitCIHigh,Average))
      write.table(outdat,paste0(MVSE_output_tag,'Year ',as.character(i),'Smooth.csv'),col.names=TRUE,row.names=FALSE,sep=',')
    }
    dev.off()
    par(mfrow=c(1,1))
  } else {

    result<- matrix(unlist(bestfit),ncol=nSamples,byrow = FALSE)
    resultAverage<- matrix(unlist(averageP),ncol=nSamples,byrow = FALSE)
    meanFit<- rowMeans(result,na.rm = T)
    meanFitCILow<- rowMeans(result,na.rm = T)-apply(result,1,cierr,CIlevel)
    meanFitCIHigh<- rowMeans(result,na.rm = T)+apply(result,1,cierr,CIlevel)
    AllFitCILow<- apply(result,1,quantile,(1-CIlevel)/2,na.rm=TRUE)
    AllFitCIHigh<- apply(result,1,quantile,1-(1-CIlevel)/2,na.rm=TRUE)
    Average<- rowMeans(resultAverage,na.rm = T)

    miny<- min(Average,na.rm = T)
    maxy<- max(Average,na.rm = T)
    pdf('ExpectedPosteriorAll.pdf',w=10,h=6.5,bg='white')
    par(mar=c(4, 4, 1.5, 4),cex=0.9)
    plot(thisDates,meanFit, type='l',ylab='Index P',xlab='',xaxt="n",ylim=c(miny,maxy), main ='All years')
    axis.Date(1,thisDates,format= "%y-%m", las = 2,xpd=TRUE,srt=45)
    box(lwd=2.5)
    lines(thisDates,Average,col='tomato3')
    box(lwd=1.5)
    legend('topleft',legend=c('Filtered P','Averaged P'),fill=c('black','tomato3'))
    dev.off()
    outdat<- data.frame(cbind(meanFit,meanFitCILow,meanFitCIHigh,AllFitCILow,AllFitCIHigh,Average))
    write.table(outdat,paste0(MVSE_output_tag,'AllYearsSmooth.csv'),col.names=TRUE,row.names=FALSE,sep=',')
  }

  return(bestfit)
}


#' Choose samples from the posterior and find index P peaks by years.
#'
#' @param nSamples An integer defining the number of samples to apply the function on. Default is all the samples in the input.
#' @param CIlevel A numeric value between 0 and 1, indicating the level of confidence for generated confidence intervals (CI).
#' @return Returns a list with the max peak estimation, a low CI bound for the peak time, and
#'         a high CI bound for the peak time. This list is also saved as a csv file. Additionally, it produces a pdf barplot
#'         of peak times colored by year.
#'
#'
distributionPeak<-function(nSamples=NA, CIlevel=0.95){

  thisIndexP<- MVSE_indexP

  if(is.na(nSamples)) nSamples<-dim(thisIndexP)[1]

  allpeaks<-list()
  thisDates<-as.Date(MVSE_data$date)

  takecol<-function(x,column){
    return(x[,column])
  }

  averageP<-list()
  if(is.null(nSamples)){
    nSamples<- length(thisIndexP)
  }

  breakdate<- thisDates[1]
  breakdatem<- strftime(breakdate,'%m')
  breakdated<- strftime(breakdate,'%d')
  Yearinds<- which(strftime(thisDates,'%m')==breakdatem & strftime(thisDates,'%d')==breakdated )
  numYears<- length(Yearinds)
  rowsout<- max(diff(Yearinds))
  splitvec<- rep(1:length(Yearinds),times=diff(c(Yearinds,length(thisDates)+1)))

  for(i in 1:nSamples){
    if(length(dim(thisIndexP))>0){
      values<- thisIndexP[i,]
    } else {
      values<- thisIndexP
    }

    valueslist<- split(values,splitvec)
    averageP[[i]]<- matrix(ncol=numYears,nrow=rowsout)

    for(j in 1:numYears){
      lvaluescurr<- valueslist[[j]]
      if(length(lvaluescurr)!=length(averageP[[i]][,j])){
        length(lvaluescurr)<- length(averageP[[i]][,j])
      }
      averageP[[i]][,j]<- lvaluescurr
    }
  }

  finalmat<- (matrix(0,ncol=numYears,nrow = rowsout))

  for( i in 1:nSamples){
    Peaktime<- apply(averageP[[i]],2,which.max)
    for(j in 1:numYears){
      finalmat[Peaktime[j],j]=finalmat[Peaktime[j],j]+1
    }
  }

  finalmat<- t(finalmat)
  colnames(finalmat)<- as.character(1:dim(finalmat)[2])
  rownames(finalmat)<- paste('Year',as.character(1:dim(finalmat)[1]))

  cl<- colors(distinct=T)
  cls<- cl[sample(1:length(cl),numYears,replace = FALSE)]

  pdf('Peaks.pdf',w=10,h=6,bg='white')
  par(mar=c(4, 4, 1.5, 4),cex=0.9)
  alldates<- unique(strftime(thisDates,format='%m-%d'))
  barplot(names.arg = alldates[1:min(rowsout,length(alldates))],height = finalmat,xlim=c(0, ncol(finalmat) + 10),
          xlab='Peak time',ylab=paste0('Frequency (of ',as.character(nSamples),' simulations)'),
          col=cls,
          legend.text=TRUE,
          args.legend=list(x='top',
                           bty = "n")
          ,xpd = FALSE,offset = 0)
  box(lwd=1.5)
  dev.off()

  combYears<- matrix(NA,ncol=numYears,nrow=nSamples)
  for(i in 1:numYears){
    resultAverage<- matrix(unlist(lapply(averageP,takecol,i)),ncol=nSamples,nrow=rowsout)
    allpeaks[[i]]<-  .getPeak(resultAverage, CIlevel)
    finalvalues<- apply(resultAverage,2,which.max)
    finalvalues[sapply(finalvalues, function(x) length(x)==0)] <- NA
    finalvalues<- unlist(finalvalues)
    combYears[,i]<- finalvalues
  }
  combYears<- as.vector(combYears)

  low<- quantile(combYears,(1-CIlevel)/2,na.rm = TRUE)
  high<- quantile(combYears,1-(1-CIlevel)/2,na.rm = TRUE)
  Mean<- mean(combYears,na.rm = TRUE)
  allpeaks[[length(allpeaks)+1]]<- c(Mean,low,high)
  names(allpeaks)[1:(length(allpeaks)-1)]<- sapply(1:numYears,FUN=function(x){paste0('Year_',as.character(x))})
  names(allpeaks)[length(allpeaks)]='CombinedYears'

  peakdat<- data.frame(matrix(unlist(allpeaks),ncol=3,byrow =TRUE))
  names(peakdat)<- c('Average','CILow','CIHigh')
  peakdat$Year<- names(allpeaks)
  write.table(peakdat,paste0(MVSE_output_tag,'Peaks.csv'),col.names=TRUE,row.names=FALSE,sep=',')

  return(allpeaks)
}




#' Choose samples from the posterior and find the largest index P values.
#'
#' @param timethreshold A non-negative integer value indicating the maximum time units P can drop below 1
#'         while the interval of the 'suitable season' of P>=Pthreshold remains continuous.
#' @param Pthreshold A numeric value determining the threshold of P wanted to determine a suuitability season.
#' @return Returns a data frame with years as columns, time as rows, and 1 indicating  P>=Pthreshold .
#'        This data frame is also saved as a csv file. Additionally a pdf containing a plot of the peak times is produced.
#'
suitableSeason<-function(timethreshold=0,Pthreshold=1){

  thisIndexP<- MVSE_indexP
  thisDates<- MVSE_data$date

  dates<- as.Date(thisDates)

  breakdate<- dates[1]
  breakdatem<- strftime(breakdate,'%m')
  breakdated<- strftime(breakdate,'%d')
  Yearinds<- which(strftime(dates,'%m')==breakdatem & strftime(dates,'%d')==breakdated )
  numYears<- length(Yearinds)
  rowsout<- max(diff(Yearinds))

  splitvec<- rep(1:length(Yearinds),times=diff(c(Yearinds,length(dates)+1)))
  values<- thisIndexP
  valueslist<- split(values,splitvec)

  averageP<- matrix(ncol=numYears,nrow=rowsout)

  for(j in 1:numYears){

    lvaluescurr=valueslist[[j]]
    if(length(lvaluescurr)!=length(averageP[,j])){
      length(lvaluescurr)=length(averageP[,j])
    }

    averageP[,j]=lvaluescurr

  }

  allseasons=list()
  pdf('Seasons.pdf',w=10,h=6,bg='white')

  for ( i in 1:numYears){
    currdates<- dates[Yearinds[i]:(c(Yearinds,length(dates)+1)[i+1]-1)]

    allseasons[[i]]<- .getSeasons(averageP[,i], timethreshold,Pthreshold)

    if(i==1){
      par(mar=c(4, 4, 1.5, 4),cex=0.9)
      layout(matrix(c(2,2,1,1),2,2,byrow=T),widths=c(1,1), heights=c(4,6))
      plot(rep(i,length(allseasons[[i]]))~allseasons[[i]], pch=15,xlim=c(1,dim(averageP)[1]),ylim=c(-1,numYears),
           xlab='',xaxt="n",ylab='',yaxt="n")
      mylabels=seq.Date(min(currdates),max(currdates)+1,by='month')

      axis(side = 1, at = seq(from=1,to=length(currdates),length.out = 13), labels = FALSE)
      text(x =seq(from=1,to=length(currdates),length.out = 13) , y = -2, labels = format(mylabels, "%d/%m"),
           srt = 45, pos = 1, xpd = TRUE)

      epiyear=paste(unique(strftime(currdates,'%Y')),collapse = '-')
      axis(2, at=i, labels=FALSE)
      text(y=i, x=-2,labels=epiyear, srt=45, adj=1, xpd=TRUE)

      box(lwd=1.5)
      legend('bottom',legend=paste0('Dates with P>',as.character(Pthreshold)))


    } else {

      points(rep(i,length(allseasons[[i]]))~allseasons[[i]], pch=15)
      epiyear=paste(unique(strftime(currdates,'%Y')),collapse = '-')
      axis(2, at=i, labels=FALSE)
      text(y=i, x=-2,labels=epiyear, srt=45, adj=1, xpd=TRUE)

    }
  }

  seasondat<- data.frame(matrix(0,ncol=numYears,nrow = rowsout))

  for (i in 1:numYears){
    seasondat[allseasons[[i]],i]=1
  }

  seasondat[which(is.na(averageP[,dim(averageP)[2]])),dim(seasondat)[2]]=NA

  names(seasondat)=paste0('Year',as.character(1:numYears))
  write.table(seasondat,paste0(MVSE_output_tag,'Seasons.csv'),col.names=TRUE,row.names=FALSE,sep=',')

  par(mar=c(0,4,1,4),cex=0.9)
  # par(mar=c(4, 4, 1.5, 4),cex=0.9)
  barplot(rowMeans(seasondat,na.rm=T)*100,xlim=c(0.5,(rowsout)-0.5),space=0,ylim=c(0,100),
          ylab=paste0('% of years with P>',as.character(Pthreshold)))
  box(lwd=1.5)
  dev.off()
  return(seasondat)
}





#' Transform daily date into averaged weekly values.
#'
#' @param indexP A numeric matrix of Simulated index P, with rows as samples.
#' @param dates A character vector containing dates in a yyyy-mm-dd format. Should have the same length
#'        as the number of columns in indexP.
#' @return Returns a list with the the averaged index P as the first element and a vector of the corresponding
#'         dates as the second element. The data is broken by years so that every year starts at
#'         the same data as the first input date. This causes a loss of 1-2 days per year.
#'
weeklyAverage<-function(indexP,dates,nSamples=dim(indexP)[1]){

  dates=as.Date(dates)
  breakdate=dates[1]
  breakdatem=strftime(breakdate,'%m')
  breakdated=strftime(breakdate,'%d')
  Yearinds=which(strftime(dates,'%m')==breakdatem & strftime(dates,'%d')==breakdated )
  numYears=length(Yearinds)
  splitvec=rep(1:length(Yearinds),times=diff(c(Yearinds,length(dates)+1)))
  averageP=matrix(NA,nrow = nSamples,ncol=numYears*52)

  for(i in 1:nSamples){
    if(length(dim(indexP))>0){
      values=indexP[i,]
    } else {
      values=indexP }
    valueslist=split(values,splitvec)

    for(j in 1:numYears){

      lvaluescurr=valueslist[[j]]

      length(lvaluescurr)=min(7*ceiling(length(lvaluescurr)/7),364)

      valuesmat=matrix(lvaluescurr,nrow=7,byrow = FALSE)
      valuesweek=colMeans(valuesmat,na.rm=T)
      lvaluescurr=valuesweek
      lvaluescurr=na.omit(lvaluescurr)

      if(length(lvaluescurr)!=52){
        length(lvaluescurr)=52
      }

      averageP[i,(1+(j-1)*52):(j*52)]=lvaluescurr


    }

  }
  newdates=c()
  for(i in 1:numYears){
    newdates=c(newdates,seq.Date(dates[Yearinds[i]],by=7,length.out = 52))

  }
  newdates=as.Date(newdates,origin="1970-01-01")
  retlis=list()
  retlis[[1]]=averageP
  retlis[[2]]=newdates

  return(retlis)
}


###########################################
#### INTERNAL FUNCTIONS
###########################################

## WRAPPER TO CALL THE MCMC ON RHO AND ETA ESTIMATION, ALSO CLIPS OUTPUT
## ACCORDING TO BURNIN PERIOD
.estimateFactorsRhoEta<- function(nMCMC=10000, bMCMC=0.5, cRho=0.5, cEta=2.0, gauJump=0.05){
  accepted<- .MVSErunMCMC_RhoEta(nMCMC, cRho, cEta, gauJump)
  RHOSacc<- accepted[,1]
  ETASacc<- accepted[,2]
  RHOS<- RHOSacc[(length(RHOSacc)*bMCMC):length(RHOSacc)]
  ETAS<- ETASacc[(length(ETASacc)*bMCMC):length(ETASacc)]
  return(list(RHOS, ETAS))
}

##ESTIMATES ALPHA AND SAVES G IN MVSE_DATA (INCUBATION TO INFECTION RATE) IN GLOBAL
.estimateFactorAlpha<- function(N=10000){
  gammaV<- 1/MVSE_prior_ic_rdist(N)
  gammaV<- gammaV[which(gammaV>0 & gammaV<1)]
  alphas<- unlist(lapply(gammaV, FUN= function(X){ mean(X/getTempEffIncPer()) }))
  return(alphas)
}

##BUILDS AN OUTPUT PATH NAME THAT INCLUDES VALUES USED IN PRIORS BY DEFAULT
.MVSEbuildOutPath<- function(){
  # return(paste(MVSE_output_tag,MVSE_prior_lev_mean,MVSE_prior_lev_sd,MVSE_prior_ic_mean,MVSE_prior_ic_sd,MVSE_prior_a_mean,MVSE_prior_a_sd,sep='_'))
  return(paste(MVSE_output_tag,sep='_'))
}

##RUNS THE MCMC ON ETA AND RHO
.MVSErunMCMC_RhoEta<- function(nMCMC, cRho, cEta, gauJump){
  ##proposal variables
  pL<- NA
  pRho<- NA
  pEta<- NA
  ##
  cL<-0.00001
  ##
  accepted<- matrix(rep(NA,2*nMCMC),ncol=2)
  countAccepted<- 0
  ##
  print('Running MCMC...')
  pb <- txtProgressBar(min=0, max = nMCMC-1, style = 3)
  for(ii in 1:nMCMC){
    while(TRUE){
      pRho<- rnorm(1, mean=cRho, sd=gauJump)
      if(pRho>0 & pRho<10) break;
    }
    while(TRUE){
      pEta<- rnorm(1, mean=cEta, sd=gauJump)
      if(pEta>0 & pEta<10) break;
    }
    pMu<- 1/mean(pEta*getTempEffLifespan()*(1+getHumEffLifespan())^pRho)
    pA<- mean(MVSE_prior_a_mean*(1+MVSE_data$A)^pRho)
    ##address acceptance probability / likelihood
    p1<- MVSE_prior_lev_ddist(pMu)
    p2<- MVSE_prior_a_ddist(pA)
    #pL<- p1 + p2; pAccept<- exp(pL - cL) ##acceptance 'prob', log like
    pL<- p1 * p2; pAccept<- pL/cL ##acceptance 'prob'
    ##decision on step
    event<- runif(1, min=0, max=1)
    if(event< pAccept){
      ##accept
      accepted[ii,]<- c(pRho, pEta)
      cRho<- pRho
      cEta<- pEta
      cL<- pL
      countAccepted<- countAccepted+ 1
    }else{
      ##reject
      accepted[ii,]<- c(cRho, cEta)
    }
    setTxtProgressBar(pb, ii)
  }
  colnames(accepted)<- c('Rho','Eta')
  print(' ')
  print(paste('Acceptance %:',round(100*countAccepted/nMCMC,3)))
  close(pb)
  return(accepted)
}

##CALCULATES EPSVH (PROB TRANS VECTOR TO HUMAN) DEPENDENT ON TEMPERATURE IN TIME
.MVSEtemp_epsVH<- function(T){
  ep<- 0.001044*T *(T-12.286)*((32.461-T)^(0.5))
  ep[which(ep<0)]<- 0 #fix negative numbers as bio -> trim to zero
  limitingTemp<- 32 #fix asymptote
  ep[which(T>=limitingTemp)]<- 0.001044*limitingTemp *(limitingTemp-12.286)*((32.461-limitingTemp)^(0.5)) #fix NaN for > max T allowed
  return(ep)
}

##CALCULATES MOSQ DEATH RATE COMPONENT DEPENDENT ON TEMPERATURE IN TIME
.MVSEtemp_effect_muV<- function(T){
  ef<- ( 0.8692 -0.159*T +(0.01116*(T^2)) -0.0003408*(T^3) +0.000003809*(T^4) )
  ef[which(ef<0)]<- 0 #fix negative numbers as bio -> trim to zero
  return(ef)
}

##CALCULATES MOSQ DEATH RATE COMPONENT DEPENDENT ON HUMIDITY IN TIME
.MVSEhum_effect_muV<- function(U,meanU){
  ef<- meanU- (U-meanU)/sqrt(1+(U-meanU)^2)
  return(ef)
  ##no need to trim the function for bio meaning
}

##CALCULATES MOSQ INCUBATION TO INFECTION RATE DEPENDENT ON TEMPERATURE IN TIME
.MVSEtemp_effect_gammaV<- function(T){
  Tk<- T+273.15
  R<- 1.987
  ef<- (24.0*( 0.003359* (Tk/298.) * exp((15000./R)*(1/298.-1./Tk)) / (1.+ exp((6.203*(10^21)/R)*(1./(-2.176*(10^30))-1./Tk))) ))
  ef[which(ef<0)]<- 0 #fix negative numbers as bio -> trim to zero
  return(ef)
}

##CALCULATES MOSQ biting RATE COMPONENT DEPENDENT ON HUMIDITY IN TIME
.MVSEhum_effect_aV<- function(U,meanU){
  ef<- (U-meanU)/sqrt(1+(U-meanU)^2)
  return(ef)
  ##no need to trim the function for bio meaning
}

##CALCULATES MOSQ aquatic mortality rate COMPONENT DEPENDENT ON TEMP IN TIME
.MVSEtemp_effect_muA<- function(T){
  ef<- (2.13-0.3797*MVSE_data$T+0.02457*MVSE_data$T^2-0.0006778*MVSE_data$T^3+0.000006794*MVSE_data$T^4)
  return(ef)
 }

##CALCULATES MOSQ aquatic dev rate COMPONENT DEPENDENT ON TEMP IN TIME
.MVSEtemp_effect_epsA<- function(T){
  ef<-  0.131-0.05723*T+0.01164*T^2-0.001341*T^3+0.00008723*T^4-0.000003017*T^5+5.153*10^(-8)*T^6-3.42*10^(-10)*T^7
  return(ef)
 }

 ##CALCULATES MOSQ oviposition rate COMPONENT DEPENDENT ON TEMP IN TIME
.MVSEtemp_effect_theta<- function(T){
  ef<-  -5.4+1.8*T-0.2124*T^2+0.01015*T^3-0.0001515*T^4
  return(ef)
 }

 ##CALCULATES MOSQ egg success rate COMPONENT DEPENDENT ON TEMP IN TIME
.MVSEtemp_effect_C<- function(T){
  ef<-  (-184.8+27.94*T-0.9254*T^2+0.009226*T^3)/100
  return(ef)
 }


 ##CALCULATES MOSQ egg success rate effect COMPONENT DEPENDENT ON RAINFALL IN TIME
.MVSErain_effect_C<- function(R,meanR){
    ef<- (R-meanR)/sqrt(1+(R-meanR)^2)
    return(ef)
 }


# ##SIMULATES A GAUSSIAN WITH CERTAIN MEAN AND SD
# .MVSEpriorGaussian<- function(N, expMean, expSD){
#   return(rnorm(N, mean=expMean, sd=expSD))
# }
#
# ##SIMULATES A GAMMA WITH CERTAIN PARAMETERS
# .MVSEpriorGamma<- function(N, shape, rate){
#   return(rgamma(N, shape=shape, rate=rate))
# }
#
# ##SIMULATES A EXPONENTIAL WITH CERTAIN PARAMETER
# .MVSEpriorExp<- function(N, lambda){
#   return(rexp(N, rate=lambda))
# }

# t-BASED CI
.cierr<- function(dat,conf){
  return(qt(1-conf/2, df=length(dat)-1)*sd(dat)/sqrt(length(dat)))
}

#TAKE COLUMN FROM DATAFRAME OR MATRIX
.takecol<- function(x,colnum){
  return(x[,colnum])
}

#SMOOTHS A TIME SERIES WITH +- N STEPS AVERAGE PER POINT
.smoothUDSeries<- function(series, n){
  return(as.numeric(filter(series, rep(1/n, n), sides=2, circular=TRUE)))
}

#FILLS IN AREA BELOW A CURVE (BOTTOM LIMIT TAKEN AS Y=ZERO)
.plotbellow <- function(X,Y, colIn, colLines){
	line1<- Y
	line2<- rep(0,length(Y))
	xscale <- c(X, rev(X))
	yscale <- c(line1, rev(line2))
	polygon(xscale, yscale, col = colIn, border = colLines, lwd=1.5)
}

#FILLS IN AREA BETWEEN Y1 AND Y2 SERIES
.plotbetweenY <- function(X,Y1,Y2, colIn, colLines){
	line1<- Y1
	line2<- Y2
	xscale <- c(X, rev(X))
	yscale <- c(line1, rev(line2))
	polygon(xscale, yscale, col = colIn, border = colLines, lwd=1.5)
}

##SAMPLES THE POSTERIOR OF BITING RATE
##DOES IT IN TIME, FOR EVERY AVAILABLE TIME STEP
.calculatePostA<- function(Ns, ttt){
  print('Obtaining posterior of biting rate...')
  A<- matrix(0, ncol=ttt, nrow=Ns)
  pb <- txtProgressBar(min=0, max=ttt-1, style = 3)
  for(t in 1:ttt){
      A[,t]<- getPosteriorAVForTimet(mRho=Ns, t=t)
      setTxtProgressBar(pb, t)
  }
  close(pb)
  return(A)
}


##SAMPLES THE POSTERIORS OR / AND PRIORS OF HUMAN VECTOR PARAMETERS
##DOES IT IN TIME, FOR EVERY AVAILABLE TIME STEP
.calculatePostsVLS_EIP_IIP_HIP<- function(Ns, ttt){
  print('Obtaining posteriors of human and vector parameters...')
  VLS<- matrix(0, ncol=ttt, nrow=Ns)
  EIP<- matrix(0, ncol=ttt, nrow=Ns)
  IIP<- matrix(0, ncol=ttt, nrow=Ns)
  HIP<- matrix(0, ncol=ttt, nrow=Ns)
  pb <- txtProgressBar(min=0, max=ttt-1, style = 3)
  for(t in 1:ttt){
      hIncuPer<- replicate(Ns, MVSE_prior_gammaH_rdist(1))
      hInfePer<- replicate(Ns, MVSE_prior_deltaH_rdist(1))
      vIncuPer<- getPosteriorGammaVForTimet(mAlpha=Ns, t=t)
      vLifSPer<- getPosteriorLifespanVForTimet(mEtaRho=Ns, t=t)
      VLS[,t]<- (vLifSPer)
      EIP[,t]<- (vIncuPer)
      IIP[,t]<- (hIncuPer)
      HIP[,t]<- (hInfePer)
      setTxtProgressBar(pb, t)
  }
  close(pb)
  return(list(VLS,EIP,IIP,HIP))
}

#' Find the time of largest index P value per year.
.getPeak<-function(oneyear, CIlevel=0.95){
  ## we can permit other parameters into the function for flexibility
  maxvals=apply(oneyear,2,which.max)
  maxvals[sapply(maxvals, function(x) length(x)==0)] <- NA
  maxvals=unlist(maxvals)
  low=quantile(maxvals,(1-CIlevel)/2,na.rm = TRUE)
  high=quantile(maxvals,1-(1-CIlevel)/2,na.rm = TRUE)
  Mean=mean(maxvals,na.rm = TRUE)
  return(c(Mean,low,high))

}


#' Find the time of highest index P.
.getSeasons<-function(oneyear,timethreshold=0,Pthreshold=1){
  season=which(oneyear>=Pthreshold)
  if(length(season)==0){
    return(NULL)
  }
  diffseason=diff(season,1)
  if(length(season)==0){
    return(season)
  }
  ## find which diff are smaller than timethreshold and fill them in
  fixable=which(diffseason<=timethreshold & diffseason>1)
  for( i in fixable){
    season=union(season,c(season[i]:season[i+1]))
  }
  season=sort(season)
  return(season)
}


#' Compute trend filter with optimal lambda using CV and 1-SE rule
.getFit<-function(values){
  trendfilt1=trendfilter(values,ord=2,approx=T)
  trendfilt1cv =cv.trendfilter(trendfilt1,verbose=F,approx=T)
  return(predict.genlasso(trendfilt1 ,lambda=trendfilt1cv$lambda.1se))

}



###########################################
#### HELP FUNCTIONS
###########################################


#' Shows detailed information on the required csv format for the climatic input
#' data series.
#'
help_inputFormat<- function(){
  cat("
      The climate data is imported from a CSV file (comma separated). It needs to have
      3 columns named H, T and date. H is the humidity series, T the temperature series,
      and date is the date of H and T in the format 'yyyy-mm-dd'. Precipitation data is
      optional and the column name should be P. If no column P exists, calculation of
      Q and V0 is not possible and will be switched off.
      ")
}

#' Shows detailed information global variables used by MVSE.
#'
help_globalVariables<- function(){
  cat("
      This function outputs a long string of characters summarizing each global variable
      used in this package. All global variables have been defined to ease share of content
      between functions but are accessible to the user. The format below is as follows:

      [Name of variable] (name of function where it is set / created) description.

      [MVSE_results] (estimateEcoCoefficients, simulateEmpiricalIndexP, simulateTheoreticalIndexP,
      simulateGenerationTime). List of results. Functions will add elements to this list as they
      are called. Eventually list contains the chains (values) for the posteriors of Rho, Eta
      and Alpha (acessible by name, e.g. MVSE_results$alpha), and also the index P estimations
      (acessible by name, e.g. MVSE_results$indexP, MVSE_results$indexP_theoretical).
      To see available elements use names(MVSE_results).

      [MVSE_data] (setEmpiricalClimateSeries) Stores the empirical climated data
      imported and digested. Columns are 'T', 'H', 'oH', 'year', 'month', 'day', 'date', 'A',
      'Z', 'Y', 'K', 'G'. In order: temperature, humidity (transformed), humidity (original),
      year, month, day, date, effect of humidity on biting rate, effect of temperature on
      lifespan, effect of humidity on lifespan, effect of temperature on probability of
      transmission per bite from V to H, effect of temperature on mosquito incubation
      period. This variable also stores 'B', 'C', 'D', 'E', 'F'. In order:
      effect of temperature on aquatic lifespan, effect of temperature on aquatic development,
      effect of temperature on oviposition, effect of temperature on egg hatch success, effect of
      humidity on egg hatch success.

      [MVSE_data_theoretical] (setTheoreticalClimateSeries) See the entry [MVSE_data] on
      this output. The functionality is the same but in this case it involves the theoretical
      climate variables, not the empirical ones.

      [MVSE_filepath_input] (setEmpiricalClimateSeries) Keeps the filepath to the input file.

      [MVSE_output_tag] (setOutputFilePathAndTag) Keeps the character string that is added to
      all output files.

      [MVSE_prior_epsilonHV_xxx] (setHumanInfPerPrior) 'xxx' can be mean, sd, dist, rdist.
      Each of those variables stores a property of the prior for the human ifectious period.
      mean stores the mean, sd stores the standard deviation, dist stores the name of the
      distribution, rdist the function for random numbers.

      [MVSE_prior_deltaH_xxx] (setHumanMosqTransProbPrior) 'xxx' can be mean, sd, dist, rdist.
      Each of those variables stores a property of the prior for the probability of transmission
      from human to vector per infectious bite. mean stores the mean, sd stores the standard
      deviation, dist stores the name of the distribution, rdist the function for random numbers.

      [MVSE_prior_gammaH_xxx] (setHumanIncPerPrior) 'xxx' can be mean, sd, dist, rdist.
      Each of those variables stores a property of the prior for the human incubation period.
      mean stores the mean, sd stores the standard deviation, dist stores the name of the
      distribution, rdist the function for random numbers.

      [MVSE_prior_muH_xxx] (setHumanLifeExpPrior) 'xxx' can be mean, sd, dist, rdist.
      Each of those variables stores a property of the prior for the human life-span.
      mean stores the mean, sd stores the standard deviation, dist stores the name of the
      distribution, rdist the function for random numbers.

      [MVSE_indexP] (simulateEmpiricalIndexP) Stores all simulated index P time series in a
      matrix: each row a simulation, each column a time step. This is the empirical index P
      simulated using the observed range of climatic variables.

      [MVSE_Q] (simulateEmpiricalIndexP) Stores all simulated Q time series in a
      matrix: each row a simulation, each column a time step. This is the empirical Q
      simulated using the observed range of climatic variables.

      [MVSE_V0] (simulateEmpiricalIndexP) Stores all simulated V0 time series in a
      matrix: each row a simulation, each column a time step. This is the empirical V0
      simulated using the observed range of climatic variables.

      [MVSE_indexP_theoretical] (simulateTheoreticalIndexP) Stores all simulated index P time series in a
      matrix: each row a simulation, each column a time step. This is the theoretical index P
      simulated using a theoretical range of climate variables.

      [MVSE_prior_lev_xxx] (setMosqLifeExpPrior) 'xxx' can be mean, sd, dist, rdist, ddist.
      Each of those variables stores a property of the prior for the mosquito lifespan.
      mean stores the mean, sd stores the standard deviation, dist stores the name of the
      distribution, rdist the function for random numbers, ddist the function for probability
      density function.

      [MVSE_prior_ic_xxx] (setMosqIncPerPrior) 'xxx' can be mean, sd, dist, rdist, ddist.
      Each of those variables stores a property of the prior for the mosquito incubation period.
      mean stores the mean, sd stores the standard deviation, dist stores the name of the
      distribution, rdist the function for random numbers, ddist the function for probability
      density function.

      [MVSE_prior_a_xxx] (setMosqBitingPrior) 'xxx' can be mean, sd, dist, rdist, ddist.
      Each of those variables stores a property of the prior for the mosquito biting rate.
      mean stores the mean, sd stores the standard deviation, dist stores the name of the
      distribution, rdist the function for random numbers, ddist the function for probability
      density function.

      ")

}
