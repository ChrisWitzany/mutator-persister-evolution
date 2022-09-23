#!/usr/bin/env Rscript

#-------------------------------------------------
# mutator persister evolution (main model)
#-------------------------------------------------

#librarys:
require(adaptivetau)
require(tidyverse)
require(plyr)
require(reshape2)
require(lhs) #as the name suggests, needed for lhs generation
require(KScorrect) #needed to generate uniform log distributions via the dlunif() function

#-------------------------------------------------

run <- F #if true and you run the whole code the model will run and get saved
name <- "test" #if you are running this whole thing you need a name for the created data

#for running via cluster:
#args = commandArgs(trailingOnly=TRUE)
#name <- paste0(args[1])


#-------------------------------------------------
#parameters:
#-------------------------------------------------

params <- c(#pharmacokinetics:
            A_max = 10,         #maximal antimicrobial concentration, (this exact value is only used for plot_single_simulation()/layered_plots())
            k = 0.1,            #AM degradation rate
            tau = 24,           #time period between treatments
            t_max = 24*8,       #maximum run time of simulation
            
            #pharmacodynamics:
            MIC = 1,            #Minimal Inhibitory Concentration
            K = 10^12,          #Carrying capacity
            psimin = -5,        #minimal net growth if antimicorbials present
            psimax = 1,         #max net growth if antimicrobials are absent
            kappa = 1.5,        #Hill Parameter / steepness of the kill curve
            c_R = 0.1,            #cost of resistance
            c_M = 0.02869152,  #cost growth of Hypermutators calculated from data from Montanori et al. 2007
            d = 0.01,             #background death rate of growing cells
            B0 = 10**9,         #size of starting population (N0+P0)
            
            #immunesystem:
            beta = 0.1,  #clearance rate by the immunesyste (only for bacteriostatic ABs)
            
            #mutation rates:
            p_R = 1*10**-6,       #mutation rate of susceptibles to resistant/hypermutators/highpersisters
            to_res = 1,         #resistance possible=1, resistance not possible=0
            h_p = 230,      #fold mutation rate increase of hypermutators
            p_M =  1*10**-5,   #mutation rate of susceptibles to hypermutators
            p_P =  1*10**-5,   #mutation rate of susceptibles to highpersisters
            con_M = 1,         #for turning off M->U mutations
            con_P = 1,         #for turning off P->U mutations

            #switching:
            s_F = 1*1.2*10**-6,  #switching rate of susceptibles into persisterstate
            s_B = 1*0.1,         #switching rate of persisterstate back to susceptible
            h_F = 833.3333,     #fold switching rate increase of highpersisters into persisterstate
            h_B = 10**-3,   # fold switching rate decrease of higherpersister in persisterstate back to normal
            Rp = 0,            #1 = Resistance populations have persistence, 0 = they dont
            
            #A_max range parameters and number of simulations
            min_con = 0,        #lowest concentration of AM
            max_con = 50,       #highest concentration of AM
            increase_by = 0.1,  #stepsize from min to max con
            sim_per_con = 100  #number of simulations per step
            )


#-------------------------------------------------
#starting populations dependent on params:
#-------------------------------------------------

calculate_y <- function(params){
  
  a = params["s_F"]
  b = params["s_B"]
 
  lambda_plus <- as.numeric(params["psimax"] - a - b  + sqrt((-params["psimax"] + a + b)^2 + 4 * params["psimax"] * b))
  v <- as.numeric(c((lambda_plus + b)/a, 1))
  P0 = as.numeric(params["B0"])/(1 + v[1])
  N0 = as.numeric(params["B0"]) - P0   #see Rodriguez-Rojas et al., 2021 for details
  
  y <- c(N     = ceiling(N0), #WT-Susceptible
         N_p   = ceiling(P0), #WT-Susceptible in persisterstate 
         N_R     = 0,         #WT-Resistant
         P    = 0,            #Highpersisters 
         P_p  = 0,            #Highpersisters in persisterstate
         P_R  = 0,            #Highpersisters with acquired Resistance
         M    = 0,            #Hypermutators
         M_p  = 0,            #Hypermutators in persierstate
         M_R  = 0,            #Hypermutators which acquired Resistance
         U   = 0,             #"Mutator-Persisters"
         U_p = 0,             #"Mutator-Persisters" in persisterstate
         U_R = 0,             #"Mutator-Persisters" which acquired resistance
         N_R_p = 0,           #WT-Resistant persister population
         P_R_p = 0,           #resistant Highpersister persister population
         M_R_p = 0,           #resistant Hypermutator persister population
         U_R_p = 0,           #"Mutator-Persisters" persister population
         track_M_U = 0,       #tracking population-whenever M mutates to U it gets the same increase
         track_P_U = 0)       #tracking population - whenever P mutates to U it gets the same increase
         
  return(y)
}

y <- calculate_y(params)

#--------------------------------------------------------------
## transitions and rate function for bactericidal AB treatment 
#--------------------------------------------------------------

transitions = list( #TRANSITIONS AFFECTING N/P:
                    c(N = +1),           #trans1: growth of N
                    c(N = -1, N_p = +1), #trans2: switching from N to P
                    c(N = +1, N_p = -1), #trans3: switching from N_p to N 
                    c(N = -1),           #trans4: effect of antimicrobials on N
                    c(N = -1),           #trans5: capacity limitation
                    
                    #TRANSITIONS AFFECTING N_R:
                    c(N_R = +1),         #trans6: growth of N_R
                    c(N_R = -1),         #trans7: effect of antimicrobials on N_R
                    c(N_R = -1),         #trans8: capacity limitation
                    
                    #TRANSITIONS AFFECTING P/P_p:
                    c(P = +1),           #trans9: growth of P
                    c(P = -1, P_p = +1), #trans10: switching from P to P_p
                    c(P = +1, P_p = -1), #trans11: switching from P_p to P
                    c(P = -1),           #trans12: effect of antimicrobials on P
                    c(P = -1),           #trans13: capacity limitaiton
                    
                    #TRANSITIONS AFFECTING P_R:
                    c(P_R = +1),         #trans14: growth of P_R
                    c(P_R = -1),         #trans15: effect of antimicrobials on P_R
                    c(P_R = -1),         #trans16: capacity limitation
                    
                    
                    #TRANSITIONS AFFECTING M/M_p:
                    c(M = +1),           #trans17: growth of M
                    c(M = -1, M_p = +1), #trans18: switching from M to M_p
                    c(M = +1, M_p = -1), #trans19: switching from M_p to M
                    c(M = -1),           #trans20: effect of antimicrobials on M
                    c(M = -1),           #trans21: capacity limitation
                    
                    #TRANSITIONS AFFECTING M_R:
                    c(M_R = +1),         #trans22: growth of M_R
                    c(M_R = -1),         #trans23: effect of antimicrobials on M_R
                    c(M_R = -1),         #trans24: capacity limitation
                    
                    #TRANSITIONS AFFECTING U/U_p:
                    c(U = +1),           #trans25: growth of U
                    c(U = -1, U_p = +1), #trans26: switching from U to U_p
                    c(U = +1, U_p = -1), #trans27: switching from U_p to U
                    c(U = -1),           #trans28: effect of AM on U
                    c(U = -1),           #trans29: capacity limitation
                    
                    #TRANSITIONS AFFECTING U_R:
                    c(U_R = +1),         #trans30: growth of U_R
                    c(U_R = -1),         #trans31: effect of AM on U_R
                    c(U_R = -1),         #trans32: capacity limitation
                    
                    #MUTATIONS:
                    #"TO"
                    c(N_R = +1, N = -1),               #trans33: mutation N->N_R 
                    c(P = +1, N = -1),                 #trans34: mutation N->P
                    c(M = +1, N = -1),                 #trans35: mutation N->M
                    c(P_R = +1, P = -1),               #trans36: mutation P->P_R
                    c(M_R = +1, M = -1),               #trans37: mutation M->M_R
                    c(U = +1, M = -1, track_M_U = +1), #trans38: mutation M->U
                    c(U = +1, P = -1, track_P_U = +1), #trans39: mutation P->U
                    c(U_R = +1, U = -1),               #trans40: backmutation U_R->U

                    #Resistant persister populations
                    c(N_R = -as.integer(params["Rp"]), N_R_p = +as.integer(params["Rp"])), #trans41: switching of N_R to N_R_p
                    c(N_R = +as.integer(params["Rp"]), N_R_p = -as.integer(params["Rp"])), #trans42: switching of N_R_p to N_R
                    c(P_R = -as.integer(params["Rp"]), P_R_p = +as.integer(params["Rp"])), #trans43: switching of P_R to P_R_p
                    c(P_R = +as.integer(params["Rp"]), P_R_p = -as.integer(params["Rp"])), #trans44: switching of P_p to P_R_p
                    c(M_R = -as.integer(params["Rp"]), M_R_p = +as.integer(params["Rp"])), #trans45: switching of Hm_r to M_R_p
                    c(M_R = +as.integer(params["Rp"]), M_R_p = -as.integer(params["Rp"])), #trans46: switching of Hm_r_p to M_R
                    c(U_R = -as.integer(params["Rp"]), U_R_p = +as.integer(params["Rp"])), #trans47: switching of U_R to U_R_p
                    c(U_R = +as.integer(params["Rp"]), U_R_p = -as.integer(params["Rp"]))  #trans48: switching of U_R_p to U_R
                  
                    )


rate_function <- function(y, params, t){
  
    with(as.list(c(params)),{
      
      A = A_max*exp(-k*(t-((t %/% tau)*tau))) #pharmakokinetics
      
      eN = ((psimax-psimin)*((A/MIC)^kappa))/((A/MIC)^kappa-(psimin/psimax)) #effect of antimicrobial for N
      eR = ((psimax*(1-c_R)-psimin)*((A/(MIC*10))^kappa))/((A/(MIC*10))^kappa-(psimin/(psimax*(1-c_R))))
      #note that eP and eP_R are identical to eN and eR and thus not explicitly needed
      
      eM = ((psimax*(1-c_M)-psimin)*((A/MIC)^kappa))/((A/MIC)^kappa-(psimin/psimax*(1-c_M))) #effect of antimicrobial for M
      eM_R = ((psimax*(1-c_M)*(1-c_R)-psimin)*((A/(10*MIC))^kappa))/((A/(10*MIC))^kappa-(psimin/(psimax*(1-c_M)*(1-c_R)))) #effect of antimicrobial for N_R, 10times as resistant, see MIC
      #note that eU and eU_R are identical to eM and eM_R 
      
      total = sum(y[1:12])
      
      return(c(
        #TRANSITIONS AFFECTING N/N_p:
        psimax*y["N"],                                    #trans1: growth of N
        s_F*y["N"],                                       #trans2: switching from N to P
        s_B*y["N_p"],                                     #trans3: switching from N_p to N
        eN*y["N"],                                        #trans4: effect of antimicrobials on N
        (d+(1-p_R-p_M-p_P)*psimax*(total/K))*y["N"],      #trans5: capacity limitation + background death
        
        #TRANSITIONS AFFECTING N_R:
        psimax*(1-c_R)*y["N_R"],                          #trans6: growth of N_R
        eR*y["N_R"],                                      #trans7: effect of antimicrobials on N_R
        (d+(1-c_R)*psimax*(total/K))*y["N_R"], #trans8: capactiy limitation + background death
        
        #TRANSITIONS AFFECTING P/P_p:
        psimax*y["P"],                                    #trans9: growth of P
        h_F*s_F*y["P"],                                   #trans10: switching from P to P_p
        h_B*s_B*y["P_p"],                                 #trans11: switching from P_p to P
        eN*y["P"],                                        #trans12: effect of antimicrobials on P
        (d+psimax*(1-p_R-p_M)*(total/K))*y["P"],   #trans13: capacity limitation + background death
        
        #TRANSITIONS AFFECTING P_R:
        psimax*(1-c_R)*y["P_R"],                          #trans14: growth of P_R
        eR*y["P_R"],                                      #trans15: effect of antimicrobials on P_R
        (d+(1-c_R)*psimax*(total/K))*y["P_R"], #trans16: capacity limitation + background death
        
        #TRANSITIONS AFFECTING M/M_p:
        psimax*(1-c_M)*y["M"],                                            #trans17: growth of M
        s_F*y["M"],                                                       #trans18: switching from M to M_p
        s_B*y["M_p"],                                                     #trans19: switching from M_p to M
        eM*y["M"],                                                        #trans20: effect of antimicrobials on M
        (d+(1-c_M)*(1-h_p*(p_R+p_P))*psimax*(total/K))*y["M"], #trans21: capacity limit + background death
        
        #TRANSITIONS AFFECTING M_R:
        psimax*(1-c_R)*(1-c_M)*y["M_R"],                                  #trans22: growth of M_R
        eM_R*y["M_R"],                                                    #trans23: effect of antimicrobials on M_R
        (d+(1-c_R)*(1-c_M)*psimax*(total/K))*y["M_R"],     #trans24: capacity limitaiton + background death
        
        #TRANSITIONS AFFECTING U/U_p:
        psimax*(1-c_M)*y["U"],                                            #trans25: growth of U
        h_F*s_F*y["U"],                                                   #trans26: switching from U to U_p
        h_B*s_B*y["U_p"],                                                 #trans27: switching from U_p to U
        eM*y["U"],                                                        #trans28: effect of AM on U
        (d+(1-c_M)*(1-h_p*p_R)*psimax*(total/K))*y["U"],   #trans29: capacity limitaition + background death
        
        #TRANSITIONS AFFECTING U_R:
        psimax*(1-c_R)*(1-c_M)*y["U_R"],                                  #trans30: growth of U_R
        eM_R*y["U_R"],                                                    #trans31: effect of AM on U_R
        (d+(1-c_R)*(1-c_M)*psimax*(total/K))*y["U_R"],     #trans32: capacity limitation + background death
        #MUTATIONS:
        
        #"TO"
        to_res*p_R*psimax*y["N"],               #trans33: mutation N->N_R 
        p_P*psimax*y["N"],                      #trans34: mutation N->P
        p_M*psimax*y["N"],                      #trans35: mutation N->M
        p_R*psimax*to_res*y["P"],               #trans36: mutation P->P_R
        h_p*p_R*psimax*(1-c_M)*to_res*y["M"],   #trans37: mutation M->M_R
        h_p*p_P*psimax*(1-c_M)*y["M"]*con_M,    #trans38: mutation M->U
        p_M*psimax*y["P"]*con_P,                #trans39: mutation P->U
        h_p*p_R*psimax*(1-c_M)*to_res*y["U"],   #trans40: mutation U->U_R

        #Resistant persister populations:
        Rp*s_F*y["N_R"],                        #trans41: switching of N_R to N_R_p
        Rp*s_B*y["N_R_p"],                      #trans42: switching of N_R_p to N_R
        Rp*h_F*s_F*y["P_R"],                    #trans43: switching of P_R to P_R_p
        Rp*h_B*s_B*y["P_R_p"],                  #trans44: switching of P_p to P_R_p
        Rp*s_F*y["M_R"],                        #trans45: switching of Hm_r to M_R_p
        Rp*s_B*y["M_R_p"],                      #trans46: switching of Hm_r_p to M_R
        Rp*h_F*s_F*y["U_R"],                    #trans47: switching of U_R to U_R_p
        Rp*h_B*s_B*y["U_R_p"]                   #trans48: switching of U_R_p to U_R
        
      ))
    })
}

#---------------------------------------------------------------
# transitions and rate function for bacteriostatic AB treatment 
#---------------------------------------------------------------

transitions_bacteriostatic = list( 
  #TRANSITIONS AFFECTING N/P:
  c(N = +1),           #trans1: growth of N (includes AB effect)
  c(N = -1, N_p = +1), #trans2: switching from N to P
  c(N = +1, N_p = -1), #trans3: switching from N_p to N 
  c(N = -1),           #trans4: clearance of N by the immune system
  c(N = -1),           #trans5: capacity limitation
  
  #TRANSITIONS AFFECTING N_R:
  c(N_R = +1),         #trans6: growth of N_R (includes AB effect)
  c(N_R = -1),         #trans7: clearance of N_R by the immune system
  c(N_R = -1),         #trans8: capacity limitation
  
  #TRANSITIONS AFFECTING P/P_p:
  c(P = +1),           #trans9: growth of P (includes AB effect)
  c(P = -1, P_p = +1), #trans10: switching from P to P_p
  c(P = +1, P_p = -1), #trans11: switching from P_p to P
  c(P = -1),           #trans12: clearance of P/P_p by the immune system
  c(P = -1),           #trans13: capacity limitation
  
  #TRANSITIONS AFFECTING P_R:
  c(P_R = +1),         #trans14: growth of P_R (includes AB effect)
  c(P_R = -1),         #trans15: clearance of P_R by the immune system
  c(P_R = -1),         #trans16: capacity limitation
  
  
  #TRANSITIONS AFFECTING M/M_p:
  c(M = +1),           #trans17: growth of M (includes AB effect)
  c(M = -1, M_p = +1), #trans18: switching from M to M_p
  c(M = +1, M_p = -1), #trans19: switching from M_p to M
  c(M = -1),           #trans20: clearance of M by the immune system
  c(M = -1),           #trans21: capacity limitation
  
  #TRANSITIONS AFFECTING M_R:
  c(M_R = +1),         #trans22: growth of M_R (includes AB effect)
  c(M_R = -1),         #trans23: clearance of M_R by the immune system
  c(M_R = -1),         #trans24: capacity limitation
  
  #TRANSITIONS AFFECTING U/U_p:
  c(U = +1),           #trans25: growth of U (includes AB effect)
  c(U = -1, U_p = +1), #trans26: switching from U to U_p
  c(U = +1, U_p = -1), #trans27: switching from U_p to U
  c(U = -1),           #trans28: clearance of U by the immune system
  c(U = -1),           #trans29: capacity limitation
  
  #TRANSITIONS AFFECTING U_R:
  c(U_R = +1),         #trans30: growth of U_R (includes AB effect)
  c(U_R = -1),         #trans31: clearance of U_R by the immune system
  c(U_R = -1),         #trans32: capacity limitation
  
  #MUTATIONS:
  #"TO"
  c(N_R = +1, N = -1),               #trans33: mutation N->N_R 
  c(P = +1, N = -1),                 #trans34: mutation N->P
  c(M = +1, N = -1),                 #trans35: mutation N->M
  c(P_R = +1, P = -1),               #trans36: mutation P->P_R
  c(M_R = +1, M = -1),               #trans37: mutation M->M_R
  c(U = +1, M = -1, track_M_U = +1), #trans38: mutation M->U
  c(U = +1, P = -1, track_P_U = +1), #trans39: mutation P->U
  c(U_R = +1, U = -1),               #trans40: backmutation U_R->U
  
  #Resistant persister populations
  c(N_R = -as.integer(params["Rp"]), N_R_p = +as.integer(params["Rp"])), #trans41: switching of N_R to N_R_p
  c(N_R = +as.integer(params["Rp"]), N_R_p = -as.integer(params["Rp"])), #trans42: switching of N_R_p to N_R
  c(P_R = -as.integer(params["Rp"]), P_R_p = +as.integer(params["Rp"])), #trans43: switching of P_R to P_R_p
  c(P_R = +as.integer(params["Rp"]), P_R_p = -as.integer(params["Rp"])), #trans44: switching of P_p to P_R_p
  c(M_R = -as.integer(params["Rp"]), M_R_p = +as.integer(params["Rp"])), #trans45: switching of Hm_r to M_R_p
  c(M_R = +as.integer(params["Rp"]), M_R_p = -as.integer(params["Rp"])), #trans46: switching of Hm_r_p to M_R
  c(U_R = -as.integer(params["Rp"]), U_R_p = +as.integer(params["Rp"])), #trans47: switching of U_R to U_R_p
  c(U_R = +as.integer(params["Rp"]), U_R_p = -as.integer(params["Rp"]))  #trans48: switching of U_R_p to U_R
  
)


rate_function_bacteriostatic <- function(y, params, t){

  with(as.list(c(params)),{
    
    A = A_max*exp(-k*(t-((t %/% tau)*tau))) #pharmakokinetics
    
    eN = ((psimax-d-psimin)*((A/MIC)^kappa))/((A/MIC)^kappa-(psimin/(psimax-d))) #effect of antimicrobial for N
    eR = ((psimax*(1-c_R)-d-psimin)*((A/(MIC*10))^kappa))/((A/(MIC*10))^kappa-(psimin/(psimax*(1-c_R)-d)))
    #note that eP and eP_R are identical to eN and eR and thus not explicitly needed
    
    eM = ((psimax*(1-c_M)-d-psimin)*((A/MIC)^kappa))/((A/MIC)^kappa-(psimin/psimax*(1-c_M))) #effect of antimicrobial for M
    eM_R = ((psimax*(1-c_M)*(1-c_R)-d-psimin)*((A/(10*MIC))^kappa))/((A/(10*MIC))^kappa-(psimin/(psimax*(1-c_M)*(1-c_R)-d))) #effect of antimicrobial for N_R, 10times as resistant, see MIC
    #note that eU and eU_R are identical to eM and eM_R 
    
    total = sum(y[1:12])
    
    return(c(
      #TRANSITIONS AFFECTING N/N_p:
      (psimax-eN)*y["N"],                                #trans1: growth of N with bacteriostatic AB effect
      s_F*y["N"],                                        #trans2: switching from N to P
      s_B*y["N_p"],                                      #trans3: switching from N_p to N
      beta*y["N"],                                       #trans4: clearance of N by the immune system
      (d+(1-p_R-p_M-p_P)*(psimax-eN)*(total/K))*y["N"],  #trans5: capacity limitation + background death
      
      #TRANSITIONS AFFECTING N_R:
      (psimax*(1-c_R)-eR)*y["N_R"],                          #trans6: growth of N_R with bacteriostatic AB effect
      beta*y["N_R"],                                         #trans7: clearance of N_R by the immune system
      (d+((1-c_R)*psimax-eR)*(total/K))*y["N_R"], #trans8: capactiy limitation + background death
      
      #TRANSITIONS AFFECTING P/P_p:
      (psimax-eN)*y["P"],                                    #trans9: growth of P with bacteriostatic AB effect
      h_F*s_F*y["P"],                                        #trans10: switching from P to P_p
      h_B*s_B*y["P_p"],                                      #trans11: switching from P_p to P
      beta*y["P"],                                           #trans12: clearance of P by the immune system
      (d+(psimax-eN)*(1-p_R-p_M)*(total/K))*y["P"],   #trans13: capacity limitation + background death
      
      #TRANSITIONS AFFECTING P_R:
      (psimax*(1-c_R)-eR)*y["P_R"],                          #trans14: growth of P_R with bacteriostatic AB effect
      beta*y["P_R"],                                           #trans15: clearance of P_R by the immune system
      (d+(psimax*(1-c_R)-eR)*(total/K))*y["P_R"], #trans16: capacity limitation + background death
      
      #TRANSITIONS AFFECTING M/M_p:
      (psimax*(1-c_M)-eM)*y["M"],                                            #trans17: growth of M with bacteriostatic AB effect
      s_F*y["M"],                                                            #trans18: switching from M to M_p
      s_B*y["M_p"],                                                          #trans19: switching from M_p to M
      beta*y["M"],                                                           #trans20: clearance of M by the immune system
      (d+(1-h_p*(p_R+p_P))*((1-c_M)*psimax-eM)*(total/K))*y["M"], #trans21: capacity limit + background death
      
      #TRANSITIONS AFFECTING M_R:
      (psimax*(1-c_R)*(1-c_M)-eM_R)*y["M_R"],                                #trans22: growth of M_R
      beta*y["M_R"],                                                           #trans23: effect of antimicrobials on M_R
      (d+((1-c_R)*(1-c_M)*psimax-eM_R)*(total/K))*y["M_R"],   #trans24: capacity limitaiton + background death
      
      #TRANSITIONS AFFECTING U/U_p:
      (psimax*(1-c_M)-eM)*y["U"],                                            #trans25: growth of U with bacteriostatic AB effect
      h_F*s_F*y["U"],                                                        #trans26: switching from U to U_p
      h_B*s_B*y["U_p"],                                                      #trans27: switching from U_p to U
      beta*y["U"],                                                           #trans28: clearance of U by the immune system
      (d+(1-h_p*p_R)*(psimax*(1-c_M)-eM)*(total/K))*y["U"],   #trans29: capacity limitaition + background death
      
      #TRANSITIONS AFFECTING U_R:
      (psimax*(1-c_R)*(1-c_M)-eM_R)*y["U_R"],                                #trans30: growth of U_R 
      beta*y["U_R"],                                                         #trans31: clearance of U_R by the immune system
      (d+((1-c_R)*(1-c_M)*psimax-eM_R)*(total/K))*y["U_R"],   #trans32: capacity limitation + background death
      
      #MUTATIONS:
      to_res*p_R*(psimax-eN)*y["N"],               #trans33: mutation N->N_R 
      p_P*(psimax-eN)*y["N"],                      #trans34: mutation N->P
      p_M*(psimax-eN)*y["N"],                      #trans35: mutation N->M
      p_R*(psimax-eN)*to_res*y["P"],               #trans36: mutation P->P_R
      h_p*p_R*(psimax*(1-c_M)-eM)*to_res*y["M"],   #trans37: mutation M->M_R
      h_p*p_P*(psimax*(1-c_M)-eM)*y["M"]*con_M,    #trans38: mutation M->U
      p_M*(psimax-eN)*y["P"]*con_P,                #trans39: mutation P->U
      h_p*p_R*(psimax*(1-c_M)-eM)*to_res*y["U"],   #trans40: mutation U->U_R
      
      #Resistant persister populations:
      Rp*s_F*y["N_R"],                        #trans41: switching of N_R to N_R_p
      Rp*s_B*y["N_R_p"],                      #trans42: switching of N_R_p to N_R
      Rp*h_F*s_F*y["P_R"],                    #trans43: switching of P_R to P_R_p
      Rp*h_B*s_B*y["P_R_p"],                  #trans44: switching of P_p to P_R_p
      Rp*s_F*y["M_R"],                        #trans45: switching of Hm_r to M_R_p
      Rp*s_B*y["M_R_p"],                      #trans46: switching of Hm_r_p to M_R
      Rp*h_F*s_F*y["U_R"],                    #trans47: switching of U_R to U_R_p
      Rp*h_B*s_B*y["U_R_p"]                   #trans48: switching of U_R_p to U_R
      
    ))
  })
}

#-------------------------------------------------
# meta colors values for the distinct populations
#-------------------------------------------------

# colors for N
color_N <- "royalblue1"    # "#4876FF"
color_N_p <- "royalblue1"  # "#4876FF"
color_N_R <- "green"#"royalblue4"        # "#00EE00"
color_N_R_p <- "royalblue4"        # "#00EE00"

# color for P
color_P <- "pink"     # "#FFC0CB"
color_P_p <- "pink"   # "#FFC0CB"
color_P_R <- "violet" # "#EE82EE"
color_P_R_p <- "violet" # "#EE82EE"

# colors M
color_M <- "tan1"         # "#FFA54F"
color_M_p <- "tan1"       # "#FFA54F"
color_M_R <- "darkorange" # "#FF8C00"
color_M_R_p<- "darkorange" # "#FF8C00"
 
# colors U
color_U <- "mediumpurple1"   # "#AB82FF"
color_U_p <- "mediumpurple1" # "#AB82FF"
color_U_R <- "purple3"       # "#7D26CD"
color_U_R_p<- "purple3"       # "#7D26CD"


#-------------------------------------------------
#short function to run and plot one model as a oneliner
#-------------------------------------------------

#NOTE: default rate function and transitions are for bactericidal AB treatment
plot_single_simulation <- function(timeframe = c(0, params["t_max"]), show_total = FALSE, trans = transitions, rate_fun = rate_function, legend = FALSE){
  
  sim_results <- data.frame(ssa.adaptivetau(init.values = y, 
                                            transitions = trans, 
                                            rateFunc = rate_fun,
                                            params = params,
                                            tf = as.integer(params["t_max"]),
                                            tl.params = list(epsilon = 0.005)))
  
  sim_results$total <- rowSums(sim_results[ ,2:17], na.rm = TRUE )
  
  par(mar = c(5.1,4.1,4.1,2.1), xpd = TRUE)
    matplot(sim_results[ ,"time"], log10(sim_results[ ,c("N","N_p","N_R","P","P_p","P_R","M","M_p","M_R","U","U_p","U_R","N_R_p","P_R_p","M_R_p","U_R_p",if(show_total == TRUE){"total"} )]), type = "l",
            xlab = "Time (h)", ylab = "Counts (log10)",bty = "L",
            col = c(color_N, color_N_p, color_N_R, color_P, color_P_p,color_P_R,color_M,color_M_p,color_M_R,color_U, color_U_p,color_U_R,color_N_R_p,color_P_R_p,color_M_R_p,color_U_R_p,if(show_total == TRUE){"cyan"}),
            lty = c(1,2,1,1,2,1,1,2,1,1,2,1,2,2,2,2), 
            xlim = timeframe,
            lwd = 1)
    if(legend == TRUE){legend("top", inset = c(0,-0.1), legend = c("N","N_p","N_R","P","P_p","P_R","M","M_p","M_R","U","U_p","U_R","N_R_p","P_R_p","M_R_p","U_R_p",if(show_total == TRUE){"total"}),
           col = c(color_N, color_N_p, color_N_R, color_P, color_P_p,color_P_R,color_M,color_M_p,color_M_R,color_U, color_U_p,color_U_R,color_N_R_p,color_P_R_p,color_M_R_p,color_U_R_p,if(show_total == TRUE){"cyan"}),
           lty = c(1,2,1,1,2,1,1,2,1,1,2,1,2,2,2,2), 
           bty = "n", ncol = 8, cex = 1)}
    par(mar = c(5.1,4.1,4.1,2.1), xpd = FALSE)
  
}


#-------------------------------------------------
# conc_kill_curve function: 
#-------------------------------------------------

#NOTE: default rate function and transitions are for bactericidal AB treatment
conc_kill_curve <- function(min_con = params["min_con"], max_con = params["max_con"], increase_by = params["increase_by"], 
                       sim_per_con = params["sim_per_con"], trans = transitions, rate_fun = rate_function, y = calculate_y(params), 
                       saving, save_every_run = FALSE){
  
  n_increase <- (max_con - min_con) / increase_by #how many increases must happen e.g. number of outer loop runs
  
  for (n in 0:n_increase){
    
    params["A_max"] <- min_con + n * increase_by
    
    for(i in 1:sim_per_con){
      
      df_results <- data.frame(ssa.adaptivetau(init.values = y, 
                                               transitions = trans,
                                               rateFunc = rate_fun,
                                               params = params,
                                               tf = as.integer(params["t_max"]),
                                               tl.params = list(epsilon = 0.005)))
      
      if(save_every_run != FALSE){ #saving every single run in total in addition to only the end state
        
        if(exists("starting_dic") == FALSE){
          starting_dic <- getwd()
          dir.create(save_every_run) #working directory
          setwd(paste0(starting_dic,"/", save_every_run)) #go to that directory
          dir.create("every_run") #create folder for every run
          setwd(paste0(starting_dic,"/", save_every_run,"/every_run")) #change to that folder
          starting_dic <- paste0(starting_dic,"/", save_every_run)
          }
        
        df_results$A_max <- params["A_max"] #add the A_max conc. to the run data
        print("writing run...")
        saveRDS(object = df_results, file = paste0("run_inc_",n,"_no_",i,".csv"))
        print("done.")
        
      }
      
      #readout
      last_line <- tail(df_results, 1) 
      
      last_line$total <- sum(last_line["N"],last_line["N_p"],last_line["N_R"],
                             last_line["P"],last_line["P_p"],last_line["M"], last_line["M_p"],
                             last_line["M_R"],last_line["P_R"],last_line["U"],
                             last_line["U_p"],last_line["U_R"],
                             last_line["N_R_p"],last_line["P_R_p"],last_line["M_R_p"],last_line["U_R_p"])
      
      last_line$A_max <- params["A_max"]
      
      #read out time of clearance
      last_line$time_of_clearance <- ifelse(length(head(df_results[rowSums(df_results[,2:17]) == 0,], 1)$time)==0, 
                                  NA,
                                  head(df_results[rowSums(df_results[,2:17]) == 0,], 1)$time) 
      
      #save to df 
      if(i == 1 & n == 0){
        
        save_df <- last_line #create save_df from the readout if first run 

        }else{
        
        save_df[nrow(save_df) + 1 , ] = last_line #add data to new row of save_df
        
        }
    }
    print(paste(n/n_increase*100, "% of total simulations done"))
  }
  
  if(saving != FALSE){
    
    if(exists("starting_dic") == TRUE){setwd(starting_dic)}else{starting_dic <- getwd()} #to avoid conflict with every single run
    
    dir.create(saving) #create folder
    setwd(paste0(starting_dic,"/",saving)) 
    
    saveRDS(object = save_df, file = paste0(saving,".Rds"))
    write.csv(y, paste0(saving, "_y.csv"), row.names = TRUE)
    write.csv(params, paste0(saving, "_params.csv"), row.names = TRUE)
    
  }
  
  return(save_df)
  
}


#-------------------------------------------------
#sensitivity analysis
#-------------------------------------------------

sensitivity_analysis <- function(runs, name, trans = transitions, rate_fun = rate_function){
  
  params_space <- matrix( c( A_max_min = 10, A_max_max = 35, A_max_dist = "unif",
                             s_F_min = 10**-9, s_F_max = 1, s_F_dist = "log", 
                             s_B_min = 10**-9, s_B_max = 1, s_B_dist = "log",
                             p_R_min = 10**-12, p_R_max = 10**-5, p_R_dist = "log"),
                          ncol = 3, byrow = TRUE)
  
  #latin hypersquare sampling
  lhs <- randomLHS(n = runs,                # n = row = number of samples 
                   k = nrow(params_space))  # k = col =  number of parameters
  
  lhs_params <- matrix(nrow = nrow(lhs), ncol = ncol(lhs))   #each row will consist of one sample of the parameter space.
  colnames(lhs_params) <- c("A_max", "s_F","s_B","p_R") #easier indexing
  
  #transform lhs via params_space to actual parameter values
  for (i in 1:nrow(params_space)){
    if(params_space[i,3]=="unif"){
      lhs_params[,i] <- qunif(lhs[,i], 
                              min = as.numeric(params_space[i,1]),
                              max = as.numeric(params_space[i,2])) 
    } else if(params_space[i,3]=="log"){
      lhs_params[,i] <- qlunif(lhs[,i],
                               min = as.numeric(params_space[i,1]), 
                               max = as.numeric(params_space[i,2]))
      
    }}
  
  for(i in 1:runs){
    
    #change params according to lhs_params
    params["s_F"] = lhs_params[i,"s_F"]
    params["A_max"] = lhs_params[i,"A_max"]
    params["s_B"] = lhs_params[i,"s_B"]
    params["p_R"] = lhs_params[i,"p_R"]
    
    #adapt dependent parameters accordingly
    params["p_P"] = 10*params["p_R"]
    params["p_M"] = 10*params["p_R"]
    
    #adapt y to the random parameters, i.e. avoid N/N_p ratio to be not realistic for the respective switching rates
    y <- calculate_y(params = params)
    
    #run the simulation with this set of parameters (and y)
    df_results <- data.frame(ssa.adaptivetau(init.values = y, 
                                             transitions = trans,
                                             rateFunc = rate_fun,
                                             params = params,
                                             tf = as.integer(params["t_max"]),
                                             tl.params = list(epsilon = 0.005)))
    
    last_line <- tail(df_results,1) #main readout
    
    #add total number of cells
    last_line$total <- sum(last_line["N"],last_line["N_p"],last_line["N_R"],
                           last_line["P"],last_line["P_p"],last_line["M"], last_line["M_p"],
                           last_line["M_R"],last_line["P_R"],last_line["U"],
                           last_line["U_p"],last_line["U_R"],
                           last_line["N_R_p"],last_line["P_R_p"],last_line["M_R_p"],last_line["U_R_p"])
    
    #read out Time of clearance
    last_line$time_of_clearance <- ifelse(length(head(df_results[rowSums(df_results[,2:17]) == 0,], 1)$time)==0, 
                                          NA,
                                          head(df_results[rowSums(df_results[,2:17]) == 0,], 1)$time)
    
    #add parameters to last_line and create save_df
    if(i == 1){
      save_df <- cbind(last_line,data.frame(as.list(params))) #create save_df from the readout if first run
    }else{
      save_df[nrow(save_df) + 1 , ] = cbind(last_line,data.frame(as.list(params))) #add data to new row of save_df
    }
    
    #progress report:
    print(paste0(i*100/runs,"% of runs done"))
    
  }
  
  if(saving != FALSE){
    
    starting_dic <- getwd()
    dir.create(saving)
    setwd(paste0(starting_dic,"/",saving))
    
    saveRDS(object = save_df, file = paste0(saving,".Rds"))
    write.csv(y, paste0(saving, "_y.csv"), row.names = TRUE)
    write.csv(params_space, paste0(saving, "_params_space.csv"), row.names = TRUE)
    
  }
  
  return(save_df)
  
}



#-------------------------------------------------
#run simulation by hand / opposed to function:
#-------------------------------------------------

# #bactericidal:
# sim_results = ssa.adaptivetau(init.values = y,
#                               transitions = transitions,
#                               rateFunc = rate_function,
#                               params = params,
#                               tf = as.integer(params["t_max"]),
#                               tl.params = list(epsilon = 0.005))

# #bacteriostatic:
# params["psimin"] = -0.01
# sim_results = ssa.adaptivetau(init.values = y,
#                               transitions = transitions_bacteriostatic,
#                               rateFunc = rate_function_bacteriostatic,
#                               params = params,
#                               tf = as.integer(params["t_max"]),
#                               tl.params = list(epsilon = 0.005))

#---------------------------------------------------
# functionality check and running main simulations
#---------------------------------------------------

#NOTE: default rate function and transitions are for bactericidal AB treatment

#bactericidal:
#single simulation for params
plot_single_simulation()
#A_max variation
if(run == TRUE){run <- conc_kill_curve(saving = name)}

# #bacteriostatic:
# params["psimin"] = -0.01 #set maximal killing close to zero for bacteriostatic simulations!
# plot_single_simulation(trans = transitions_bacteriostatic, rate_fun = rate_function_bacteriostatic)
# #A_max variation
# if(run == TRUE){run <- conc_kill_curve(saving = name, trans = transitions_bacteriostatic, rate_fun = rate_function_bacteriostatic)}
