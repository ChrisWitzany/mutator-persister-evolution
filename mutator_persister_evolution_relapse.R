#!/usr/bin/env Rscript

#-------------------------------------------------
# mutator persister evolution (relapse simulations)
#-------------------------------------------------

#librarys:
require(adaptivetau)
require(tidyverse)
require(plyr)
require(reshape2)

#-------------------------------------------------

# relapse simulations takes in simulation data during treatment as starting points

name <- "test" #name of the .Rds file (the data of the treatment simulation)

#read in
data <- readRDS(paste0(name,".Rds"))
data <- rapply(data, as.factor, classes = "character", how = "replace")
data.table::setDF(data)

#for cluster:
#args = commandArgs(trailingOnly=TRUE)
#name <- paste0(args[1])

run = FALSE # if true and you run the whole code relapse runs and gets saved
name <- "test" #necessary for saving


y <- deframe(read.csv(paste0(name,"_y.csv"))) #read in y (starting numbers)
params <- deframe(read.csv(paste0(name,"_params.csv"))) #and parameters used in the treatment simulation


#-------------------------------------------------
#transitions and rate function are needed for the adaptivetau simulation:
#-------------------------------------------------

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
  #"BACK"
  c(N_R = -as.integer(params["BM"]), N = +as.integer(params["BM"])), #trans41: backmutation N_R->N 
  c(P = -as.integer(params["BM"]), N = +as.integer(params["BM"])),   #trans42: backmutation P->N
  c(M = -as.integer(params["BM"]), N = +as.integer(params["BM"])),   #trans43: backmutation M->N
  c(P_R = -as.integer(params["BM"]), P = +as.integer(params["BM"])), #trans44: backmutation P_R->P
  c(M_R = -as.integer(params["BM"]), M = +as.integer(params["BM"])), #trans45: backmutation M_R->M
  c(U = -as.integer(params["BM"]), M = +as.integer(params["BM"])),   #trans46: backmutation U->M
  c(U = -as.integer(params["BM"]), P = +as.integer(params["BM"])),   #trans47: backmutation U->P
  c(U_R = -as.integer(params["BM"]),U = +as.integer(params["BM"])),  #trans48: backmutation U_R->U
  
  #Resistant persister populations
  c(N_R = -as.integer(params["Rp"]), N_R_p = +as.integer(params["Rp"])), #trans49: switching of N_R to N_R_p
  c(N_R = +as.integer(params["Rp"]), N_R_p = -as.integer(params["Rp"])), #trans50: switching of N_R_p to N_R
  c(P_R = -as.integer(params["Rp"]), P_R_p = +as.integer(params["Rp"])), #trans51: switching of P_R to P_R_p
  c(P_R = +as.integer(params["Rp"]), P_R_p = -as.integer(params["Rp"])), #trans52: switching of P_p to P_R_p
  c(M_R = -as.integer(params["Rp"]), M_R_p = +as.integer(params["Rp"])), #trans53: switching of Hm_r to M_R_p
  c(M_R = +as.integer(params["Rp"]), M_R_p = -as.integer(params["Rp"])), #trans54: switching of Hm_r_p to M_R
  c(U_R = -as.integer(params["Rp"]), U_R_p = +as.integer(params["Rp"])), #trans55: switching of U_R to U_R_p
  c(U_R = +as.integer(params["Rp"]), U_R_p = -as.integer(params["Rp"]))  #trans56: switching of U_R_p to U_R
  
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
          (d+(1-c_R)*(1-p_R*BM)*psimax*(total/K))*y["N_R"], #trans8: capactiy limitation + background death
          
          #TRANSITIONS AFFECTING P/P_p:
          psimax*y["P"],                                    #trans9: growth of P
          h_F*s_F*y["P"],                                   #trans10: switching from P to P_p
          h_B*s_B*y["P_p"],                                 #trans11: switching from P_p to P
          eN*y["P"],                                        #trans12: effect of antimicrobials on P
          (d+psimax*(1-p_R-p_M-p_P*BM)*(total/K))*y["P"],   #trans13: capacity limitation + background death
          
          #TRANSITIONS AFFECTING P_R:
          psimax*(1-c_R)*y["P_R"],                          #trans14: growth of P_R
          eR*y["P_R"],                                      #trans15: effect of antimicrobials on P_R
          (d+(1-c_R)*(1-p_R*BM)*psimax*(total/K))*y["P_R"], #trans16: capacity limitation + background death
          
          #TRANSITIONS AFFECTING M/M_p:
          psimax*(1-c_M)*y["M"],                                            #trans17: growth of M
          s_F*y["M"],                                                       #trans18: switching from M to M_p
          s_B*y["M_p"],                                                     #trans19: switching from M_p to M
          eM*y["M"],                                                        #trans20: effect of antimicrobials on M
          (d+(1-c_M)*(1-h_p*(p_R+p_P)-h_p*p_M*BM)*psimax*(total/K))*y["M"], #trans21: capacity limit + background death
          
          #TRANSITIONS AFFECTING M_R:
          psimax*(1-c_R)*(1-c_M)*y["M_R"],                                  #trans22: growth of M_R
          eM_R*y["M_R"],                                                    #trans23: effect of antimicrobials on M_R
          (d+(1-c_R)*(1-c_M)*(1-h_p*p_R*BM)*psimax*(total/K))*y["M_R"],     #trans24: capacity limitaiton + background death
          
          #TRANSITIONS AFFECTING U/U_p:
          psimax*(1-c_M)*y["U"],                                            #trans25: growth of U
          h_F*s_F*y["U"],                                                   #trans26: switching from U to U_p
          h_B*s_B*y["U_p"],                                                 #trans27: switching from U_p to U
          eM*y["U"],                                                        #trans28: effect of AM on U
          (d+(1-c_M)*(1-h_p*(p_R-(p_M+p_P)*BM))*psimax*(total/K))*y["U"],   #trans29: capacity limitaition + background death
          
          #TRANSITIONS AFFECTING U_R:
          psimax*(1-c_R)*(1-c_M)*y["U_R"],                                  #trans30: growth of U_R
          eM_R*y["U_R"],                                                    #trans31: effect of AM on U_R
          (d+(1-c_R)*(1-c_M)*(1-h_p*p_R*BM)*psimax*(total/K))*y["U_R"],     #trans32: capacity limitation + background death
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
          
          #"BACK"
          BM*p_R*(1-c_R)*psimax*y["N_R"],         #trans41: backmutation N_R->N 
          BM*p_P*psimax*y["P"],                   #trans42: backmutation P->N
          BM*p_M*h_p*(1-c_M)*psimax*y["M"],       #trans43: backmutation M->N
          BM*p_R*(1-c_R)*psimax*y["P_R"],         #trans44: backmutation P_R->P
          BM*p_R*(1-c_R)*(1-c_M)*psimax*y["M_R"], #trans45: backmutation M_R->M
          BM*p_M*h_p*(1-c_M)*psimax*y["U"],       #trans46: backmutation U->M
          BM*p_P*h_p*(1-c_M)*psimax*y["U"],       #trans47: backmutation U->P
          BM*p_R*h_p*psimax*(1-c_M-c_R)*y["U_R"], #trans48: backmutation U_R->U
          
          #Resistant persister populations:
          Rp*s_F*y["N_R"],                        #trans49: switching of N_R to N_R_p
          Rp*s_B*y["N_R_p"],                      #trans50: switching of N_R_p to N_R
          Rp*h_F*s_F*y["P_R"],                    #trans51: switching of P_R to P_R_p
          Rp*h_B*s_B*y["P_R_p"],                  #trans52: switching of P_p to P_R_p
          Rp*s_F*y["M_R"],                        #trans53: switching of Hm_r to M_R_p
          Rp*s_B*y["M_R_p"],                      #trans54: switching of Hm_r_p to M_R
          Rp*h_F*s_F*y["U_R"],                    #trans55: switching of U_R to U_R_p
          Rp*h_B*s_B*y["U_R_p"]                   #trans56: switching of U_R_p to U_R
          
        ))
      })
}
 
#-------------------------------------------------
# relapse function
#-------------------------------------------------

relapse <- function(data, relapse_duration = 87600, detection_limit = 10**5, relapse_value = 10**6){

  relapse_params <- params #creating stand-in params variable #actually this is probably not optimal.
  
  #relapse is assumed to be in the absence of additional ABs, however ABs are still present at the end of the previous treatment
  #we assume treatment ends at the point in time where you would give AB but dont, therefore degradation happened for a duration equal to tau
  #hence the initial leftover AB conc. for relapse will be calculated with:
  A <- function(t, current_A_max){with(as.list(c(relapse_params)), { current_A_max*exp(-k*(t-((t%/% tau)*tau))) })}
  
  #vectorizing as much as possible:
  data$relapse_type <- ifelse(data$total>=detection_limit, "acute_failure", NA)
  data$relapse_type <- ifelse(data$total==0, "previous_treatment_success", data$relapse_type)
  
  #now get only the ones which are small non-zero populations:
  relapse_pools <- subset(data, data$total < detection_limit & data$total != 0)
  
  #to avoid multiple indexing
  total_iterations <- length(relapse_pools[,1])
  percent <- total_iterations/100
  
  #-----------------------------
  #create the readout variable here:
  save_matrix <- data.frame(matrix(ncol = 22, nrow = total_iterations, dimnames = list(NULL, c("time","N","N_p","N_R","P","P_p","P_R","M","M_p","M_R","U","U_p","U_R","N_R_p","P_R_p","M_R_p","U_R_p",
                                                                                               "track_M_U","track_P_U","total","A_max","relapse_type"))))
  print("loop start") #just a quick check to see that it works.
  
  for(i in 1:length(relapse_pools$time)){ #check previous end states and categorize accordingly
    
    #to avoid multiple indexing:
    current <- relapse_pools[i,]
    current_total <- current$total
    
    #calculate leftover A_max concentration for this (current) state
    relapse_params["tau"] = params["tau"] #reset to initial tau for correct calculation!
    relapse_params["A_max"] = A(t = relapse_params["tau"], current_A_max = current$A_max)  #set this as the max_AB for during the relapse, in combination with tau=very large this equals one time appliance of leftover_A
    relapse_params["tau"] = 10**14 #setting tau very,very high to ensure no AB is given at any point
    
    #run simulation with these as input:
    relapse_sim <- data.frame(ssa.adaptivetau(init.values = c(unlist(current[, 2:19])), #2:19 because I do not need time, A_max, time_of_clearance, relapse_type - using unlist to keep the names!
                                                    transitions = transitions, 
                                                    rateFunc = rate_function,
                                                    params = relapse_params,
                                                    tf = as.integer(relapse_duration),
                                                    tl.params = list(epsilon = 0.005)))
    
    #final total value of the sim:
    endstate_relapse_sim <- as.numeric(rowSums(tail(relapse_sim[,2:17], 1)))
      
    #no relapse happened, i.e. no growth
    if(endstate_relapse_sim == current_total){
        
        relapse_line <- tail(relapse_sim, 1)
        #calculating total:
        relapse_line$total <- sum(relapse_line[2:17])
        #need to retain previous max AM conc.
        relapse_line$A_max <- current$A_max
        #add a label to keep track what happened here:
        relapse_line$relapse_type <- "no_relapse"
        #save and go to next i
        save_matrix[i, ] = as.list(c(relapse_line))
        if((i%%round(percent))==0){print(paste0((i/total_iterations)*100, "% of relapse sims done"))}
        next()

    }
      
    #they got eradicated by the leftover concentration
    if(endstate_relapse_sim < current_total){
        
        relapse_line <- tail(relapse_sim, 1)
        #calculating total:
        relapse_line$total <- sum(relapse_line[2:17])
        #need to retain previous max AM conc.
        relapse_line$A_max <- current$A_max
        #add label:
        relapse_line$relapse_type <- ifelse(endstate_relapse_sim==0,"leftover_clearance", "leftover_decline")
        #save and go to next i
        save_matrix[i, ] = as.list(c(relapse_line))
        if((i%%round(percent))==0){print(paste0((i/total_iterations)*100, "% of relapse sims done"))}
        next()

      }
      
    #relapse, i.e. growth
    if(endstate_relapse_sim > current_total){
        
        #calculating total:
        relapse_sim$total <- rowSums(relapse_sim[,2:17])
        
        #get the time of real relapse e.g. when total pop reaches the relapse value
        #because a possibility is that there is growth but the relapse value is not reached, we have to make this a bit more complicated...
        ifelse(max(relapse_sim$total) >= relapse_value,
               time_of_relapse <- head(relapse_sim[relapse_sim["total"] >= relapse_value,], 1)$time,
               time_of_relapse <- relapse_duration)  
       
        #now cut off the data at that point in time
        relapse_line <- relapse_sim[relapse_sim["time"] == time_of_relapse, ]
        
        #need to retain previous max AM conc.
        relapse_line$A_max <- current$A_max
        
        # add relapse label
        relapse_line$relapse_type <- "relapse"
        
        #save and go to next i
        save_matrix[i, ] = as.list(c(relapse_line))
        if((i%%round(percent))==0){print(paste0(round((i/total_iterations)*100), "% of relapse sims done"))}
        next()

    }
  
    
  }
  
  #now that the loop is done, remerge the save_matrix with the acute treatment failures and treatment successes (which we ignored till know)
  success_and_acute_fail <- subset(data[,-c(22)], data$relapse_type=="acute_failure"|data$relapse_type=="previous_treatment_success") #only the columns which relapse_line also has
  merged_save_matrix <- rbind(success_and_acute_fail, save_matrix)
  
  return(merged_save_matrix)
}

#---------------------------------
# run it:
#---------------------------------


if(run == TRUE){
 
  data_relapse <- relapse(data) #corresponds to 10 years

  #save the data
  saveRDS(data_relapse, paste0(name,"_relapse.Rds"))

}



