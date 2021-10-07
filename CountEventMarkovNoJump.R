# IDEAL OBSERVER ESTIMATING HIDDEN TRANSITION PROBABILITES THAT GENERATE THE SEQUENCE #
#                                                                                     #
### ### === ORIGINAL MATLAB FILES: CountEventInMemory.m and MarkovNoJump.m  === ### ###
### ### === === === === by Florent Meyniel & Maxime Maheu 2016  === === === === ### ###
#                                                                                     #
### ### = LINK: https://github.com/florentmeyniel/MinimalTransitionProbsModel = ### ###
#                                                                                     #
### ### ===  Partly rewritten in R and modified by Daniel Schlunegger 2019  === ### ###

#######################################################################################

# INPUT:                                                                              #
#     - s: A binary sequence consisting of 0's and 1's (0 == A & 1 == B)              #
#     - omega: A leak / decay factor (default: omega = Inf)                           #
#                                                                                     #
# OUTPUT:                                                                             #
#     - a dataframe containing the following varibles:                                #
#         - trial:        trial number                                                #
#         - AorB:         original binary sequence provided by the caller (s)         #
#         - pXgivenY:     transition probability of X given Y                         #
#                         (maximum a posteriori transition probability of X given Y)  #
#         - predXgivenY:  predictive likelihood of the next stimulus being X given Y  #
#                         ( = mean posterior transition probability)                  #
#         - stdTPgivenX:  standard deviation of the posterior estimate                #
#         - surpXgivenY:  surprise level when the next observation is an X given Y    #
#         - tpPredNextA:  predictive likelihood of the next stimulus being A          #
#                                                                                     #
# NOTE: - If omega = Inf (default) - the observer counts transitions without decay    #
#       - The first event in the sequence is considered arbitrary (drawn randomly),   #
#         This means that analytic posterior values are computed (Beta priors [1,1])  #

#######################################################################################

library(tidyverse)

# Add sub-function to shift down vector n rows and replace with NA
shiftdown <- function(x, n){
  c(rep(NA, n), head(x, -n))
}

# Main function
CountEventMarkov <- function(s, omega = Inf){
  
  # Store input sequence in dataframe
  df <- enframe(s, name = "trial", value = "AorB")
  
  # Drop NA's:
  df <- df %>% 
    drop_na()
  
  # Shift down vector for comparison (DEBUG 09.09.2019: replace s with AorB) 
  df <- df %>% 
    mutate(compareTP = shiftdown(AorB, 1))
  
  # Count the occurence of every transition type
  df <- df %>% 
    mutate(AgivenAcount = ifelse(test = (AorB == 0 & compareTP == 0), yes = 1, no = 0),
           AgivenBcount = ifelse(test = (AorB == 0 & compareTP == 1), yes = 1, no = 0),
           BgivenAcount = ifelse(test = (AorB == 1 & compareTP == 0), yes = 1, no = 0),
           BgivenBcount = ifelse(test = (AorB == 1 & compareTP == 1), yes = 1, no = 0))
  
  # AgivenA
  df <- df %>% 
    group_by(AgivenAcount) %>% 
    mutate(AgivenA = sequence(n())) %>% 
    mutate(AgivenA = ifelse(test = (AgivenAcount == 1), yes = sequence(n()), no = 0))
  
  # AgivenB
  df <- df %>% 
    group_by(AgivenBcount) %>% 
    mutate(AgivenB = sequence(n())) %>% 
    mutate(AgivenB = ifelse(test = (AgivenBcount == 1), yes = sequence(n()), no = 0))
  
  # BgivenA
  df <- df %>% 
    group_by(BgivenAcount) %>% 
    mutate(BgivenA = sequence(n())) %>% 
    mutate(BgivenA = ifelse(test = (BgivenAcount == 1), yes = sequence(n()), no = 0))
  
  # BgivenB
  df <- df %>% 
    group_by(BgivenBcount) %>% 
    mutate(BgivenB = sequence(n())) %>% 
    mutate(BgivenB = ifelse(test = (BgivenBcount == 1), yes = sequence(n()), no = 0))
  
  # Ungroup df
  df <- ungroup(df)
  
  # Replace NA's with 0's
  df[is.na(df)] <- 0
  
  # AgivenA
  df <- df %>% 
    mutate(AgivenA = replace(AgivenA, cumsum(AgivenA != 0) > 0 & AgivenA == 0, NA)) %>% 
    fill(AgivenA)
  
  # BgivenA
  df <- df %>% 
    mutate(BgivenA = replace(BgivenA, cumsum(BgivenA != 0) > 0 & BgivenA == 0, NA)) %>% 
    fill(BgivenA)
  
  # AgivenB
  df <- df %>% 
    mutate(AgivenB = replace(AgivenB, cumsum(AgivenB != 0) > 0 & AgivenB == 0, NA)) %>% 
    fill(AgivenB)
  
  # BgivenB
  df <- df %>% 
    mutate(BgivenB = replace(BgivenB, cumsum(BgivenB != 0) > 0 & BgivenB == 0, NA)) %>% 
    fill(BgivenB)
  
  # Preallocation for the counting variables
  NAgA <- NULL
  NBgA <- NULL
  NAgB <- NULL
  NBgB <- NULL
  
  for (k in 1:length(df$AorB)){
    
    # Calculate decay factor / leaky factor
    decay <- NULL
    for (i in 0:k-1){
      decay[i+1] <- exp(1)**(-i/omega)
    }
    
    # difference to Matlab Code: !!! here 1 and 0! i.E: 1 is same && 2 == 0 !!!
    # Make subs4trn
    subs4trn <- diff(df$AgivenA[1:k] + df$BgivenA[1:k])
    
    # Make trn
    trn <- diff(df$AgivenB[1:k] + df$BgivenA[1:k])
    trnIndex <- df$AgivenBcount[2:k]
    trn <- ifelse(trnIndex == 1, -1, trn)
    
    # Store decay und "turn it around" using rev
    decay <- decay[2:k]
    decay <- rev(decay)
    
    subs4trnIndexGA <- subs4trn == 1
    subs4trnIndexGB <- subs4trn == 0
    
    # Indexing
    trnIndexNAgA <- trn == 0 
    trnIndexNBgA <- trn == 1
    trnIndexNAgB <- trn == -1
    trnIndexNBgB <- trn == 0
    
    leftsideIndexNAgA <- trnIndexNAgA[subs4trnIndexGA]
    leftsideIndexNBgA <- trnIndexNBgA[subs4trnIndexGA]
    leftsideIndexNAgB <- trnIndexNAgB[subs4trnIndexGB]
    leftsideIndexNBgB <- trnIndexNBgB[subs4trnIndexGB]
    
    trnDecayGA <- decay[subs4trnIndexGA]
    trnDecayGB <- decay[subs4trnIndexGB]
    
    # Counting (sum up the weights)
    NAgA[k] <- sum(trnDecayGA[leftsideIndexNAgA])
    NBgA[k] <- sum(trnDecayGA[leftsideIndexNBgA])
    NAgB[k] <- sum(trnDecayGB[leftsideIndexNAgB])
    NBgB[k] <- sum(trnDecayGB[leftsideIndexNBgB])
    
  }
  
  # Compute transition probabilities (Maximum a posteriori transition probabilities)
  df <- df %>% 
    mutate(pAgivenA = NAgA / (NAgA + NBgA),
           pBgivenA = NBgA / (NBgA + NAgA), # == 1 - pAgivenA
           pAgivenB = NAgB / (NAgB + NBgB),
           pBgivenB = NBgB / (NBgB + NAgB)) # == 1 - pAgivenB
  
  # Compute predictive likelihood of what simulus will follow next
  # (Which is also the mean of the distribution)
  df <- df %>%
    mutate(predAgivenA = (NAgA + 1) / (NAgA + NBgA + 2),
           predBgivenA = (NBgA + 1) / (NBgA + NAgA + 2),
           predAgivenB = (NAgB + 1) / (NAgB + NBgB + 2),
           predBgivenB = (NBgB + 1) / (NBgB + NAgB + 2))
  
  # Compute standard deviation of the posterior estimate
  df <- df %>% 
    mutate(stdTPgivenA = sqrt((NBgA + 1)*(NAgA + 1) / ((NBgA + NAgA + 2)**2*(NBgA + NAgA + 3))),
           stdTPgivenB = sqrt((NAgB + 1)*(NBgB + 1) / ((NAgB + NBgB + 2)**2*(NAgB + NBgB + 3))))
  
  # Compute theoretical surprise levels 
  df <- df %>% 
    mutate(surpAgivenA = -log2(predAgivenA),
           surpBgivenA = -log2(predBgivenA),
           surpAgivenB = -log2(predAgivenB),
           surpBgivenB = -log2(predBgivenB))
  
  # Extract the actual predictive likelihood
  df <- df %>% 
    mutate(tpBeforeShiftdown = ifelse(AgivenAcount == 1 | BgivenAcount == 1, predAgivenA, predAgivenB))
  
  df <- df %>% 
    mutate(tpPredNextA = shiftdown(tpBeforeShiftdown, 1))
  
  # Extract the actual standard deviation of the posterior estimate
  df <- df %>% 
    mutate(tpStdBeforeShiftdown = ifelse(AgivenAcount == 1 | BgivenAcount == 1, stdTPgivenA, stdTPgivenB))
  
  df <- df %>% 
    mutate(tpStdNextA = shiftdown(tpStdBeforeShiftdown, 1))

  # Extract the actual surprise level
  df <- df %>% 
    mutate(tpSurpriseBeforeShiftdown = surpAgivenA,
           tpSurpriseBeforeShiftdown = ifelse(BgivenAcount == 1, surpBgivenA, tpSurpriseBeforeShiftdown),
           tpSurpriseBeforeShiftdown = ifelse(AgivenBcount == 1, surpAgivenB, tpSurpriseBeforeShiftdown),
           tpSurpriseBeforeShiftdown = ifelse(BgivenBcount == 1, surpBgivenB, tpSurpriseBeforeShiftdown))
  
  # For the current purpose this is not necessary
  df <- df %>% 
    mutate(tpSurpriseNextA = shiftdown(tpSurpriseBeforeShiftdown, 1))
  
  # Remove some columns before returning the dataframe
  df <- df %>%
    select(
           -compareTP,
           -pAgivenA, -pBgivenA, -pAgivenB, -pBgivenB,
           -AgivenAcount, -AgivenBcount, -BgivenAcount, -BgivenBcount,
           -AgivenA, -AgivenB, -BgivenA, -BgivenB,
           -predAgivenA, -predBgivenA, -predAgivenB, -predBgivenB,
           -stdTPgivenA, -stdTPgivenB,
           -surpAgivenA, -surpBgivenA, -surpAgivenB, -surpBgivenB,
           -tpBeforeShiftdown, -tpStdBeforeShiftdown, 
           #-tpSurpriseBeforeShiftdown,
           #-tpSurpriseNextA,
           -tpStdNextA)
  
  return(df)  
  
}

















