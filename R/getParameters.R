getParameters = function(){
df = list(
    wInf = 100, # Asmptotic weight
    a = 0.35, # Physiological mortality
    h = 20, # Maximum intake
    #kappa = 1, # u
    Rmax = 1000, # close to the trait based ratio for a 25g fish
    years = 1960:2015, # years to simulate
    Fdev = 0.13, # Fishing mortality variation
    Mdev = 0.1, # Natural mortality variation
    nsize = 20, # number of weight bins in survey
    surv.sd = 0.5, # Observation error on survey
    sd.catch = 0.1, # Observation error on catch
    q = 1e-5 # Catchability
    
)
A <- df$h*0.6*0.4 # growth rate 

df$A <- A  #return(df)
return(df)
}
# 
#     wInf = 100 # Asmptotic weight
#     a = 0.35 # Physiological mortality
#     h = 20 # Maximum intake
#     #kappa = 1, # u
#     Rmax = 1000 # close to the trait based ratio for a 25g fish
#     years = 1960:2015 # years to simulate
#     Fdev = 0.13 # Fishing mortality variation
#     Mdev = 0.1 # Natural mortality variation
#     nsize = 20 # number of weight bins in survey
#     surv.sd = 0.5 # Observation error on survey
#     sd.catch = 0.1 # Observation error on catch
