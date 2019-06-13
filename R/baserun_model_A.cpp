
// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
 // Integer data
  DATA_INTEGER(iTimeMax);
  DATA_INTEGER(wlength); // length of w
  DATA_INTEGER(nlength); // Number of weight classes in survey
  DATA_INTEGER(nobs); // number of total survey observations
  DATA_INTEGER(widxM);
  DATA_INTEGER(widxF);

  // Vector data
  DATA_VECTOR(w); // weight vector
  DATA_VECTOR(dw); // dw vector
  DATA_VECTOR(phiMat); // Maturity
  DATA_VECTOR(gg); // growth
  DATA_VECTOR(Sequ); // Initial size distribution
  DATA_VECTOR(Amat); // A matrix for McKendrick Von Voerster
  DATA_VECTOR(Mpredin);// Base mortality
  DATA_VECTOR(mbins); // bins for the survey
  //DATA_VECTOR(surveysel);
  // Index vector
  DATA_IVECTOR(yearidx); // To index years in the binning

  // Scalar data
  DATA_SCALAR(alphaEgg);
  DATA_SCALAR(Zbase);
  DATA_SCALAR(dt);
  //DATA_SCALAR(Fstart);
  //DATA_SCALAR(catchability);
  DATA_SCALAR(SDlog);
  DATA_SCALAR(wInf);
//  DATA_SCALAR(logSDzero);

  // Data for objective function
  DATA_VECTOR(survey); // size spectrum survey
  DATA_VECTOR(Catch); // Observed Catch

  // Parameter integers
  PARAMETER(Rmaxlog);
  PARAMETER(SDlogsurv);
  PARAMETER(logSDM);
  PARAMETER(logSDF);
  PARAMETER(lognF);
  PARAMETER(logmuF);
  PARAMETER(lFstart);
  PARAMETER(lognFs);
  //PARAMETER(logphi);
  PARAMETER(logmuFs);
  PARAMETER(logQ);

 // PARAMETER(logSDzero);

//  PARAMETER(logQ);
  //PARAMETER_VECTOR(Mdevin);
//
  PARAMETER_ARRAY(U);
// Re-transform from log
  Type Rmax = exp(Rmaxlog);
  Type SD = exp(SDlog);
  Type SDsurv = exp(SDlogsurv);
  Type SDF = exp(logSDF);
  Type SDM = exp(logSDM);
  Type nF = exp(lognF);
  Type muF = exp(logmuF);
  Type Fstart = exp(lFstart);
  Type nFs = exp(lognFs);
  //Type phizero = exp(logphi);
  Type muFs = exp(logmuFs);
  Type catchability = exp(logQ); // Catchability

// Initialize vector structures
  array<Type> S(wlength,iTimeMax);
  vector<Type> B(wlength);
  array<Type> Z(wlength,iTimeMax);
  vector<Type> N(wlength);
  vector<Type> Ninit(wlength);
  array<Type> Mpred(wlength,iTimeMax);
  array<Type> Fin(wlength,iTimeMax);
  vector<Type> R(iTimeMax);
  vector<Type> Catchest(iTimeMax);
  vector<Type> Bio(iTimeMax);
  vector<Type> Rp(iTimeMax);
  vector<Type> SSB(iTimeMax);

 // Initialize structures to save
  array<Type> Nsave(wlength,iTimeMax);

  // Save F and M and two sizes
  vector<Type>Msave(iTimeMax);
  vector<Type>Fsave(iTimeMax);


// Random effects (F and M)
  vector<Type> logF(iTimeMax);
  vector<Type> logM(iTimeMax);
//  vector<Type> logR(iTimeMax);

//
//
logF(0) = 0;
logM(0) = 0;

for(int j=0;j<iTimeMax;j++){
      logF(j)=U(0,j);
}

for(int j=0;j<iTimeMax;j++){
      logM(j)=U(1,j);
}

//




// Fill the matrices with zeros - there is one liner for this. Find later
for(int j=1;j<iTimeMax;j++){ // Calculate the total mortality
    for(int i=0;i<wlength;i++){
        S(i,j) = Type(0.0);
        B(i) = Type(0.0);
    }
}


vector<Type> Fsel(wlength); //  Fisheries selectivity

for(int i=0;i<wlength;i++){ // Initialize vectors
  // Fisheries selectivity
 Fsel(i) =pow( (Type(1.0)+pow( (w(i)/(nF * wInf)),-Type(muF))),Type(-1));

}

vector<Type> surveysel(wlength); //  Survey selectivity

for(int i=0;i<wlength;i++){ // Initialize vectors
  // survey selectivity
 surveysel(i) =pow( (Type(1.0)+pow( (w(i)/(nFs * wInf)),-Type(muFs))),Type(-1));

}


for(int i=0;i<wlength;i++){ // Initialize vectors

  N(i) = Sequ(i)*Rmax;//*Type(3.0);

  // Initialize Mortality
  Mpred(i,0) = Mpredin(i)*exp(logM(0));//phizero;
  Fin(i,0) = Fsel(i)*Fstart*exp(logF(0)); // Initial fishing
  Z(i,0) = Zbase + Fin(i,0) + Mpred(i,0);
}




Fsave(0) = Fin(widxF-1,0);
Msave(0) = Mpred(widxM-1,0);

//
//Type cAR = exp(logSDzero);

for(int i=1;i<iTimeMax;i++){ // Calculate the total mortality


    for(int j=0;j < wlength;j++){

      Mpred(j,i) = Mpred(j,i-1)*exp(logM(i)); // Natural

      Fin(j,i) = Fin(j,i-1)*exp(logF(i)); // Fishing

      Z(j,i) = Zbase + Fin(j,i) + Mpred(j,i); // Total
    }
Fsave(i) = Fin(widxF-1,i);
Msave(i) = Mpred(widxM-1,i);
}

//
for(int i=0;i<iTimeMax;i++){ // start time loop

    for (int j=1;j<wlength;j++) { // Calculate the S part of the matrix (A is constant for now)
    S(j,i) = N(j);
    B(j) = 1.0 + gg(j)*dt/dw(j) + Z(j,i)*dt;

    }

   Type Rptemp=0.0; //    Recruitment
    for(int j=0;j<wlength;j++){
    SSB(i) += phiMat(j)*N(j)*w(j)*dw(j);
    }

    Rptemp = SSB(i)*alphaEgg;
    Rp(i) = Rptemp;

    Type Rtemp  = 0.0;
    Rtemp = (Rmax*Rptemp/(Rmax+Rptemp)); //;
    R(i) = Rtemp;

    B(0) = Type(1.0) + gg(0)*dt/dw(0) + Z(0,i)*dt;
    N(0) = (N(0)+ Rtemp*dt/dw(0))/B(0);

    for(int j=1;j<wlength;j++){ // Invert matrices
      N(j) = (S(j,i)-Amat(j)*N(j-1))/B(j);
    }

    Type Btemp = 0.0;
    Type Ctemp = 0.0;
     for(int j=0;j<wlength;j++){
        Btemp += N(j)*w(j)*dw(j);
        Ctemp += Fin(j,i)*N(j)*w(j)*dw(j);
        Nsave(j,i) = N(j);
     }

     Catchest(i) = Ctemp;
     Bio(i) = Btemp;
//



}  // End time loop
//
//
// Calculate stuff
vector<Type> Ntemp(nlength);
vector<Type> Nest(nobs);
Type Q = catchability;
for(int k=0;k<(iTimeMax);k++){ // bin the data

   for(int j=0;j<(nlength);j++){
        Ntemp(j) = 0.0;

        for(int i=0;i<(wlength);i++){

        Type leftEdge = mbins(j);
        Type rightEdge = mbins(j+1);

        if (w(i) >= leftEdge && w(i)<= rightEdge){
            Ntemp(j) += Nsave(i,k)*surveysel(i)*Q; // Fix catchability later

        }

    }
}

int idxtmp;
for(int u=0;u<(nlength);u++){
    idxtmp = yearidx(k)+u;
    Nest(idxtmp) = Ntemp(u);
}

}


Type ans=0.0;
Type kappa = 1e-15; // To prevent log errors
using namespace density;
//
  for(int i=0;i<iTimeMax;i++){
     ans+= -dnorm(logF(i),Type(0.0),SDF,TRUE); // F-Process likelihood
  }
//
// for(int i=1;i<iTimeMax;i++){
//     ans+= -dnorm(logM(i),logM(i-1), SDM, TRUE); // M-Process likelihood
//  }
//
   for(int i=0;i<iTimeMax;i++){
     ans+= -dnorm(logM(i),Type(0.0), SDM, TRUE); // M-Process likelihood
  }



 vector<Type> RBH(iTimeMax); //  stock recruitment

  for(int i=0;i<iTimeMax;i++){
     RBH(i) = (Rmax*Rp(i))/(Rmax+Rp(i));
     ans+= -dnorm(log(R(i)),log(RBH(i)),Type(0.6),TRUE); // R-Process likelihood
  }


for (int i=0;i<nobs;i++){
    ans += -dnorm(log(Nest(i)+kappa), log(survey(i)+kappa), SDsurv, TRUE); // lognormal fit
    }


for (int i=0;i<iTimeMax;i++){
    ans += -dnorm(log(Catchest(i)+kappa), log(Catch(i)+kappa), SD, TRUE); // lognormal fit
    }

//

// Report calculations
ADREPORT(Catchest)
ADREPORT(Bio)
ADREPORT(Fsave)
ADREPORT(Msave)
ADREPORT(R)
ADREPORT(Fsel)
//ADREPORT(Rp)
//ADREPORT(Nest)
//ADREPORT(dw)
//ADREPORT(A)
//ADREPORT(Z)
//ADREPORT(S)
 // Calculate the likelihood

  return ans;
}



