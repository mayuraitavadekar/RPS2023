#include "RPS2023Run.hh"
#include "RPS2023DetectorConstruction.hh"
#include "RPS2023PrimaryGeneratorAction.hh"
#include "RPS2023HistoManager.hh"

#include "G4ProcessTable.hh"
#include "G4Radioactivation.hh"
#include "G4TwoVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023Run::RPS2023Run(RPS2023DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023Run::~RPS2023Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle; // set particle definition
  fEkin = energy; // set initial KE
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023Run::ParticleCount(G4String name, G4double Ekin, G4int iVol)
{

  // layer 1
  if (iVol == 1) {
   std::map<G4String, ParticleData>::iterator it1 = fParticleDataMap1.find(name);
   if ( it1 == fParticleDataMap1.end()) {
     fParticleDataMap1[name] = ParticleData(1, Ekin, Ekin, Ekin);
   }
   else {
     ParticleData& data = it1->second;
     data.fCount++;
     data.fEmean += Ekin;
     //update min max
     G4double emin = data.fEmin;
     if (Ekin < emin) data.fEmin = Ekin;
     G4double emax = data.fEmax;
     if (Ekin > emax) data.fEmax = Ekin; 
   }   
  }
  
  
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023Run::Merge(const G4Run* run)
{
  const RPS2023Run* localRun = static_cast<const RPS2023Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  

  //map: created particles in layer 1   
  std::map<G4String,ParticleData>::const_iterator itc1;
  for (itc1 = localRun->fParticleDataMap1.begin(); 
       itc1 != localRun->fParticleDataMap1.end(); ++itc1) {

    G4String name = itc1->first;
    const ParticleData& localData = itc1->second;   
    if ( fParticleDataMap1.find(name) == fParticleDataMap1.end()) {
      fParticleDataMap1[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap1[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  

  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023Run::EndOfRun() 
{
  // some variables for structured prints on terminal
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  // run condition
  //   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through : ";
  

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
  
  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "********** PARTICLE COUNTS PER LAYER **********" << G4endl;
    
  // particles count in layer 1
  G4cout << "\n List of generated particles in layer 1 GRAIN:" << G4endl;
     
  std::map<G4String,ParticleData>::iterator itc1;               
  for (itc1 = fParticleDataMap1.begin(); itc1 != fParticleDataMap1.end(); itc1++) {
     G4String name = itc1->first;
     ParticleData data = itc1->second;
     G4int count = data.fCount;
     G4double eMean = data.fEmean/count;
     G4double eMin = data.fEmin;
     G4double eMax = data.fEmax;    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
  }
  
  
  
  //remove all contents
  fParticleDataMap1.clear();
  fParticleDataMap2.clear();
  fParticleDataMap3.clear();
  fParticleDataMap4.clear();
  fParticleDataMap5.clear();
  fParticleDataMap6.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
