//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

class RPS2023DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023Run : public G4Run
{
  public:
    
    RPS2023Run(RPS2023DetectorConstruction*);
   ~RPS2023Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);         
    void CountProcesses(const G4VProcess* process, G4int iVol);
    void ParticleCount(G4String, G4double, G4int); 
                          
    virtual void Merge(const G4Run*);
    void EndOfRun();    
   
  private:
  
    // c++ struct for storing data associated with each particle
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };
     
  private:
    // object of class RPS2023DetectorConstruction class
    RPS2023DetectorConstruction* fDetector;
    // instance of particleDefinition class
    G4ParticleDefinition* fParticle;
    // KE or E of particle is measure in eV
    G4double fEkin;

    
    // process counter for each layer
    std::map<G4String,G4int>        fProcCounter1;
    std::map<G4String,G4int>        fProcCounter2;
    std::map<G4String,G4int>        fProcCounter3;
    std::map<G4String,G4int>        fProcCounter4;
    std::map<G4String,G4int>        fProcCounter5;
    std::map<G4String,G4int>        fProcCounter6; // for world
        
    // particle data maps for each layer
    std::map<G4String,ParticleData> fParticleDataMap1;                    
    std::map<G4String,ParticleData> fParticleDataMap2;
    std::map<G4String,ParticleData> fParticleDataMap3;
    std::map<G4String,ParticleData> fParticleDataMap4;
    std::map<G4String,ParticleData> fParticleDataMap5;
    std::map<G4String,ParticleData> fParticleDataMap6; // for world
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

