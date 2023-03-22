// created on August 23, 2022
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef nbPrimaryGeneratorAction_h
#define nbPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"



class G4GeneralParticleSource;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class nbPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    nbPrimaryGeneratorAction();    
   ~nbPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);
    // G4ParticleGun* GetParticleGun() { return fParticleGun;} ;
            
  private:
    // G4ParticleGun*  fParticleGun;
    G4GeneralParticleSource* GPS; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

