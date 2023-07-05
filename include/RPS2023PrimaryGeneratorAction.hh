// created on August 23, 2022
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef RPS2023PrimaryGeneratorAction_h
#define RPS2023PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4GeneralParticleSource;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    RPS2023PrimaryGeneratorAction();    
   ~RPS2023PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() { return fParticleGun;};
            
  private:
    G4ParticleGun*  fParticleGun;
    G4GeneralParticleSource* GPS; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

