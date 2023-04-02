
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023RadonPhysicsList: public G4VModularPhysicsList
{
public:
  RPS2023RadonPhysicsList();
 ~RPS2023RadonPhysicsList();

public:
  virtual void ConstructParticle();
  virtual void SetCuts();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
