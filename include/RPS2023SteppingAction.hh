//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class RPS2023DetectorConstruction;
class RPS2023EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023SteppingAction : public G4UserSteppingAction
{
  public:
  
    // constructor & destructors
    RPS2023SteppingAction(RPS2023DetectorConstruction*, RPS2023EventAction*);
   ~RPS2023SteppingAction();

    // iRPS2023uilt method in steppingAction class
    virtual void UserSteppingAction(const G4Step*);
    
  private:
    
    // objects of detector construction & eventAction class
    RPS2023DetectorConstruction* fDetector;  
    RPS2023EventAction* fEventAction;    

  public:
    G4int stepCounter = 0;
    G4double emanated = -1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

