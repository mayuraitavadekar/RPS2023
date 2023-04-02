#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023EventAction : public G4UserEventAction
{
  public:
    // constructors and destructors
    RPS2023EventAction();
   ~RPS2023EventAction();

  public:
    
    // inbuilt methods for eventAction class
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    

