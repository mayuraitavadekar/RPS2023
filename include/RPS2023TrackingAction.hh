//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class RPS2023DetectorConstruction;
class RPS2023EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023TrackingAction : public G4UserTrackingAction {

  public:  
  
    // constructor & desstructor
    RPS2023TrackingAction(RPS2023DetectorConstruction*, RPS2023EventAction*);
   ~RPS2023TrackingAction() {};
   
    // iRPS2023uilt methods for trackingAction class
    virtual void PreUserTrackingAction(const G4Track*);   
    virtual void PostUserTrackingAction(const G4Track*);
    
  private:
  
    // object of detector construction class
    RPS2023DetectorConstruction* fDetector;
    RPS2023EventAction* fEvent;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

