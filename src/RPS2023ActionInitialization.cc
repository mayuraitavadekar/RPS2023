#include "RPS2023ActionInitialization.hh"
#include "RPS2023PrimaryGeneratorAction.hh"
#include "RPS2023RunAction.hh"
#include "RPS2023EventAction.hh"
#include "RPS2023TrackingAction.hh"
#include "RPS2023SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023ActionInitialization::RPS2023ActionInitialization(RPS2023DetectorConstruction* detector)
 : G4VUserActionInitialization(),
   fDetector(detector)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023ActionInitialization::~RPS2023ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023ActionInitialization::BuildForMaster() const
{
  RPS2023RunAction* runAction = new RPS2023RunAction(fDetector, 0);
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023ActionInitialization::Build() const
{
  RPS2023PrimaryGeneratorAction* primary = new RPS2023PrimaryGeneratorAction(); // particle gun : mostly used 
  SetUserAction(primary);
    
  RPS2023RunAction* runAction = new RPS2023RunAction(fDetector, primary);
  SetUserAction(runAction);
  
  RPS2023EventAction* event = new RPS2023EventAction();
  SetUserAction(event);  
  
  RPS2023TrackingAction* trackingAction = new RPS2023TrackingAction(fDetector, event);
  SetUserAction(trackingAction);
  
  RPS2023SteppingAction* steppingAction = new RPS2023SteppingAction(fDetector, event);
  SetUserAction(steppingAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

