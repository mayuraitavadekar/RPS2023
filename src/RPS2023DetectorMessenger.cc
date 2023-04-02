#include "RPS2023DetectorMessenger.hh"
#include "RPS2023DetectorConstruction.hh"

#include "G4UIdirectory.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023DetectorMessenger::RPS2023DetectorMessenger(RPS2023DetectorConstruction* det)
:fDetector(det)
{ 
  // create ui directory
  RPS2023Dir = new G4UIdirectory("/RPS2023/");
  RPS2023Dir->SetGuidance("UI commands for this detector");
  
  detDir = new G4UIdirectory("/RPS2023/det/");
  detDir->SetGuidance("detector control");

  detDir = new G4UIdirectory("/RPS2023/det/H2O/");
  detDir->SetGuidance("detector control");

  H2OCmd = new G4UIcmdWithAnInteger("/RPS2023/det/H2O/setH2O",this);  
  H2OCmd->SetGuidance("choose amount of H2O in grain");
  H2OCmd->SetParameterName("H2O",false);    
  H2OCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023DetectorMessenger::~RPS2023DetectorMessenger()
{
  delete H2OCmd;
  delete detDir;
  delete RPS2023Dir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == H2OCmd)
  { 
    // update the H2O content variable
    fDetector->updateH2OContent(H2OCmd->GetNewIntValue(newValue)); 

    // update the chemical composition of grain
    fDetector->setGrainMaterial(H2OCmd->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......    

