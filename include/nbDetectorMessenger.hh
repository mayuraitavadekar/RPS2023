#ifndef nbDetectorMessenger_h
#define nbDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class nbDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class nbDetectorMessenger: public G4UImessenger
{
public:
  
  // general construction and destruction function
  nbDetectorMessenger(nbDetectorConstruction*);
  ~nbDetectorMessenger();
  
  // inbuilt function to override
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  // object det construction class
  nbDetectorConstruction*       fDetector = nullptr ;
  
  G4UIdirectory*                nbDir;
  G4UIdirectory*                detDir;
  G4UIdirectory*                pHDir;
  
  // set concentration of H in pH material 
  G4UIcmdWithADouble* pHSetterCmdForH;
  
  // set material of each soil layer
  G4UIcmdWithAString* layer1MatCmd;
  G4UIcmdWithAString* layer2MatCmd;
  G4UIcmdWithAString* layer3MatCmd;
  G4UIcmdWithAString* layer4MatCmd;
  G4UIcmdWithAString* layer5MatCmd;
  
  // set height of each soil layer
  G4UIcmdWithADoubleAndUnit* layer1HeightCmd;
  G4UIcmdWithADoubleAndUnit* layer2HeightCmd;
  G4UIcmdWithADoubleAndUnit* layer3HeightCmd;
  G4UIcmdWithADoubleAndUnit* layer4HeightCmd;
  G4UIcmdWithADoubleAndUnit* layer5HeightCmd;
  
  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

