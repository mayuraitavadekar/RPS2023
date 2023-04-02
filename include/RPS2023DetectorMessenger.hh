#ifndef RPS2023DetectorMessenger_h
#define RPS2023DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RPS2023DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023DetectorMessenger: public G4UImessenger
{
public:
  
  // general construction and destruction function
  RPS2023DetectorMessenger(RPS2023DetectorConstruction*);
  ~RPS2023DetectorMessenger();
  
  // iRPS2023uilt function to override
  void SetNewValue(G4UIcommand*, G4String);
  
  
private:
  
  // object det construction class
  RPS2023DetectorConstruction*       fDetector = nullptr ;
  
  G4UIdirectory*                RPS2023Dir;
  G4UIdirectory*                detDir;
  G4UIdirectory*                pHDir;
  G4UIdirectory*                H2ODir;
  
  G4UIcmdWithAnInteger* H2OCmd;

  // set concentration of H in pH material 
  G4UIcmdWithADouble* pHSetterCmdForH;
  
  // set material of each soil layer
  G4UIcmdWithAString* layerMatCmd;
  
  // set height of each soil layer
  G4UIcmdWithAString* layerHeightCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

