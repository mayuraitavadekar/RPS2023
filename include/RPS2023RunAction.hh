//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RPS2023RunAction_h
#define RPS2023RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023DetectorConstruction;
class RPS2023Run;
class RPS2023PrimaryGeneratorAction;
class RPS2023HistoManager;
class G4GeneralParticleSource;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RPS2023RunAction : public G4UserRunAction
{
  public:
    
    // runAction constructor
    RPS2023RunAction(RPS2023DetectorConstruction*, RPS2023PrimaryGeneratorAction*);
    // destructor
   ~RPS2023RunAction();

  public:
    
    // iRPS2023uilt method invoked at the begining of beamOn()
    virtual G4Run* GenerateRun();  
    
    // iRPS2023uilt methods 
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
                            
  private:
    // create objects of classes
    RPS2023DetectorConstruction*      fDetector;
    RPS2023PrimaryGeneratorAction*    fPrimary;
    RPS2023Run*                       fRun;    
    RPS2023HistoManager*              fHistoManager;
    G4GeneralParticleSource*     gps;
    
  public:
    
    // root output filename
    G4String rootOutputFile;
    G4String descFile;
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


