//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class nbDetectorConstruction;
class nbRun;
class nbPrimaryGeneratorAction;
class nbHistoManager;
class G4GeneralParticleSource;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class nbRunAction : public G4UserRunAction
{
  public:
    
    // runAction constructor
    nbRunAction(nbDetectorConstruction*, nbPrimaryGeneratorAction*);
    // destructor
   ~nbRunAction();

  public:
    
    // inbuilt method invoked at the begining of beamOn()
    virtual G4Run* GenerateRun();  
    
    // inbuilt methods 
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
                            
  private:
    // create objects of classes
    nbDetectorConstruction*      fDetector;
    nbPrimaryGeneratorAction*    fPrimary;
    nbRun*                       fRun;    
    nbHistoManager*              fHistoManager;
    G4GeneralParticleSource*     gps;
    
  public:
    
    // root output filename
    G4String rootOutputFile;
    G4String descFile;
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


