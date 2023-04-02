//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"

#include "RPS2023DetectorConstruction.hh"

#include "RPS2023ActionInitialization.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "RPS2023RadonPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  //xiaohang added three lines to change the random engine seed, hopefully it will work.
  G4long seed=time(0); 
  CLHEP::HepRandom::setTheSeed(seed);
  if (argc == 4) CLHEP::HepRandom::setTheSeed(seed+atol(argv[3]));
  CLHEP::HepRandom::showEngineStatus();

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  //construct the run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();  
  if (argc==3) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
  }  

  //set mandatory initialization classes
  RPS2023DetectorConstruction* det= new RPS2023DetectorConstruction;
  runManager->SetUserInitialization(det);
  
  // physics list
  RPS2023RadonPhysicsList* phys = new RPS2023RadonPhysicsList;
  runManager->SetUserInitialization(phys);

 
 
  // include detector geometry
  runManager->SetUserInitialization(new RPS2023ActionInitialization(det));
  

  //initialize G4 kernel
  runManager->Initialize();

  //initialize visualization
  G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive;
   visManager->Initialize();
   UImanager->ApplyCommand("/control/execute vis.mac");
   ui->SessionStart();
   delete ui;
  }
  else  {
   //batch mode
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
  }

  //job termination
  delete visManager;
  delete runManager;
}
