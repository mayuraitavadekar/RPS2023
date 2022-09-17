// learning resource: https://www.slac.stanford.edu/xorg/geant4/SLACTutorial14/HandsOn1/
// learning resource: https://www.slac.stanford.edu/xorg/geant4/SLACTutorial14/Agenda.html
// learning resource: https://indico.cern.ch/event/647154/contributions/2714212/attachments/1529029/2397032/BookForApplicationDevelopers.pdf

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"

#include "nbDetectorConstruction.hh"
#include "nbPhysicsList.hh"
#include "nbActionInitialization.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//namespace {
//  void PrintUsage() {
//    G4cerr << " Usage: " << G4endl;
//    G4cerr << " neutronBhttps://www.slac.stanford.edu/xorg/geant4/SLACTutorial14/HandsOn1/#ex1aall [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
//    G4cerr << "   note: -t option is available only for multi-threaded mode."
//           << G4endl;
//  }
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

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
  nbDetectorConstruction* det= new nbDetectorConstruction;
  runManager->SetUserInitialization(det);

  nbPhysicsList* phys = new nbPhysicsList;
  runManager->SetUserInitialization(phys);

  runManager->SetUserInitialization(new nbActionInitialization(det));


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
