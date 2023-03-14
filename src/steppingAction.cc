//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "nbSteppingAction.hh"

#include "nbDetectorConstruction.hh"
#include "nbRun.hh"
#include "nbEventAction.hh"
#include "nbHistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nbSteppingAction::nbSteppingAction(nbDetectorConstruction* det, nbEventAction* event)
: G4UserSteppingAction(), fDetector(det), fEventAction(event)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nbSteppingAction::~nbSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // instance of G4Run
  nbRun* run = static_cast<nbRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
  
  // instance of analysisManager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // get track
  G4Track* aTrack = aStep->GetTrack();
  
  // get volume of the current step
  G4String pVolume = aStep->GetTrack()->GetVolume()->GetName();
  G4int iVol = 0;
  if (pVolume == fDetector->getNameOfLayer1())   iVol = 1; // grain
  if (pVolume == fDetector->getNameOfLayer2())   iVol = 2; // grain 
  
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();  
  
  G4String pName      = particle->GetParticleName();
  G4int pid           = particle->GetPDGEncoding();
  G4int Z             = particle->GetAtomicNumber();
  G4int A             = particle->GetAtomicMass();
  G4double charge     = particle->GetPDGCharge();  
  G4double x          = aStep->GetPostStepPoint()->GetPosition().x()/um;
  G4double y          = aStep->GetPostStepPoint()->GetPosition().y()/um;
  G4double z          = aStep->GetPostStepPoint()->GetPosition().z()/um;
  G4double pVelocity  = aStep->GetPostStepPoint()->GetVelocity();
  G4double pKE        = aStep->GetPostStepPoint()->GetKineticEnergy()/keV;
  G4double momX         = aStep->GetPostStepPoint()->GetMomentum().x();
  G4double momY         = aStep->GetPostStepPoint()->GetMomentum().y();
  G4double momZ         = aStep->GetPostStepPoint()->GetMomentum().z();
  G4double time       = aStep->GetPostStepPoint()->GetGlobalTime()/s;


  if(pName == "alpha") pid = 1;
  else if(pName == "gamma") pid = 2;
  else if(pName == "e-") pid = 3;
  else if(pName == "e+") pid = 4;
  else if(pName == "Rn222" || A == 222 || Z == 86)
  {
    // stotr pid
    pid = 5;

    // push in upward y
    aTrack->SetMomentumDirection(G4ThreeVector(0., 1., 0.));
  }

  // fill ntuple with id = 2
  G4int id = 2;
  analysisManager->FillNtupleDColumn(id,0, x);
  analysisManager->FillNtupleDColumn(id,1, y);
  analysisManager->FillNtupleDColumn(id,2, z);
  analysisManager->FillNtupleDColumn(id,3, pid);
  analysisManager->FillNtupleDColumn(id,4, Z);
  analysisManager->FillNtupleDColumn(id,5, A);
  analysisManager->FillNtupleSColumn(id,6, pName);
  analysisManager->FillNtupleIColumn(id,7, iVol);
  analysisManager->FillNtupleDColumn(id,8, pVelocity);
  analysisManager->FillNtupleDColumn(id,9, pKE);
  analysisManager->FillNtupleDColumn(id,10, momX);
  analysisManager->FillNtupleDColumn(id,11, momY);
  analysisManager->FillNtupleDColumn(id,12, momZ);
  analysisManager->FillNtupleIColumn(id,13, stepCounter++);
  analysisManager->FillNtupleDColumn(id,14, time);
  analysisManager->FillNtupleDColumn(id,15, charge);
  // add row
  analysisManager->AddNtupleRow(id);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
