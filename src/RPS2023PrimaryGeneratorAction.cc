// August 23, 2022
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "RPS2023PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>       /* cos */



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023PrimaryGeneratorAction::RPS2023PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   GPS(0)
{
  GPS = new G4GeneralParticleSource(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023PrimaryGeneratorAction::~RPS2023PrimaryGeneratorAction()
{
  delete GPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  GPS->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

