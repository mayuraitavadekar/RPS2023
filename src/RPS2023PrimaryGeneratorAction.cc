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
   fParticleGun(0), GPS(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // set initial properties in constructor
  fParticleGun->SetParticleEnergy(0*keV); 
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));

  // also create instance of GPS
  GPS = new G4GeneralParticleSource(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RPS2023PrimaryGeneratorAction::~RPS2023PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete GPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //when you want to use particle gun (almost everytime) ///
  //  use below snippet and comment out GPS snippet ///
  G4double PI, phi, cosTheta, sinTheta;
  G4double px, py, pz;
  G4double rRange, xPos, yPos, zPos;
  G4double minR = 977*nm;  
  G4double maxR = 1000*nm; // recoil range = maxR - minR; maxR = radius of grain

  PI = 3.14159265;
  
  // random momentum direction
  cosTheta = 2*G4UniformRand() - 1.;
  phi = 2*PI*G4UniformRand();
  sinTheta = sqrt(1. - cosTheta*cosTheta);
  px = sinTheta*cos(phi);
  py = sinTheta*sin(phi);
  pz = cosTheta;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));

  // random position between minR and maxR (recoil range)
  G4double r = minR + G4UniformRand() * (maxR - minR);
  phi = 2 * PI * G4UniformRand();
  G4double costheta = 2 * G4UniformRand() - 1;
  G4double sintheta = std::sqrt(1 - costheta * costheta);
  xPos = r * sintheta * std::cos(phi);
  yPos = r * sintheta * std::sin(phi);
  zPos = r * costheta;
  fParticleGun->SetParticlePosition(G4ThreeVector(xPos, yPos, zPos));   

  // initiate the event
  fParticleGun->GeneratePrimaryVertex(anEvent);


  // /// when you want to use GPS (when you need to find recoil range) ///
  // /// use this snippet and comment out particle gun snippet ///
  // GPS->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

