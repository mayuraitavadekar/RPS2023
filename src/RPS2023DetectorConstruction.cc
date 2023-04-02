#include "RPS2023DetectorConstruction.hh"

#include "RPS2023DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <string>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// constructor to initialization values of shellPV, overlap checks 
// set data member values from config.txt file
RPS2023DetectorConstruction::RPS2023DetectorConstruction()
 : G4VUserDetectorConstruction(),
   shellPV_1(nullptr),
   fCheckOverlaps(true),
   fractionMassForH(0.),
   fractionMassForOH(0.),
   pHValue(4.0),
   fdetectorMessenger(0)
{
  // set the values for H and OH concentration by default
  pHValue = 4;
  fractionMassForH = (pHValue*100.00)/14.00;
  fractionMassForOH = 100.00-fractionMassForH;
  
  // call detector messenger
  fdetectorMessenger = new RPS2023DetectorMessenger(this);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// destructor
RPS2023DetectorConstruction::~RPS2023DetectorConstruction()
{ 
    delete fdetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// member function of header file
G4VPhysicalVolume* RPS2023DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();

  G4cout << "soil grain size : " << grainSize << G4endl; 

  // define chemical composition maps
  DefineChemicalComps();
  fillGrainWithChemComps();

  DefinePoreChemicalComps();
  fillPoresWithChemComps();
  
  // initialize grain material
  grainMaterial = grainComp0;
  poreMaterial = poreComp0;

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023DetectorConstruction::DefineMaterials()
{   
  // Lead material defined using NIST Manager
  nistManager = G4NistManager::Instance();
  nistManager->SetVerbose(1);
  
  //
  // create all your materials in this function
  //
  
  // define pure nist materials
  H = nistManager->FindOrBuildMaterial("G4_H");
  Fe = nistManager->FindOrBuildMaterial("G4_Fe");
  Mn = nistManager->FindOrBuildMaterial("G4_Mn");
  Ra = nistManager->FindOrBuildMaterial("G4_Ra");
  Ca = nistManager->FindOrBuildMaterial("G4_Ca");
  C = nistManager->FindOrBuildMaterial("G4_C");
  N = nistManager->FindOrBuildMaterial("G4_N");
  O = nistManager->FindOrBuildMaterial("G4_O");
  Na = nistManager->FindOrBuildMaterial("G4_Na");
  Mg = nistManager->FindOrBuildMaterial("G4_Mg");
  P = nistManager->FindOrBuildMaterial("G4_P");
  Si = nistManager->FindOrBuildMaterial("G4_Si");
  K = nistManager->FindOrBuildMaterial("G4_K");
  Al = nistManager->FindOrBuildMaterial("G4_Al");
  I = nistManager->FindOrBuildMaterial("G4_I");
  S = nistManager->FindOrBuildMaterial("G4_S");
  Ti = nistManager->FindOrBuildMaterial("G4_Ti");
  H2OVapor = nistManager->FindOrBuildMaterial("G4_WATER_VAPOR");

  // define elements
  elO = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  elH = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  elC  = nistManager->FindOrBuildElement("C");
  elN  = nistManager->FindOrBuildElement("N");
  elMn = nistManager->FindOrBuildElement("Mn");
  elK = nistManager->FindOrBuildElement("K");   // Potassium
  elMg = nistManager->FindOrBuildElement("Mg"); // magnesium
  elS = nistManager->FindOrBuildElement("S");   // Sulfur
  elP = nistManager->FindOrBuildElement("P");   // Phosphorus
  elSi = new G4Element("Silicon", "Si", z=14., a= 28.0855*g/mole);
  
  // define chemical compounds
  Quartz = new G4Material("Quartz", 2.64 *g/cm3, ncomponents= 2);
  Quartz-> AddElement(elSi, natoms=1);
  Quartz-> AddElement(elO,  natoms=2);

  H2O = nistManager->FindOrBuildMaterial("G4_WATER");
  SiO2 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  Al2O3 = nistManager->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  Fe2O3 = nistManager->FindOrBuildMaterial("G4_FERRIC_OXIDE");
  CaO = nistManager->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
  MgO = nistManager->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
  TiO2 = nistManager->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");

  Air = nistManager->FindOrBuildMaterial("G4_AIR");
  
  Mn2O3 = new G4Material("ManganeseOxide", density = 4.5*g/cm3, ncomponents = 2);
  Mn2O3->AddElement(elMn, natoms=2);
  Mn2O3->AddElement(elO, natoms=3);
  
  //////////////////////////// pH section ////////////////////////////

  // create OH for material
  OH = new G4Material("Hydroxide", density = 1.33*g/cm3, ncomponents = 2);
  OH->AddElement(elO, natoms=1);
  OH->AddElement(elH, natoms=1);
  
  // pH material to add in soil
  pH = new G4Material("pH", density = 1.33*g/cm3, ncomponents = 2);
  pH->AddMaterial(H, fractionmass=fractionMassForH*perCent);
  pH->AddMaterial(OH, fractionmass=fractionMassForOH*perCent);
  
  // G4cout << "H concentration = " << fractionMassForH << G4endl;
  // G4cout << "OH concentration = " << fractionMassForOH << G4endl;
  
  //////////////////////////// end pH section ////////////////////////////
  
  // organic material
  // We assume organic soil components have 1.33 g/cm3 density. This number should be updated with
  // specifric soil type simulation. 4/12/2019

  // ref: http://passel-test.unl.edu/beta/pages/informationmodule.php?idinformationmodule=1130447039&topicorder=5&maxto=10&minto=1
  OrganicMat = new G4Material("Organic", density = 0.8*g/cm3, ncomponents = 7); 
  OrganicMat->AddElement(elC, fractionmass=2.65*perCent);    // variable from 0.7% upto 14%
  OrganicMat->AddElement(elO, fractionmass=15.6*perCent);   // less than 20.6
  OrganicMat->AddElement(elN, fractionmass=79.0*perCent); 
  OrganicMat->AddElement(elK, fractionmass=1.5*perCent);    // 0.3 to 2.5%
  OrganicMat->AddElement(elP, fractionmass=0.6*perCent);    // avarage % of P in soil
  OrganicMat->AddElement(elMg, fractionmass=0.5*perCent);   // 0.05 and 0.5% 
  OrganicMat->AddElement(elS, fractionmass=0.15*perCent);   // avarage

  coatingMaterial = new G4Material("coatingMaterial", density = 0.50*g/cm3, ncomponents = 2);
  coatingMaterial->AddMaterial(H2OVapor, fractionmass=50*perCent);
  coatingMaterial->AddMaterial(H2O, fractionmass=50*perCent);
  
  // vacuum material
  defaultMaterial = new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RPS2023DetectorConstruction::DefineVolumes()
{  
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // Geometry parameters  
    
    // World
    //
    auto worldSizeXY = 2.5 * boxSizeX;
    // auto worldSizeZ  = worldSizeXY; 
    
    auto worldS 
        = new G4Box("World",           // its name
                    worldSizeXY/2, worldSizeXY/2, worldSizeXY/2); // its size
                            
    auto worldLV
        = new G4LogicalVolume(
                    worldS,           // its solid
                    defaultMaterial,  // its material
                    "World");         // its name
                                    
    auto worldPV
        = new G4PVPlacement(
                    0,                // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    worldLV,          // its logical volume                         
                    "World",          // its name
                    0,                // its mother  volume
                    false,            // no boolean operation
                    0,                // copy number
                    fCheckOverlaps);  // checking overlaps 
    // box
    G4Box*  
    box = new G4Box("box",                            // its name
                    boxSizeX,boxSizeY,boxSizeZ);     //its size
                                            
    boxLV = new G4LogicalVolume(box,                  // its solid
                                poreMaterial,      // its material
                                "boxLV");             // its name                                
    boxPV = new G4PVPlacement(0,                      // no rotation
                                G4ThreeVector(),        // at (0,0,0)
                                boxLV,                  // its logical volume
                                "boxPV",                // its name
                                worldLV,                // its mother  volume
                                false,                  // no boolean operation
                                0);                     // copy number

    // grain
    // Define the angles for spheres
    G4double startAnglePhi = 0.0*deg;
    G4double spanningAnglePhi = 360.0*deg;
    G4double startAngleTheta = 0.0*deg;
    G4double spanningAngleTheta = 180.0*deg; //90


    // grain
    auto grainSolid = new G4Sphere("grain", 0, grainSize,
                    startAnglePhi, spanningAnglePhi, 
                    startAngleTheta, spanningAngleTheta);  
    grainLV
        = new G4LogicalVolume(
                    grainSolid,        // its solid
                    grainMaterial,     // its material/chemical composition of soil grain
                    "grainLV");        // its name
                                    
    grainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(),   // at (0,0,0)
                    grainLV,           // its logical volume                         
                    "grainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // surrounding grains
    auto surroundingGrainSolidS = new G4Sphere("surroundingGrain", 0, grainSize,
                    startAnglePhi, spanningAnglePhi, 
                    startAngleTheta, spanningAngleTheta); 

    surroundingGrainLV
            = new G4LogicalVolume(
                        surroundingGrainSolidS,        // its solid
                        grainMaterial,                 // its material/chemical composition of soil grain
                        "surroundingGrainLV");                    // its name

    // right x+
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,0,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps
    // top y+
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,2*grainSize,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps
    // left x-
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,0,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom y-
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,-2*grainSize,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // front z+
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,0,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // front left x- z+
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,0,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // front right x+ z+
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,0,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // back z-
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,0,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // back left x- z-
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,0,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // back right x+ z-
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,0,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom front 
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,-2*grainSize,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom left 
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,-2*grainSize,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom right 
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,-2*grainSize,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom back
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,-2*grainSize,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom front left 
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,-2*grainSize,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom front right 
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,-2*grainSize,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom back left
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,-2*grainSize,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // bottom back right
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,-2*grainSize,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps


    // top left
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,2*grainSize,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top right
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,2*grainSize,0),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top front
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,2*grainSize,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top back
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(0,2*grainSize,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top front left
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,2*grainSize,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top front right
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,2*grainSize,2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top back left
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(-2*grainSize,2*grainSize,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,             // its mother  volume
                    false,             // no boolean operation
                    0,                 // copy number
                    fCheckOverlaps);   // checking overlaps

    // top back right
    surroundingGrainPV = new G4PVPlacement(
                    0,                 // no rotation
                    G4ThreeVector(2*grainSize,2*grainSize,-2*grainSize),   // at (0,0,0)
                    surroundingGrainLV,           // its logical volume                         
                    "surroundingGrainPV",         // its name
                    boxLV,                        // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps

  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto grainVisAtt= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  grainVisAtt->SetVisibility(true);
  grainLV->SetVisAttributes(grainVisAtt);

  auto surroundingGrainVisAtt= new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  surroundingGrainVisAtt->SetVisibility(true);
  surroundingGrainLV->SetVisAttributes(surroundingGrainVisAtt);

  PrintLayersMaterials();
  
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RPS2023DetectorConstruction::PrintLayersMaterials()
{
  G4cout << "\t******* MATERIALS OF EACH LAYER *******" << G4endl;
  G4cout << " grain material : " << grainLV->GetMaterial()->GetName() << G4endl;
  G4cout << " surrounding grain material : " << surroundingGrainLV->GetMaterial()->GetName() << G4endl;
//   G4cout << " pore grain material : " << poreGrainLV->GetMaterial()->GetName() << G4endl;
  G4cout << " pore material : " << boxLV->GetMaterial()->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023DetectorConstruction::updateH2OContent(G4int value)
{
    // update value of variable H2O content
    H2OContent = value;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023DetectorConstruction::setGrainMaterial(G4int value) // value = h20 content
{
    if(value == 5)
    {
        grainMaterial = grainComp1;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }
        
        poreMaterial = poreComp1;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 10)
    {
        grainMaterial = grainComp2;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp2;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 15)
    {
        grainMaterial = grainComp3;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp3;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 20)
    {
        grainMaterial = grainComp4;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp4;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 25)
    {
        grainMaterial = grainComp5;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp5;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 30)
    {
        grainMaterial = grainComp6;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp6;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 35)
    {
        grainMaterial = grainComp7;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp7;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 40)
    {
        grainMaterial = grainComp8;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp8;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 45)
    {
        grainMaterial = grainComp9;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp9;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 50)
    {
        grainMaterial = grainComp10;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp10;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 55)
    {
        grainMaterial = grainComp11;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp11;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 60)
    {
        grainMaterial = grainComp12;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp12;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 65)
    {
        grainMaterial = grainComp13;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp13;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 70)
    {
        grainMaterial = grainComp14;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp14;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 75)
    {
        grainMaterial = grainComp15;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp15;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 80)
    {
        grainMaterial = grainComp16;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp16;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 85)
    {
        grainMaterial = grainComp17;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp17;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 90)
    {
        grainMaterial = grainComp18;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp18;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 95)
    {
        grainMaterial = grainComp19;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp19;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    if(value == 100)
    {
        grainMaterial = grainComp20;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
        
        if(surroundingGrainLV) { surroundingGrainLV->SetMaterial(grainMaterial); }

        poreMaterial = poreComp20;
        if(boxLV) { boxLV->SetMaterial(poreMaterial); }
    }

    // notify run manager that physics is modified
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String RPS2023DetectorConstruction::getNameOfLayer1()
{
    return "grainPV"; 
}

G4String RPS2023DetectorConstruction::getNameOfLayer2()
{
    return "boxPV"; 
}

G4String RPS2023DetectorConstruction::getNameOfLayer3()
{
    return "surroundingGrainPV";
}

G4String RPS2023DetectorConstruction::getNameOfLayer4()
{
    return "poreGrainPV";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023DetectorConstruction::DefineChemicalComps()
{
    // define chemical composition maps
    // use these maps in soil layers
    // you can add or delete any material composition
    // in any map
    // just make sure that summation = 100*perCent

    // refs:
    // organic matter %: https://www.agric.wa.gov.au/measuring-and-assessing-soils/what-soil-organic-carbon#:~:text=Organic%20matter%20makes%20up%20just,biological%20function%20of%20agricultural%20soils.
    // soil grains does not have any pore space i.e. no air

    // 0% water
    chem_composition_0 = {
        {Air, 3*perCent},
        {H, 4*perCent},
        {C, 10*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 28*perCent},
        {Na, 18*perCent},
        {Ca, 15*perCent},
        {Mg, 10*perCent},
        {S, 2*perCent},
        {P, 2*perCent},
        {Si, 1*perCent},
        {Al, 1*perCent},
        {Fe, 1*perCent}
    };

    // 5% water
    chem_composition_1 = {
        {H2O, 5*perCent},
        {Air, 2*perCent},
        {H, 3*perCent},
        {C, 9*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 26*perCent},
        {Na, 18*perCent},
        {Ca, 15*perCent},
        {Mg, 10*perCent},
        {S, 2*perCent},
        {P, 2*perCent},
        {Si, 1*perCent},
        {Al, 1*perCent},
        {Fe, 1*perCent}
    };

    // 10% water
    chem_composition_2 = {
        {H2O, 10*perCent},
        {Air, 2*perCent},
        {H, 2*perCent},
        {C, 8*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 25*perCent},
        {Na, 17*perCent},
        {Ca, 14*perCent},
        {Mg, 10*perCent},
        {S, 2*perCent},
        {P, 2*perCent},
        {Si, 1*perCent},
        {Al, 1*perCent},
        {Fe, 1*perCent}
    };

    // 15% water
    chem_composition_3 = {
        {H2O, 15*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 25*perCent},
        {Na, 15*perCent},
        {Ca, 15*perCent},
        {Mg, 10*perCent},
        {S, 2*perCent},
        {P, 2*perCent},
        {Si, 1*perCent},
        {Al, 1*perCent},
        {Fe, 1*perCent}
    };
    
    // 20% water
    chem_composition_4 = {
        {H2O, 20*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 25*perCent},
        {Na, 14*perCent},
        {Ca, 14*perCent},
        {Mg, 9*perCent},
        {S, 2*perCent},
        {P, 1*perCent},
        {Si, 1*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };
    
    // 25% water
    chem_composition_5 = {
        {H2O, 25*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 25*perCent},
        {Na, 13*perCent},
        {Ca, 13*perCent},
        {Mg, 8*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };  

    // 35% water
    chem_composition_6 = {
        {H2O, 30*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 23*perCent},
        {Na, 12*perCent},
        {Ca, 12*perCent},
        {Mg, 7*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };
    
    // 35% water
    chem_composition_7 = {
        {H2O, 35*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 21*perCent},
        {Na, 11*perCent},
        {Ca, 11*perCent},
        {Mg, 6*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };  

    // 40% water
    chem_composition_8 = {
        {H2O, 40*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 20*perCent},
        {Na, 10*perCent},
        {Ca, 9*perCent},
        {Mg, 5*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 45% water
    chem_composition_9 = {
        {H2O, 45*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 18*perCent},
        {Na, 9*perCent},
        {Ca, 8*perCent},
        {Mg, 4*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 50% water
    chem_composition_10 = {
        {H2O, 50*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 17*perCent},
        {Na, 8*perCent},
        {Ca, 7*perCent},
        {Mg, 2*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 55% water
    chem_composition_11 = {
        {H2O, 55*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 5*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 15*perCent},
        {Na, 7*perCent},
        {Ca, 6*perCent},
        {Mg, 1*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 60% water
    chem_composition_12 = {
        {H2O, 60*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 4*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 13*perCent},
        {Na, 6*perCent},
        {Ca, 5*perCent},
        {Mg, 1*perCent},
        {S, 1*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 65% water
    chem_composition_13 = {
        {H2O, 65*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 4*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 11*perCent},
        {Na, 5*perCent},
        {Ca, 4*perCent},
        {Mg, 0.5*perCent},
        {S, 0.5*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 70% water
    chem_composition_14 = {
        {H2O, 70*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 4*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 9*perCent},
        {Na, 3*perCent},
        {Ca, 3*perCent},
        {Mg, 0.5*perCent},
        {S, 0.5*perCent},
        {P, 0.5*perCent},
        {Si, 0.5*perCent},
        {Al, 0.5*perCent},
        {Fe, 0.5*perCent}
    };

    // 75% water
    chem_composition_15 = {
        {H2O, 75*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 4*perCent},
        {N, 3*perCent},
        {O, 2*perCent},
        {K, 7*perCent},
        {Na, 2*perCent},
        {Ca, 2*perCent},
        {Mg, 0.5*perCent},
        {S, 0.3*perCent},
        {P, 0.3*perCent},
        {Si, 0.3*perCent},
        {Al, 0.3*perCent},
        {Fe, 0.3*perCent}
    };
    
    // 80% water
    chem_composition_16 = {
        {H2O, 80*perCent},
        {Air, 1*perCent},
        {H, 2*perCent},
        {C, 4*perCent},
        {N, 2*perCent},
        {O, 1*perCent},
        {K, 6*perCent},
        {Na, 1*perCent},
        {Ca, 1*perCent},
        {Mg, 0.5*perCent},
        {S, 0.3*perCent},
        {P, 0.3*perCent},
        {Si, 0.3*perCent},
        {Al, 0.3*perCent},
        {Fe, 0.3*perCent}
    };

    // 85% water
    chem_composition_17 = {
        {H2O, 85*perCent},
        {Air, 1*perCent},
        {H, 1*perCent},
        {C, 3*perCent},
        {N, 1*perCent},
        {O, 1*perCent},
        {K, 4*perCent},
        {Na, 1*perCent},
        {Ca, 1*perCent},
        {Mg, 0.5*perCent},
        {S, 0.3*perCent},
        {P, 0.3*perCent},
        {Si, 0.3*perCent},
        {Al, 0.3*perCent},
        {Fe, 0.3*perCent}
    };

    // 90% water
    chem_composition_18 = {
        {H2O, 90*perCent},
        {Air, 1*perCent},
        {H, 1*perCent},
        {C, 3*perCent},
        {N, 1*perCent},
        {O, 1*perCent},
        {K, 1*perCent},
        {Na, 1*perCent},
        {Ca, 1*perCent}
    };

    // 95% water
    chem_composition_19 = {
        {H2O, 95*perCent},
        {H, 1*perCent},
        {C, 3*perCent},
        {N, 1*perCent}
    };

    // 100% water
    chem_composition_20 = {
        {H2O, 100*perCent},
    };
}

void RPS2023DetectorConstruction::fillGrainWithChemComps()
{   
    // grain comp1 has particle density = 2.66
    // ref: http://passel-test.unl.edu/beta/pages/informationmodule.php?idinformationmodule=1130447039&topicorder=5&maxto=10&minto=1
    grainComp0 = new G4Material("0W", density = 0.6*g/cm3, ncomponents=chem_composition_0.size());
    for (it = chem_composition_0.begin(); it != chem_composition_0.end(); it++) {
            grainComp0->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp1 = new G4Material("5W", density = 0.6*g/cm3, ncomponents=chem_composition_1.size());
    for (it = chem_composition_1.begin(); it != chem_composition_1.end(); it++) {
            grainComp1->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp2 = new G4Material("10W", density = 0.6*g/cm3, ncomponents=chem_composition_2.size());
    for (it = chem_composition_2.begin(); it != chem_composition_2.end(); it++) {
            grainComp2->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp3 = new G4Material("15W", density = 0.6*g/cm3, ncomponents=chem_composition_3.size());
    for (it = chem_composition_3.begin(); it != chem_composition_3.end(); it++) {
            grainComp3->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp4 = new G4Material("20W", density = 0.6*g/cm3, ncomponents=chem_composition_4.size());
    for (it = chem_composition_4.begin(); it != chem_composition_4.end(); it++) {
            grainComp4->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp5 = new G4Material("25W", density = 0.6*g/cm3, ncomponents=chem_composition_5.size());
    for (it = chem_composition_5.begin(); it != chem_composition_5.end(); it++) {
            grainComp5->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp6 = new G4Material("30W", density = 0.6*g/cm3, ncomponents=chem_composition_6.size());
    for (it = chem_composition_6.begin(); it != chem_composition_6.end(); it++) {
            grainComp6->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp7 = new G4Material("35W", density = 0.6*g/cm3, ncomponents=chem_composition_7.size());
    for (it = chem_composition_7.begin(); it != chem_composition_7.end(); it++) {
            grainComp7->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp8 = new G4Material("40W", density = 0.6*g/cm3, ncomponents=chem_composition_8.size());
    for (it = chem_composition_8.begin(); it != chem_composition_8.end(); it++) {
            grainComp8->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp9 = new G4Material("45W", density = 0.6*g/cm3, ncomponents=chem_composition_9.size());
    for (it = chem_composition_9.begin(); it != chem_composition_9.end(); it++) {
            grainComp9->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp10 = new G4Material("50W", density = 0.6*g/cm3, ncomponents=chem_composition_10.size());
    for (it = chem_composition_10.begin(); it != chem_composition_10.end(); it++) {
            grainComp10->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp11 = new G4Material("55W", density = 0.6*g/cm3, ncomponents=chem_composition_11.size());
    for (it = chem_composition_11.begin(); it != chem_composition_11.end(); it++) {
            grainComp11->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp12 = new G4Material("60W", density = 0.6*g/cm3, ncomponents=chem_composition_12.size());
    for (it = chem_composition_12.begin(); it != chem_composition_12.end(); it++) {
            grainComp12->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp13 = new G4Material("65W", density = 0.6*g/cm3, ncomponents=chem_composition_13.size());
    for (it = chem_composition_13.begin(); it != chem_composition_13.end(); it++) {
            grainComp13->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp14 = new G4Material("70W", density = 0.6*g/cm3, ncomponents=chem_composition_14.size());
    for (it = chem_composition_14.begin(); it != chem_composition_14.end(); it++) {
            grainComp14->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp15 = new G4Material("75W", density = 0.6*g/cm3, ncomponents=chem_composition_15.size());
    for (it = chem_composition_15.begin(); it != chem_composition_15.end(); it++) {
            grainComp15->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp16 = new G4Material("80W", density = 0.6*g/cm3, ncomponents=chem_composition_16.size());
    for (it = chem_composition_16.begin(); it != chem_composition_16.end(); it++) {
            grainComp16->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp17 = new G4Material("85W", density = 0.6*g/cm3, ncomponents=chem_composition_17.size());
    for (it = chem_composition_17.begin(); it != chem_composition_17.end(); it++) {
            grainComp17->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp18 = new G4Material("90W", density = 0.6*g/cm3, ncomponents=chem_composition_18.size());
    for (it = chem_composition_18.begin(); it != chem_composition_18.end(); it++) {
            grainComp18->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp19 = new G4Material("95W", density = 0.6*g/cm3, ncomponents=chem_composition_19.size());
    for (it = chem_composition_19.begin(); it != chem_composition_19.end(); it++) {
            grainComp19->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp20 = new G4Material("100W", density = 0.6*g/cm3, ncomponents=chem_composition_20.size());
    for (it = chem_composition_20.begin(); it != chem_composition_20.end(); it++) {
            grainComp20->AddMaterial(it->first, fractionmass=it->second);
    }

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RPS2023DetectorConstruction::DefinePoreChemicalComps()
{   
    // pore_chem_composition_0 = {
    //     {H2O, 1*perCent}, // mass fractions cannot be 0
    //     {H2OVapor, 99*perCent}
    // };

    pore_chem_composition_0 = {
        {Air, 100*perCent}
    };
    
    pore_chem_composition_1 = {
        {H2O, 5*perCent},
        {Air, 95*perCent}
    };
    
    pore_chem_composition_2 = {
        {H2O, 10*perCent},
        {Air, 90*perCent}
    };

    pore_chem_composition_3 = {
        {H2O, 15*perCent},
        {Air, 85*perCent}
    };

    pore_chem_composition_4 = {
        {H2O, 20*perCent},
        {Air, 80*perCent}
    };

    pore_chem_composition_5 = {
        {H2O, 25*perCent},
        {Air, 75*perCent}
    };

    pore_chem_composition_6 = {
        {H2O, 30*perCent},
        {Air, 70*perCent}
    };

    pore_chem_composition_7 = {
        {H2O, 35*perCent},
        {Air, 65*perCent}
    };

    pore_chem_composition_8 = {
        {H2O, 40*perCent},
        {Air, 60*perCent}
    };

    pore_chem_composition_9 = {
        {H2O, 45*perCent},
        {Air, 55*perCent}
    };

    pore_chem_composition_10 = {
        {H2O, 50*perCent},
        {Air, 50*perCent}
    };

    pore_chem_composition_11 = {
        {H2O, 55*perCent},
        {Air, 45*perCent}
    };

    pore_chem_composition_12 = {
        {H2O, 60*perCent},
        {Air, 40*perCent}
    };

    pore_chem_composition_13 = {
        {H2O, 65*perCent},
        {Air, 35*perCent}
    };

    pore_chem_composition_14 = {
        {H2O, 70*perCent},
        {Air, 30*perCent}
    };

    pore_chem_composition_15 = {
        {H2O, 75*perCent},
        {Air, 25*perCent}
    };

    pore_chem_composition_16 = {
        {H2O, 80*perCent},
        {Air, 20*perCent}
    };

    pore_chem_composition_17 = {
        {H2O, 85*perCent},
        {Air, 15*perCent}
    };

    pore_chem_composition_18 = {
        {H2O, 90*perCent},
        {Air, 10*perCent}
    };

    pore_chem_composition_19 = {
        {H2O, 95*perCent},
        {Air, 5*perCent}
    };

    pore_chem_composition_20 = {
        {H2O, 100*perCent}
    };
}

void RPS2023DetectorConstruction::fillPoresWithChemComps()
{
    poreComp0 = new G4Material("P0W", density = 0.00120479*g/cm3, ncomponents=pore_chem_composition_0.size());
    for (it = pore_chem_composition_0.begin(); it != pore_chem_composition_0.end(); it++) {
            poreComp0->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp1 = new G4Material("P5W", density = 0.051*g/cm3, ncomponents=pore_chem_composition_1.size());
    for (it = pore_chem_composition_1.begin(); it != pore_chem_composition_1.end(); it++) {
            poreComp1->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp2 = new G4Material("P10W", density = 0.10*g/cm3, ncomponents=pore_chem_composition_2.size());
    for (it = pore_chem_composition_2.begin(); it != pore_chem_composition_2.end(); it++) {
            poreComp2->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp3 = new G4Material("P15W", density = 0.15*g/cm3, ncomponents=pore_chem_composition_3.size());
    for (it = pore_chem_composition_3.begin(); it != pore_chem_composition_3.end(); it++) {
            poreComp3->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp4 = new G4Material("P20W", density = 0.20*g/cm3, ncomponents=pore_chem_composition_4.size());
    for (it = pore_chem_composition_4.begin(); it != pore_chem_composition_4.end(); it++) {
            poreComp4->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp5 = new G4Material("P25W", density = 0.25*g/cm3, ncomponents=pore_chem_composition_5.size());
    for (it = pore_chem_composition_5.begin(); it != pore_chem_composition_5.end(); it++) {
            poreComp5->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp6 = new G4Material("P30W", density = 0.30*g/cm3, ncomponents=pore_chem_composition_6.size());
    for (it = pore_chem_composition_6.begin(); it != pore_chem_composition_6.end(); it++) {
            poreComp6->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp7 = new G4Material("P35W", density = 0.35*g/cm3, ncomponents=pore_chem_composition_7.size());
    for (it = pore_chem_composition_7.begin(); it != pore_chem_composition_7.end(); it++) {
            poreComp7->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp8 = new G4Material("P40W", density = 0.40*g/cm3, ncomponents=pore_chem_composition_8.size());
    for (it = pore_chem_composition_8.begin(); it != pore_chem_composition_8.end(); it++) {
            poreComp8->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp9 = new G4Material("P45W", density = 0.45*g/cm3, ncomponents=pore_chem_composition_9.size());
    for (it = pore_chem_composition_9.begin(); it != pore_chem_composition_9.end(); it++) {
            poreComp9->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp10 = new G4Material("P50W", density = 0.50*g/cm3, ncomponents=pore_chem_composition_10.size());
    for (it = pore_chem_composition_10.begin(); it != pore_chem_composition_10.end(); it++) {
            poreComp10->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp11 = new G4Material("P55W", density = 0.55*g/cm3, ncomponents=pore_chem_composition_11.size());
    for (it = pore_chem_composition_11.begin(); it != pore_chem_composition_11.end(); it++) {
            poreComp11->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp12 = new G4Material("P60W", density = 0.60*g/cm3, ncomponents=pore_chem_composition_12.size());
    for (it = pore_chem_composition_12.begin(); it != pore_chem_composition_12.end(); it++) {
            poreComp12->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp13 = new G4Material("P65W", density = 0.65*g/cm3, ncomponents=pore_chem_composition_13.size());
    for (it = pore_chem_composition_13.begin(); it != pore_chem_composition_13.end(); it++) {
            poreComp13->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp14 = new G4Material("P70W", density = 0.70*g/cm3, ncomponents=pore_chem_composition_14.size());
    for (it = pore_chem_composition_14.begin(); it != pore_chem_composition_14.end(); it++) {
            poreComp14->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp15 = new G4Material("P75W", density = 0.75*g/cm3, ncomponents=pore_chem_composition_15.size());
    for (it = pore_chem_composition_15.begin(); it != pore_chem_composition_15.end(); it++) {
            poreComp15->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp16 = new G4Material("P80W", density = 0.80*g/cm3, ncomponents=pore_chem_composition_16.size());
    for (it = pore_chem_composition_16.begin(); it != pore_chem_composition_16.end(); it++) {
            poreComp16->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp17 = new G4Material("P85W", density = 0.85*g/cm3, ncomponents=pore_chem_composition_17.size());
    for (it = pore_chem_composition_17.begin(); it != pore_chem_composition_17.end(); it++) {
            poreComp17->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp18 = new G4Material("P90W", density = 0.90*g/cm3, ncomponents=pore_chem_composition_18.size());
    for (it = pore_chem_composition_18.begin(); it != pore_chem_composition_18.end(); it++) {
            poreComp18->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp19 = new G4Material("P95W", density = 0.95*g/cm3, ncomponents=pore_chem_composition_19.size());
    for (it = pore_chem_composition_19.begin(); it != pore_chem_composition_19.end(); it++) {
            poreComp19->AddMaterial(it->first, fractionmass=it->second);
    }

    poreComp20 = new G4Material("P100W", density = 1.00*g/cm3, ncomponents=pore_chem_composition_20.size());
    for (it = pore_chem_composition_20.begin(); it != pore_chem_composition_20.end(); it++) {
            poreComp20->AddMaterial(it->first, fractionmass=it->second);
    }
}

