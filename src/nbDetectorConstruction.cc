// This file will make use nbDetectorConstruction.hh

#include "nbDetectorConstruction.hh"

#include "nbDetectorMessenger.hh"

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
nbDetectorConstruction::nbDetectorConstruction()
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
  fdetectorMessenger = new nbDetectorMessenger(this);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// destructor
nbDetectorConstruction::~nbDetectorConstruction()
{ 
    delete fdetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// member function of header file
G4VPhysicalVolume* nbDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();

  G4cout << "soil grain size : " << grainSize << G4endl; 

  // define chemical composition maps
  DefineChemicalComps();
  fillGrainWithChemComps();
  
  // initialize grain material
  grainMaterial = grainComp1;

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::DefineMaterials()
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

G4VPhysicalVolume* nbDetectorConstruction::DefineVolumes()
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
                              defaultMaterial,      // its material
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

  
  // coating
//   auto coatingSolid = new G4Sphere("coating", 0, totalGrainSize,
// 			       startAnglePhi, spanningAnglePhi, 
// 			       startAngleTheta, spanningAngleTheta);

//   coatingLV
//     = new G4LogicalVolume(
//                  coatingSolid,        // its solid
//                  coatingMaterial,       // its material/chemical composition of soil grain
//                  "coatingLV");        // its name
                                   
//   coatingPV = new G4PVPlacement(
//                  0,                 // no rotation
//                  G4ThreeVector(),   // at (0,0,0)
//                  coatingLV,           // its logical volume                         
//                  "coatingPV",         // its name
//                  boxLV,             // its mother  volume
//                  false,             // no boolean operation
//                  0,                 // copy number
//                  fCheckOverlaps);   // checking overlaps


  // grain
  auto grainSolid = new G4Sphere("grain", 0, grainSize,
			       startAnglePhi, spanningAnglePhi, 
			       startAngleTheta, spanningAngleTheta);  
  grainLV
    = new G4LogicalVolume(
                 grainSolid,        // its solid
                 grainMaterial,        // its material/chemical composition of soil grain
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

  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto grainVisAtt= new G4VisAttributes(G4Colour(0.45,0.25,0.0));
  grainVisAtt->SetVisibility(true);
  grainLV->SetVisAttributes(grainVisAtt);

  PrintLayersMaterials();
  
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void nbDetectorConstruction::PrintLayersMaterials()
{
  G4cout << "\t******* MATERIALS OF EACH LAYER *******" << G4endl;
  G4cout << " grain material : " << grainLV->GetMaterial()->GetName() << G4endl;
  // G4cout << " layer 1 material : " << shellLV_1->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 2 material : " << shellLV_2->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 3 material : " << shellLV_3->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 4 material : " << shellLV_4->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 5 material : " << shellLV_5->GetMaterial()->GetName() << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::updateH2OContent(G4int value)
{
    // update value of variable H2O content
    H2OContent = value;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::setGrainMaterial(G4int value) // value = h20 content
{
    if(value == 15)
    {
        grainMaterial = grainComp1;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 20)
    {
        grainMaterial = grainComp2;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 25)
    {
        grainMaterial = grainComp3;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 30)
    {
        grainMaterial = grainComp4;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 35)
    {
        grainMaterial = grainComp5;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 40)
    {
        grainMaterial = grainComp6;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 45)
    {
        grainMaterial = grainComp7;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 50)
    {
        grainMaterial = grainComp8;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 55)
    {
        grainMaterial = grainComp9;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 60)
    {
        grainMaterial = grainComp10;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 65)
    {
        grainMaterial = grainComp11;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 70)
    {
        grainMaterial = grainComp12;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    if(value == 75)
    {
        grainMaterial = grainComp13;
        if(grainLV) { grainLV->SetMaterial(grainMaterial); }
    }

    // notify run manager that physics is modified
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String nbDetectorConstruction::getNameOfLayer1()
{
    return "grainPV"; // this is first layer beneath the surface
}

G4String nbDetectorConstruction::getNameOfLayer2()
{
    return "boxPV"; // the second layer below that
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::DefineChemicalComps()
{
    // define chemical composition maps
    // use these maps in soil layers
    // you can add or delete any material composition
    // in any map
    // just make sure that summation = 100*perCent

    // refs:
    // organic matter %: https://www.agric.wa.gov.au/measuring-and-assessing-soils/what-soil-organic-carbon#:~:text=Organic%20matter%20makes%20up%20just,biological%20function%20of%20agricultural%20soils.
    // soil grains does not have any pore space i.e. no air

    // first layer composition
    chem_composition_1 = {
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
    
    // increased H2O and decreased quartz and organic matter
    chem_composition_2 = {
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
    
    // increased H20 and decreased al2O3, organic, Mn2O3, Al2O3 and quartz
    chem_composition_3 = {
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

    chem_composition_4 = {
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
    
    chem_composition_5 = {
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
    
    chem_composition_6 = {
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

    chem_composition_7 = {
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

    chem_composition_8 = {
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

    chem_composition_9 = {
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

    chem_composition_10 = {
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

    chem_composition_11 = {
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

    chem_composition_12 = {
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

    chem_composition_13 = {
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
}

void nbDetectorConstruction::fillGrainWithChemComps()
{   
    // grain comp1 has particle density = 2.66
    // ref: http://passel-test.unl.edu/beta/pages/informationmodule.php?idinformationmodule=1130447039&topicorder=5&maxto=10&minto=1
    grainComp1 = new G4Material("15W", density = 0.6*g/cm3, ncomponents=chem_composition_1.size());
    for (it = chem_composition_1.begin(); it != chem_composition_1.end(); it++) {
            grainComp1->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp2 = new G4Material("20W", density = 0.6*g/cm3, ncomponents=chem_composition_2.size());
    for (it = chem_composition_2.begin(); it != chem_composition_2.end(); it++) {
            grainComp2->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp3 = new G4Material("25W", density = 0.6*g/cm3, ncomponents=chem_composition_3.size());
    for (it = chem_composition_3.begin(); it != chem_composition_3.end(); it++) {
            grainComp3->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp4 = new G4Material("30W", density = 0.6*g/cm3, ncomponents=chem_composition_4.size());
    for (it = chem_composition_4.begin(); it != chem_composition_4.end(); it++) {
            grainComp4->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp5 = new G4Material("35W", density = 0.6*g/cm3, ncomponents=chem_composition_5.size());
    for (it = chem_composition_5.begin(); it != chem_composition_5.end(); it++) {
            grainComp5->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp6 = new G4Material("40W", density = 0.6*g/cm3, ncomponents=chem_composition_6.size());
    for (it = chem_composition_6.begin(); it != chem_composition_6.end(); it++) {
            grainComp6->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp7 = new G4Material("45W", density = 0.6*g/cm3, ncomponents=chem_composition_7.size());
    for (it = chem_composition_7.begin(); it != chem_composition_7.end(); it++) {
            grainComp7->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp8 = new G4Material("50W", density = 0.6*g/cm3, ncomponents=chem_composition_8.size());
    for (it = chem_composition_8.begin(); it != chem_composition_8.end(); it++) {
            grainComp8->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp9 = new G4Material("55W", density = 0.6*g/cm3, ncomponents=chem_composition_9.size());
    for (it = chem_composition_9.begin(); it != chem_composition_9.end(); it++) {
            grainComp9->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp10 = new G4Material("60W", density = 0.6*g/cm3, ncomponents=chem_composition_10.size());
    for (it = chem_composition_10.begin(); it != chem_composition_10.end(); it++) {
            grainComp10->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp11 = new G4Material("65W", density = 0.6*g/cm3, ncomponents=chem_composition_11.size());
    for (it = chem_composition_11.begin(); it != chem_composition_11.end(); it++) {
            grainComp11->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp12 = new G4Material("70W", density = 0.6*g/cm3, ncomponents=chem_composition_12.size());
    for (it = chem_composition_12.begin(); it != chem_composition_12.end(); it++) {
            grainComp12->AddMaterial(it->first, fractionmass=it->second);
    }

    grainComp13 = new G4Material("75W", density = 0.6*g/cm3, ncomponents=chem_composition_13.size());
    for (it = chem_composition_13.begin(); it != chem_composition_13.end(); it++) {
            grainComp13->AddMaterial(it->first, fractionmass=it->second);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......