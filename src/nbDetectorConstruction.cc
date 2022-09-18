// Neutron ball detector construction
//
// Created on 10/13/2021: hexc, Mayur, Tien, Weisen, Austin
// neutron ball detector is constructed as a shell with a predeinfed inner and outer radia.
//
// Updated on 10/20/2021: hexc, Mayur, Tien and Weisen
//    Add code for reading in the detector configuration file
//
// Updated on 11/10/2021: hexc, Mayur, Tien, Jarvious
//   Add soil material with vaiable water contents
//   Make multiple soil layers
//
// Updated on 22 Aug, 2022
// commented unused variables and added some more comments


// This file will make use nbDetectorConstruction.hh

#include "nbDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
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

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <string>
#include <stdio.h>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// set messenger value
G4ThreadLocal G4GlobalMagFieldMessenger* nbDetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// constructor to initialization values of shellPV, overlap checks 
// set data member values from config.txt file
nbDetectorConstruction::nbDetectorConstruction()
 : G4VUserDetectorConstruction(),
   shellPV(nullptr),
   fCheckOverlaps(true)
{

  std::string line;
  // char varName[10];
  G4double inputVal[6];
  std::ifstream nbConfig("config.txt");
  G4int i = 0;
  if (nbConfig.is_open())
  {
	while (getline(nbConfig, line))
	{
	  char* token = strtok((char *)line.c_str(), ",\t");
	  int num = 0;
	  while (token != NULL)
	    {
	      G4cout << token << G4endl;
	      if (num++ != 0) inputVal[i++] = atoi(token);
	      token = strtok(NULL,",\t");
		/*
		  G4cout << line << G4endl;
		  strtok(line.c_str(), " ,\t");
		  inputVal[i++] = atoi(strtok(NULL, " ,\t"));
		*/
	    }
	}
  }

  inner_r = inputVal[0]*cm;
  outer_r = inputVal[1]*cm;
  matType = (G4int)inputVal[2];
  matType_1 = (G4int)inputVal[3];
  matType_2 = (G4int)inputVal[4];
  matType_3 = (G4int)inputVal[5];
  G4cout << "inner_r: " << inner_r << "   outer_r: " << outer_r << "   matType: " << matType << G4endl;
  G4cout << "inner_r: " << inner_r << "   outer_r: " << outer_r << "   matType_1: " << matType_1 << G4endl;
  G4cout << "inner_r: " << inner_r << "   outer_r: " << outer_r << "   matType_2: " << matType_2 << G4endl;
  G4cout << "inner_r: " << inner_r << "   outer_r: " << outer_r << "   matType_3: " << matType_3 << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// destructor
nbDetectorConstruction::~nbDetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// member function of header file
G4VPhysicalVolume* nbDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::DefineMaterials()
{ 
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Liquid argon material
  // new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  //        // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* nbDetectorConstruction::DefineVolumes()
{
  // some variables to create user defined materials
  G4String symbol;              // periodic table symbol        
  G4double a;                   // mass of a mole;
  G4double z;                   // z=mean number of protons;  
  G4double density;             // density of material
  G4int ncomponents, natoms;    // parameters
  G4double fractionmass;        // parameters    

  // instance of NIST manager which contains inbuilt materials and elements
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  
  // define pure nist materials
  G4Material* H = man->FindOrBuildMaterial("G4_H");
 
  // define nist compounds
  G4Material* H2O  = man->FindOrBuildMaterial("G4_WATER");
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV); // you may or may not include this line, overwrite computed meanExcitationEnergy with ICRU recommended value 
  
  G4Material* SiO2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Al2O3 = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  G4Material* Fe2O3 = man->FindOrBuildMaterial("G4_FERRIC_OXIDE");
  G4Material* CaO = man->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
  G4Material* MgO = man->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
  G4Material* TiO2 = man->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
  
  // define elements
  G4Element* elO = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* elH = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* elC  = man->FindOrBuildElement("C");
  G4Element* elN  = man->FindOrBuildElement("N");
  
  // create OH for material
  G4Material* OH = new G4Material("Hydroxide", density = 1.33*g/cm3, ncomponents = 2);
  OH->AddElement(elO, natoms=1);
  OH->AddElement(elH, natoms=1);
  
  // pH material to add in soil
  G4Material* pH = new G4Material("pH", density = 1.33*g/cm3, ncomponents = 2);
  pH->AddMaterial(H, fractionmass=80.0*perCent);
  pH->AddMaterial(OH, fractionmass=20.0*perCent);
  
  // organic material
  // We assume organic soil components have 1.33 g/cm3 density. This number should be updated with
  // specifric soil type simulation. 4/12/2019
  G4Material* OrganicMat = new G4Material("Organic", density = 1.33*g/cm3, ncomponents = 4);
  OrganicMat->AddElement(elC, natoms=50);
  OrganicMat->AddElement(elO, natoms=42);
  OrganicMat->AddElement(elH, natoms=5);
  OrganicMat->AddElement(elN, natoms=3);
  
  
  // create soil layers
  // Based on: http://gfnun.unal.edu.co/fileadmin/content/gruposdeinvestigacion/fisicanuclear/Tesis/DanielAndrade_TG.pdf
  // layer 1
  auto soilOne = new G4Material("DrySoil", density = 0.6*g/cm3, ncomponents=8);
  soilOne->AddMaterial(SiO2, fractionmass=46.3*perCent);   
  soilOne->AddMaterial(pH, fractionmass=15.0*perCent);      //pH
  soilOne->AddMaterial(Al2O3, fractionmass=13.0*perCent);
  soilOne->AddMaterial(Fe2O3, fractionmass=2.5*perCent);
  soilOne->AddMaterial(CaO, fractionmass=1.6*perCent);
  soilOne->AddMaterial(MgO, fractionmass=0.7*perCent);
  soilOne->AddMaterial(TiO2, fractionmass=0.6*perCent);
  soilOne->AddMaterial(OrganicMat, fractionmass=20.3*perCent);
 
  // layer 2
  // 10% moisture content
  auto soilOne10W = new G4Material("DrySoil10W", density = 0.6*g/cm3, ncomponents=9);
  soilOne10W->AddMaterial(SiO2, fractionmass=40.17*perCent);
  soilOne10W->AddMaterial(pH, fractionmass=15.0*perCent);          //pH
  soilOne10W->AddMaterial(Al2O3, fractionmass=11.7*perCent);
  soilOne10W->AddMaterial(Fe2O3, fractionmass=2.25*perCent);
  soilOne10W->AddMaterial(CaO, fractionmass=1.44*perCent);
  soilOne10W->AddMaterial(MgO, fractionmass=0.63*perCent);
  soilOne10W->AddMaterial(TiO2, fractionmass=0.54*perCent);
  soilOne10W->AddMaterial(OrganicMat, fractionmass=18.27*perCent);
  soilOne10W->AddMaterial(H2O, fractionmass=10.0*perCent);
  
  // layer 3
  // 20% moisture content. Need new density value?
  auto soilOne20W = new G4Material("DrySoil20W", density = 0.6*g/cm3, ncomponents=9);
  soilOne20W->AddMaterial(SiO2, fractionmass=34.04*perCent);
  soilOne20W->AddMaterial(pH, fractionmass=15.0*perCent);          // pH
  soilOne20W->AddMaterial(Al2O3, fractionmass=10.4*perCent);
  soilOne20W->AddMaterial(Fe2O3, fractionmass=2.0*perCent);
  soilOne20W->AddMaterial(CaO, fractionmass=1.28*perCent);
  soilOne20W->AddMaterial(MgO, fractionmass=0.56*perCent);
  soilOne20W->AddMaterial(TiO2, fractionmass=0.48*perCent);
  soilOne20W->AddMaterial(OrganicMat, fractionmass=16.24*perCent);
  soilOne20W->AddMaterial(H2O, fractionmass=20.0*perCent);
 
  // layer 4
  // 30% moisture content. Need new density value?
  auto soilOne30W = new G4Material("DrySoil30W", density = 0.6*g/cm3, ncomponents=9);
  soilOne30W->AddMaterial(SiO2, fractionmass=27.91*perCent);
  soilOne30W->AddMaterial(pH, fractionmass=15.0*perCent);
  soilOne30W->AddMaterial(Al2O3, fractionmass=9.1*perCent);
  soilOne30W->AddMaterial(Fe2O3, fractionmass=1.75*perCent);
  soilOne30W->AddMaterial(CaO, fractionmass=1.12*perCent);
  soilOne30W->AddMaterial(MgO, fractionmass=0.49*perCent);
  soilOne30W->AddMaterial(TiO2, fractionmass=0.42*perCent);
  soilOne30W->AddMaterial(OrganicMat, fractionmass=14.21*perCent);
  soilOne30W->AddMaterial(H2O, fractionmass=30.0*perCent);
    
  // layer 5 (inner most)
  // 40% moisture content. Need new density value?
  auto soilOne40W = new G4Material("DrySoil40W", density = 0.6*g/cm3, ncomponents=9);
  soilOne40W->AddMaterial(SiO2, fractionmass=17.91*perCent);
  soilOne40W->AddMaterial(pH, fractionmass=15.0*perCent);
  soilOne40W->AddMaterial(Al2O3, fractionmass=9.1*perCent);
  soilOne40W->AddMaterial(Fe2O3, fractionmass=1.75*perCent);
  soilOne40W->AddMaterial(CaO, fractionmass=1.12*perCent);
  soilOne40W->AddMaterial(MgO, fractionmass=0.49*perCent);
  soilOne40W->AddMaterial(TiO2, fractionmass=0.42*perCent);
  soilOne40W->AddMaterial(OrganicMat, fractionmass=14.21*perCent);
  soilOne40W->AddMaterial(H2O, fractionmass=40.0*perCent);
  
  // Set materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  
  // assign constructed layers to detector materials which will be inserted in logical volumes
  shellMaterial = soilOne;
  shellMaterial_1 = soilOne10W;
  shellMaterial_2 = soilOne20W;
  shellMaterial_3 = soilOne30W;
  shellMaterial_4 = soilOne40W;

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
 
  // Geometry parameters  


  // World
  //
  auto worldSizeXY = 2.5 * outer_r;
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
  
  //                               
  // Make neutron ball detector (i.e., multiple layers)
  //

  // Define the angles for spheres
  G4double startAnglePhi = 0.0*deg;
  G4double spanningAnglePhi = 360.0*deg;
  G4double startAngleTheta = 0.0*deg;
  G4double spanningAngleTheta = 180.0*deg; //90

  // Outer shell
  auto solidShell = new G4Sphere("Shell", inner_r, outer_r,
			       startAnglePhi, spanningAnglePhi, 
			       startAngleTheta, spanningAngleTheta);  
  auto shellLV
    = new G4LogicalVolume(
                 solidShell,     // its solid
                 shellMaterial,  // its material
                 "shellLV");   // its name
                                   
  shellPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 shellLV,          // its logical volume                         
                 "shellPV",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  // 1st layer (A/Ap 0-6 inches)
  auto solidShell_1 = new G4Sphere("Shell_1", 1.001*inner_r,0.94*outer_r,
				 startAnglePhi, spanningAnglePhi, 
				 startAngleTheta, spanningAngleTheta);
  
  auto shellLV_1
    = new G4LogicalVolume(
                 solidShell_1,     // its solid
                 shellMaterial_1,  // its material
                 "shellLV_1");   // its name
                                   
  shellPV_1 = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 shellLV_1,          // its logical volume                         
                 "shellPV_1",    // its name
                 shellLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // 2nd layer (E 6-9 inches - 100-9)
  auto solidShell_2 = new G4Sphere("Shell_2", 1.002*inner_r, 0.91*outer_r,
				   startAnglePhi, spanningAnglePhi, 
				   startAngleTheta, spanningAngleTheta);
  
  auto shellLV_2
    = new G4LogicalVolume(
			  solidShell_2,     // its solid
			  shellMaterial_2,  // its material
			  "shellLV_2");   // its name
  
  shellPV_2 = new G4PVPlacement(
				0,                // no rotation
				G4ThreeVector(),  // at (0,0,0)
				shellLV_2,          // its logical volume                         
				"shellPV_2",    // its name
				shellLV_1,          // its mother  volume
				false,            // no boolean operation
				0,                // copy number
				fCheckOverlaps);  // checking overlaps 
  
  
  // 3rd layer (Be, Bt, Bc 9-53 inches - 100-53 = 47)
  auto solidShell_3 = new G4Sphere("Shell_3", 1.003*inner_r, 0.47*outer_r,
				   startAnglePhi, spanningAnglePhi, 
				   startAngleTheta, spanningAngleTheta);
  
  auto shellLV_3
    = new G4LogicalVolume(
			  solidShell_3,     // its solid
			  shellMaterial_3,  // its material
			  "shellLV_3");   // its name
  
  shellPV_3 = new G4PVPlacement(
				0,                // no rotation
				G4ThreeVector(),  // at (0,0,0)
				shellLV_3,          // its logical volume                         
				"shellPV_3",    // its name
				shellLV_2,          // its mother  volume
				false,            // no boolean operation
				0,                // copy number
				fCheckOverlaps);  // checking overlaps 
  
  
  // 4th layer (C 53-80 inches - 100-80 = 20)
  auto solidShell_4 = new G4Sphere("Shell_4", 1.003*inner_r, 0.20*outer_r,
                                   startAnglePhi, spanningAnglePhi, 
                                   startAngleTheta, spanningAngleTheta);
  
  auto shellLV_4
    = new G4LogicalVolume(
                          solidShell_4,     // its solid
                          shellMaterial_4,  // its material
                          "shellLV_4");   // its name
  
  shellPV_4 = new G4PVPlacement(
                                0,                // no rotation
                                G4ThreeVector(),  // at (0,0,0)
                                shellLV_4,          // its logical volume
                                "shellPV_4",    // its name
                                shellLV_3,          // its mother  volume
                                false,            // no boolean operation
                                0,                // copy number
                                fCheckOverlaps);  // checking overlaps
  
  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The shell dimension:  inner_r " << inner_r  << "   outer_r " 
    << outer_r  << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleShellVisAtt= new G4VisAttributes(G4Colour(0.3,0.4,0.5));
  simpleShellVisAtt->SetVisibility(true);
  shellLV->SetVisAttributes(simpleShellVisAtt);

  auto simpleShellVisAtt_1= new G4VisAttributes(G4Colour(0.3,0.4,0.1));
  simpleShellVisAtt_1->SetVisibility(true);
  shellLV_1->SetVisAttributes(simpleShellVisAtt_1);

  auto simpleShellVisAtt_2= new G4VisAttributes(G4Colour(0.3,0.3,0.2));
  simpleShellVisAtt_2->SetVisibility(true);
  shellLV_2->SetVisAttributes(simpleShellVisAtt_2);

  auto simpleShellVisAtt_3= new G4VisAttributes(G4Colour(0.3,0.1,0.2));
  simpleShellVisAtt_3->SetVisibility(true);
  shellLV_3->SetVisAttributes(simpleShellVisAtt_3);
  
  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void nbDetectorConstruction::ConstructSDandField()
//{ 
//  // Create global magnetic field messenger.
//  // Uniform magnetic field is then created automatically if
//  // the field value is not zero.
//  G4ThreeVector fieldValue;
//  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
//  fMagFieldMessenger->SetVerboseLevel(1);
//  
//  // Register the field messenger for deleting
//  G4AutoDelete::Register(fMagFieldMessenger);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String nbDetectorConstruction::getNameOfLayer1()
{
    return "shellPV";
}

G4String nbDetectorConstruction::getNameOfLayer2()
{
    return "shellPV_1";
}

G4String nbDetectorConstruction::getNameOfLayer3()
{
    return "shellPV_2";
}

G4String nbDetectorConstruction::getNameOfLayer4()
{
    return "shellPV_3";
}

G4String nbDetectorConstruction::getNameOfLayer5()
{
    return "shellPV_4";
}

G4String nbDetectorConstruction::getNameOfLayer6()
{
    return "World";
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
