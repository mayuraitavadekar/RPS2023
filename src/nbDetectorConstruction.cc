// This file will make use nbDetectorConstruction.hh

#include "nbDetectorConstruction.hh"

#include "nbDetectorMessenger.hh"

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

  // this is core material data map
  // all the materials you defined in DefineMaterials() should be
  // mentioned here
  StringToMaterialMapper = {
        {"Fe", Fe},
        {"Fe2O3", Fe2O3},
        {"Al2O3", Al2O3},
        {"Mn", Mn},
        {"Mn2O3", Mn2O3},
        {"Ra", Ra},
        {"Air", Air},
        {"H2O", H2O},
        {"OH",OH},
        {"MgO", MgO},
        {"TiO2", TiO2},
        {"CaO", CaO},
        {"pH", pH},
        {"Galactic", defaultMaterial}
  };
  
  // read the data through configuration file
  std::ifstream infile ("mainConfig.txt");

  std::string line;

  while (std::getline(infile, line))
  {
    vector<string> row_values;

    split(line, ',', row_values);
    
    if(row_values[0] == "r1") r1 = stod(row_values[1])*cm;    
    else if(row_values[0] == "r2") r2 = stod(row_values[1])*cm; 
    else if(row_values[0] == "r3") r3 = stod(row_values[1])*cm; 
    else if(row_values[0] == "r4") r4 = stod(row_values[1])*cm; 
    else if(row_values[0] == "r5") r5 = stod(row_values[1])*cm; 
  }

  G4cout << "r1 : " << r1 << G4endl;
  G4cout << "r2 : " << r2 << G4endl; 
  G4cout << "r3 : " << r3 << G4endl; 
  G4cout << "r4 : " << r4 << G4endl; 
  G4cout << "r5 : " << r5 << G4endl; 

  
  // fill soil layers
  DefineSoilLayerMaps();
  
  FillSoilLayersWithMaps();
  
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

    // define elements
  elO = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  elH = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  elC  = nistManager->FindOrBuildElement("C");
  elN  = nistManager->FindOrBuildElement("N");
  elMn = nistManager->FindOrBuildElement("Mn");
  elSi = new G4Element("Silicon", "Si", z=14., a= 28.0855*g/mole);

 
  // define nist compounds
//   H2O = new G4Material("Water", density= 1.0*g/cm3, ncomponents=2);
//   H2O->AddElement(elH, natoms=2);
//   H2O->AddElement(elO, natoms=1);
//   // overwrite computed meanExcitationEnergy with ICRU recommended value 
//   H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

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

  
  // create air
  Air = nistManager->FindOrBuildMaterial("G4_AIR");
  // Air = new G4Material("Air", density = 0.1*kg/m3, ncomponents = 2,kStateGas,3000.15*kelvin);
  // Air->AddElement(elN, fractionmass=70.0*perCent);
  // Air->AddElement(elO, fractionmass=30.0*perCent);
  // Air = new G4Material("Air",density= 0.00120479*g/cm3, 1273.15*kelvin, ncomponents=2);
  // Air->AddElement(elN, 70.0*perCent);
  // Air->AddElement(elO, 30.0*perCent);
  
  
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
  OrganicMat = new G4Material("Organic", density = 1.33*g/cm3, ncomponents = 4);
  OrganicMat->AddElement(elC, natoms=50);
  OrganicMat->AddElement(elO, natoms=42);
  OrganicMat->AddElement(elH, natoms=5);
  OrganicMat->AddElement(elN, natoms=3); 
  
  // Vacuum
  defaultMaterial = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
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


  // initialize shell materials
//   shellMaterial_1 = soilOne40W;
//   shellMaterial_2 = soilOne30W;
//   shellMaterial_3 = soilOne20W;
//   shellMaterial_4 = soilOne20W;
//   shellMaterial_5 = soilOne10W;
  shellMaterial_1 = Air;

  // Geometry parameters  
  
  // World
  //
  auto worldSizeXY = 2.5 * r1;
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
  box = new G4Box("box",                          //its name
                   boxSize,boxSize,boxSize);//its size
                                           
  boxLV = new G4LogicalVolume(box,             //its solid
                                   defaultMaterial,                    //its material
                                   "boxLV");               //its name                                
  boxPV = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(),        //at (0,0,0)
                                 boxLV,             //its logical volume
                                 "boxPV",                //its name
                                 worldLV,                      //its mother  volume
                                 false,                  //no boolean operation
                                 0);                     //copy number

  
  //                               
  // Make neutron ball detector (i.e., multiple layers)
  //

  // Define the angles for spheres
  G4double startAnglePhi = 0.0*deg;
  G4double spanningAnglePhi = 360.0*deg;
  G4double startAngleTheta = 0.0*deg;
  G4double spanningAngleTheta = 180.0*deg; //90

  // grain
  auto solidShell_1 = new G4Sphere("grain", 0, grainSize,
			       startAnglePhi, spanningAnglePhi, 
			       startAngleTheta, spanningAngleTheta);  
  grainLV
    = new G4LogicalVolume(
                 solidShell_1,     // its solid
                 shellMaterial_1,  // its material
                 "grainLV");   // its name
                                   
  grainPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 grainLV,          // its logical volume                         
                 "grainPV",    // its name
                 boxLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  


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
  // G4cout << " layer 1 material : " << shellLV_1->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 2 material : " << shellLV_2->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 3 material : " << shellLV_3->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 4 material : " << shellLV_4->GetMaterial()->GetName() << G4endl; 
  // G4cout << " layer 5 material : " << shellLV_5->GetMaterial()->GetName() << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::updatepH(G4double value)
{
    pHValue = value;
    // write code to update value of
    // fractionMassForH and fractionMassForOH 
    fractionMassForH = ((pHValue*100)/14);
    fractionMassForOH = (100.00 - fractionMassForH);
    
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::setLayerMaterial(G4String matName, int layerNumber) 
{
    int flag = 0;
    for (const auto& itr: StringToMaterialMapper) 
    {
        if(itr.first == matName)
        {
            flag = 1;
            auto materialName = itr.second;
            
            // assign this material to proper layer
            if(layerNumber == 1)
            {
                shellMaterial_1 = materialName;
                if(shellLV_1) { shellLV_1->SetMaterial(shellMaterial_1); }
                G4RunManager::GetRunManager()->PhysicsHasBeenModified();
                G4cout<<"layer 1 is set to material: "<<matName<<G4endl;
      
            }
            
            else if(layerNumber == 2)
            {
                shellMaterial_2 = materialName;
                if(shellLV_2) { shellLV_2->SetMaterial(shellMaterial_2); }
                G4RunManager::GetRunManager()->PhysicsHasBeenModified();
                G4cout<<"layer 2 is set to material: "<<matName<<G4endl;
  
            }
            
            else if(layerNumber == 3)
            {
                shellMaterial_3 = materialName;
                if(shellLV_3) { shellLV_3->SetMaterial(shellMaterial_3); }
                G4RunManager::GetRunManager()->PhysicsHasBeenModified();
                G4cout<<"layer 3 is set to material: "<<matName<<G4endl;
        
            }
            
            else if(layerNumber == 4)
            {
                shellMaterial_4 = materialName;
                if(shellLV_4) { shellLV_4->SetMaterial(shellMaterial_4); }
                G4RunManager::GetRunManager()->PhysicsHasBeenModified();
                G4cout<<"layer 4 is set to material: "<<matName<<G4endl;
            
            }
            
            else if(layerNumber == 5)
            {
                shellMaterial_5 = materialName;
                if(shellLV_5) { shellLV_5->SetMaterial(shellMaterial_5); }
                G4RunManager::GetRunManager()->PhysicsHasBeenModified();
                G4cout<<"layer 5 is set to material: "<<matName<<G4endl;
            
            }
            
            else 
            {
                G4cout << "layerNumber : " << layerNumber << "is invalid" << G4endl;
                exit(0);
            }
            
            
            G4RunManager::GetRunManager()->PhysicsHasBeenModified();
            flag = 1;
            break;
        }
    }
    if(flag == 0)
    {
        G4cout<<matName<<" is not found in your datamap"<<G4endl;
 
    }
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

void nbDetectorConstruction::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::DefineSoilLayerMaps()
{
    // define chemical composition maps
    // use these maps in soil layers
    // you can add or delete any material composition
    // in any map
    // just make sure that summation = 100*perCent
    
    // first layer composition
    chem_composition_1 = {
        {SiO2, 46.3*perCent},
        {pH, 15.0*perCent},             // pH
        {Al2O3, 13.0*perCent},
        {Fe2O3, 2.5*perCent},
        {CaO, 1.6*perCent},
        {MgO, 0.7*perCent},
        {TiO2, 0.6*perCent},
        {OrganicMat, 20.3*perCent}      // organicMat
    };
    
    // second layer composition
    chem_composition_2 = {
        {SiO2, 40.17*perCent},
        {pH, 15.0*perCent},             // pH
        {Al2O3, 11.7*perCent},
        {Fe2O3, 2.25*perCent},
        {CaO, 1.44*perCent},
        {MgO, 0.63*perCent},
        {TiO2, 0.54*perCent},
        {OrganicMat, 18.27*perCent},
        {H2O, 10.0*perCent}
    };
    
    // third layer composition
    chem_composition_3 = {
        {SiO2, 34.04*perCent},
        {pH, 15.0*perCent},             // pH
        {Al2O3, 10.4*perCent},
        {Fe2O3, 2.0*perCent},
        {CaO, 1.28*perCent},
        {MgO, 0.56*perCent},
        {TiO2, 0.48*perCent},
        {OrganicMat, 16.24*perCent},
        {H2O, 20.0*perCent}
    };
    
    // fourth layer composition
    chem_composition_4 = {
        {SiO2, 27.91*perCent},
        {pH, 15.0*perCent},             // pH
        {Al2O3, 9.1*perCent},
        {Fe2O3, 1.75*perCent},
        {CaO, 1.12*perCent},
        {MgO, 0.49*perCent},
        {TiO2, 0.42*perCent},
        {OrganicMat, 14.21*perCent},    // organicMat
        {H2O, 30.0*perCent}
    };  
    
    // fifth layer composition
    chem_composition_5 = {
        {SiO2, 17.91*perCent},
        {pH, 15.0*perCent},
        {Al2O3, 9.1*perCent},
        {Fe2O3, 1.75*perCent},
        {CaO, 1.12*perCent},
        {MgO, 0.49*perCent},
        {TiO2, 0.42*perCent},
        {OrganicMat, 14.21*perCent},
        {H2O, 40.0*perCent}
    };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void nbDetectorConstruction::FillSoilLayersWithMaps()
{
  // create soil layers
  // Based on: http://gfnun.unal.edu.co/fileadmin/content/gruposdeinvestigacion/fisicanuclear/Tesis/DanielAndrade_TG.pdf
  // layer 1
  soilOne = new G4Material("DrySoil", density = 0.6*g/cm3, ncomponents=chem_composition_1.size());
  for (it = chem_composition_1.begin(); it != chem_composition_1.end(); it++) {
        soilOne->AddMaterial(it->first, fractionmass=it->second);
  }
  
  // layer 2
  // 10% moisture content
  soilOne10W = new G4Material("DrySoil10W", density = 0.6*g/cm3, ncomponents=chem_composition_2.size());
  for (it = chem_composition_2.begin(); it != chem_composition_2.end(); it++) {
        soilOne10W->AddMaterial(it->first, fractionmass=it->second);
  }
  
  
  // layer 3
  // 20% moisture content. Need new density value?
  soilOne20W = new G4Material("DrySoil20W", density = 0.6*g/cm3, ncomponents=chem_composition_3.size());
  for (it = chem_composition_3.begin(); it != chem_composition_3.end(); it++) {
        soilOne20W->AddMaterial(it->first, fractionmass=it->second);
  }
  
  // layer 4
  // 30% moisture content. Need new density value?
  soilOne30W = new G4Material("DrySoil30W", density = 0.6*g/cm3, ncomponents=chem_composition_4.size());
  for (it = chem_composition_4.begin(); it != chem_composition_4.end(); it++) {
        soilOne30W->AddMaterial(it->first, fractionmass=it->second);
  }
    
  // layer 5 (inner most)
  // 40% moisture content. Need new density value?
  soilOne40W = new G4Material("DrySoil40W", density = 0.6*g/cm3, ncomponents=chem_composition_5.size());
  for (it = chem_composition_5.begin(); it != chem_composition_5.end(); it++) {
        soilOne40W->AddMaterial(it->first, fractionmass=it->second);
  }
}
