// Created on 10/13/2021
//
// All of the geometric parameters of the neutron ball detector is define here
//
// Updated on 10/20/2021: hexc, Mayur, Tien, Weisen
// Added code for reading on detector configuration parameter
//
// Updated on 11/10/2021: hexc, Mayur, Tien, Jarvious
// Add more layers in the detector construct with different material properties.

// Updated on 09/24/2022: Mayur
// Made code more modular by adding separate functions for defining composition, filling soil layers

#ifndef nbDetectorConstruction_h
#define nbDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include <map>

using namespace std;

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4Material;
class nbDetectorMessenger;

/// Detector construction class to define materials and geometry.
/// The detector consists of concentric spherical shells. One can have as many shell as we want.
/// In this version, we simply make single shell with inner and outer radius
/// 

class nbDetectorConstruction : public G4VUserDetectorConstruction
{
    public:

      // construction definition
      nbDetectorConstruction();
      // destructor definition
      virtual ~nbDetectorConstruction();
      
    public:

      // virtual function construct() which returns physical world volume
      virtual G4VPhysicalVolume* Construct();
     
    public:
      // methods
      //
      void DefineMaterials();
      void DefineSoilLayerMaps();
      void FillSoilLayersWithMaps();
      void PrintLayersMaterials();
      
      // function for file reading
      void split(const std::string &s, char delim, std::vector<std::string> &elems);
      G4VPhysicalVolume* DefineVolumes();
      
      // data members
      //
      static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
      // magnetic field messenger
      
      G4double r1, r2, r3, r4, r5; // radii for 5 layers
      G4int matType, matType_1, matType_2, matType_3;
      G4Material *shellMaterial_1, *shellMaterial_2, *shellMaterial_3, *shellMaterial_4, *shellMaterial_5;
      G4VPhysicalVolume *shellPV_1, *shellPV_2, *shellPV_3, *shellPV_4, *shellPV_5;    // neutron ball shell physical volume
      G4LogicalVolume *shellLV_1, *shellLV_2, *shellLV_3, *shellLV_4, *shellLV_5; 
      G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    public:
     
      // define common material variables here
      G4Material *H, *H2O, *SiO2, *Al2O3, *Fe2O3, *CaO, *MgO, *TiO2, *OH, *Mn2O3, *Fe, *Mn, *Ra;
      G4Material *pH;
      G4Material *OrganicMat;
      G4Material *defaultMaterial;
      G4Material *Air;
      G4Material *Quartz;
      
      // define common elements
      G4Element *elO, *elH, *elC, *elN, *elSi;
      G4Element *elMn, *elFe, *elRa;
      
      // define soil layer variables of type G4Material here
      G4Material *soilOne, *soilOne10W, *soilOne20W, *soilOne30W, *soilOne40W;
      
      // define object nist manager
      G4NistManager *nistManager;
      
      // some variables to create user defined materials
      G4String symbol;              // periodic table symbol        
      G4double a;                   // mass of a mole;
      G4double z;                   // z=mean number of protons;  
      G4double density;             // density of material
      G4int ncomponents, natoms;    // parameters
      G4double fractionmass;        // parameters
      
      
      // H and OH concentration in pH material is set dynamically
      // via detector messenger
      // declare their variables
      G4double fractionMassForH;
      G4double fractionMassForOH;
      G4double pHValue;
      
      // map iterator
      map<G4Material*, G4double>::iterator it;
      
      // this snippet creates static maps for soil construction layers
      // maps to store framaction masses 
      map<G4String, G4Material*> StringToMaterialMapper;
      map<G4Material*, G4double> chem_composition_1;
      map<G4Material*, G4double> chem_composition_2;
      map<G4Material*, G4double> chem_composition_3;
      map<G4Material*, G4double> chem_composition_4;
      map<G4Material*, G4double> chem_composition_5;
      
    public:
      // user defined functions to get physical volume names of each layer
      G4String getNameOfLayer1(); // shellPV
      G4String getNameOfLayer2(); // shellPV_1
      G4String getNameOfLayer3(); // shellPV_2
      G4String getNameOfLayer4(); // shellPV_3
      G4String getNameOfLayer5(); // shellPV_4
      G4String getWorld(); // world
      
    
    public:
      // create pointer to detector messenger
      nbDetectorMessenger* fdetectorMessenger;
      
      // function definitions for detector messenger
      void updatepH(G4double pHValue);
      void setLayerMaterial(G4String matName, int layerNumber);
      void setLayerHeight(G4double height, int layerNumber);        
};

#endif
