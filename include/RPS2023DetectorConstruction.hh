#ifndef RPS2023DetectorConstruction_h
#define RPS2023DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include <map>

using namespace std;

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4Material;
class RPS2023DetectorMessenger;

/// Detector construction class to define materials and geometry.
/// The detector consists of concentric spherical shells. One can have as many shell as we want.
/// In this version, we simply make single shell with inner and outer radius
/// 

class RPS2023DetectorConstruction : public G4VUserDetectorConstruction
{
    public:

      // construction definition
      RPS2023DetectorConstruction();
      // destructor definition
      virtual ~RPS2023DetectorConstruction();
      
    public:

      // virtual function construct() which returns physical world volume
      virtual G4VPhysicalVolume* Construct();
     
    public:
      // methods
      //
      void DefineMaterials();
      void DefineApproximateChemicalComps();  // use this when you need to use approximate chemical composition
      void DefineXRFChemicalComps();    // use this when you need to use XRF chemical composition 
      void FillSoilLayersWithMaps();
      void PrintLayersMaterials();
      void fillGrainWithChemComps();

      void DefinePoreChemicalComps();
      void fillPoresWithChemComps();

      G4double calculateMixtureDensity(map<G4Material*, G4double> chem_composition_map);

      G4VPhysicalVolume* DefineVolumes();
      
      G4Material* grainMaterial;
      G4Material* poreMaterial;

      G4double r1, r2, r3, r4, r5; // radii for 5 layers
      G4int matType, matType_1, matType_2, matType_3;
      G4Material *shellMaterial_1, *shellMaterial_2, *shellMaterial_3, *shellMaterial_4, *shellMaterial_5, *shellMaterial_6, *shellMaterial_7, *shellMaterial_8;
      G4VPhysicalVolume *shellPV_1, *shellPV_2, *shellPV_3, *shellPV_4, *shellPV_5, *grainPV, *boxPV, *surroundingGrainPV, *poreGrainPV;    // neutron ball shell physical volume
      G4LogicalVolume *shellLV_1, *shellLV_2, *shellLV_3, *shellLV_4, *shellLV_5, *grainLV, *boxLV, *surroundingGrainLV, *poreGrainLV; 
      G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

      // grain size
      G4double grainSize = 1000000.*nm;

      // Define dimensions of box
      G4double boxSizeX = 1*cm;
      G4double boxSizeY = 1*cm;
      G4double boxSizeZ = 1*cm;

      // // Define positions of spheres on surface of box
      // G4double spherePosX = boxSizeX/2 - grainSize - 1*cm;
      // G4double spherePosY = boxSizeY/2 - grainSize - 1*cm;
      // G4double spherePosZ = grainSize;

    public:
     
      // define common material variables here
      G4Material *H, *H2O, *SiO2, *Al2O3, *Fe2O3, *CaO, *MgO, *TiO2, *OH, *Mn2O3;
      G4Material *pH;
      G4Material *OrganicMat;
      G4Material *defaultMaterial;
      G4Material *Air;
      G4Material *Quartz;
      G4Material *Ca, *C, *N, *Na, *Mg, *P, *Si, *K, *Al, *I, *S, *Ti, *O, *V, *Cr, *Fe, *Mn, *Ra, *Co, *Ni, *Cu, *Zn, *As, *Se, *Rb, *Sr, *Y, *Zr, *Nb, *Mo, *Cd, *Sn, *Ba, *La, *Ce, *Nd, *Hg, *Pb, *Bi, *Th, *U;
      G4Material *H2OVapor;
      G4Material *coatingMaterial;
      
      // define common elements
      G4Element *elO, *elH, *elC, *elN, *elSi, *elP, *elMg, *elS, *elK;
      G4Element *elMn, *elFe, *elRa;
      
      // define soil layer variables of type G4Material here
      G4Material *soilOne, *soilOne10W, *soilOne20W, *soilOne30W, *soilOne40W;

      G4Material *grainComp0, *grainComp1, *grainComp2, *grainComp3, *grainComp4, *grainComp5, *grainComp6, *grainComp7, 
                 *grainComp8, *grainComp9, *grainComp10, *grainComp11, *grainComp12, *grainComp13, *grainComp14, *grainComp15,
                 *grainComp16, *grainComp17, *grainComp18, *grainComp19, *grainComp20;

      
      G4Material *poreComp0, *poreComp1, *poreComp2, *poreComp3, *poreComp4, *poreComp5, *poreComp6, *poreComp7, 
                 *poreComp8, *poreComp9, *poreComp10, *poreComp11, *poreComp12, *poreComp13, *poreComp14, *poreComp15,
                 *poreComp16, *poreComp17, *poreComp18, *poreComp19, *poreComp20;
      
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
      map<G4Material*, G4double> chem_composition_0;
      map<G4Material*, G4double> chem_composition_1;
      map<G4Material*, G4double> chem_composition_2;
      map<G4Material*, G4double> chem_composition_3;
      map<G4Material*, G4double> chem_composition_4;
      map<G4Material*, G4double> chem_composition_5;
      map<G4Material*, G4double> chem_composition_6;
      map<G4Material*, G4double> chem_composition_7;
      map<G4Material*, G4double> chem_composition_8;
      map<G4Material*, G4double> chem_composition_9;
      map<G4Material*, G4double> chem_composition_10;
      map<G4Material*, G4double> chem_composition_11;
      map<G4Material*, G4double> chem_composition_12;
      map<G4Material*, G4double> chem_composition_13;
      map<G4Material*, G4double> chem_composition_14;
      map<G4Material*, G4double> chem_composition_15;
      map<G4Material*, G4double> chem_composition_16;
      map<G4Material*, G4double> chem_composition_17;
      map<G4Material*, G4double> chem_composition_18;
      map<G4Material*, G4double> chem_composition_19;
      map<G4Material*, G4double> chem_composition_20;

      map<G4Material*, G4double> pore_chem_composition_0;
      map<G4Material*, G4double> pore_chem_composition_1;
      map<G4Material*, G4double> pore_chem_composition_2;
      map<G4Material*, G4double> pore_chem_composition_3;
      map<G4Material*, G4double> pore_chem_composition_4;
      map<G4Material*, G4double> pore_chem_composition_5;
      map<G4Material*, G4double> pore_chem_composition_6;
      map<G4Material*, G4double> pore_chem_composition_7;
      map<G4Material*, G4double> pore_chem_composition_8;
      map<G4Material*, G4double> pore_chem_composition_9;
      map<G4Material*, G4double> pore_chem_composition_10;
      map<G4Material*, G4double> pore_chem_composition_11;
      map<G4Material*, G4double> pore_chem_composition_12;
      map<G4Material*, G4double> pore_chem_composition_13;
      map<G4Material*, G4double> pore_chem_composition_14;
      map<G4Material*, G4double> pore_chem_composition_15;
      map<G4Material*, G4double> pore_chem_composition_16;
      map<G4Material*, G4double> pore_chem_composition_17;
      map<G4Material*, G4double> pore_chem_composition_18;
      map<G4Material*, G4double> pore_chem_composition_19;
      map<G4Material*, G4double> pore_chem_composition_20;


    public:
      // user defined functions to get physical volume names of each layer
      G4String getNameOfLayer1(); // shellPV
      G4String getNameOfLayer2(); // shellPV_1
      G4String getNameOfLayer3(); // shellPV_2
      G4String getNameOfLayer4(); // shellPV_3
      G4String getNameOfLayer5(); // shellPV_4
      G4String getWorld(); // world
      
    
    public:

      // ANALYSIS VARIABLES
      G4int H2OContent = 0; // initialized with initial H2O content

      // create pointer to detector messenger
      RPS2023DetectorMessenger* fdetectorMessenger;
      
      // function definitions for detector messenger
      void updateH2OContent(G4int);
      void setGrainMaterial(G4int);     
};

#endif
