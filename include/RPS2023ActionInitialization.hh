#ifndef RPS2023ActionInitialization_h
#define RPS2023ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class RPS2023DetectorConstruction;

/// Action initialization class.
///

class RPS2023ActionInitialization : public G4VUserActionInitialization
{
  public:
    
    // constructor & destructor
   RPS2023ActionInitialization(RPS2023DetectorConstruction* detector);
    virtual  ~RPS2023ActionInitialization();
    
    //  inbuilt methods in actionInitialization class
    virtual void BuildForMaster() const;
    virtual void Build() const;
   
  private:
    // create object of detector construction class
   RPS2023DetectorConstruction* fDetector;
};

#endif

    

