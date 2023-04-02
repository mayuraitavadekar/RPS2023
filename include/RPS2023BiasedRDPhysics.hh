/*

    seperate physics list develoepd for 
    use of G4Radioactivation process

*/

#ifndef RPS2023BiasedRDPhysics_h
#define RPS2023BiasedRDPhysics_h 1

#include "G4VPhysicsConstructor.hh"


class G4Radioactivation;

class RPS2023BiasedRDPhysics : public G4VPhysicsConstructor
{
  public: 
    RPS2023BiasedRDPhysics(G4int verbose = 1);
    RPS2023BiasedRDPhysics(const G4String& name);
    virtual  ~RPS2023BiasedRDPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

};

#endif

