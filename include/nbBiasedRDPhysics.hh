/*

    seperate physics list develoepd for 
    use of G4Radioactivation process

*/

#ifndef nbBiasedRDPhysics_h
#define nbBiasedRDPhysics_h 1

#include "G4VPhysicsConstructor.hh"


class G4Radioactivation;

class nbBiasedRDPhysics : public G4VPhysicsConstructor
{
  public: 
    nbBiasedRDPhysics(G4int verbose = 1);
    nbBiasedRDPhysics(const G4String& name);
    virtual ~nbBiasedRDPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

};

#endif

