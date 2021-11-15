// TNuclide.cxx
// Holds information about individual nuclei
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include "../inc/TAMEAccessor.h"
#include "../inc/TNuclide.h"
#include "../inc/FParseNuclideName.h"

void TNuclide::SetNuclide( TString name ){
	ParseNuclide( name, fA, fZ );
	Update();
}

void TNuclide::SetNuclide( Int_t A, Int_t Z ){
	fA = A;
	fZ = Z;
	Update();
}

void TNuclide::Update(){
	gAME_ACCESESOR->GetMassExcess( fA, fZ, fMassExcess, fMassExcessErr );
}

