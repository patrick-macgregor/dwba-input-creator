// TTransferReaction.cxx
// Defines transfer reaction parameters
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include "../inc/TTransferReaction.h"
// --------------------------------------------------------------------------------------------- //
TTransferReaction::TTransferReaction( Int_t A1, Int_t Z1, Int_t A2, Int_t Z2, Int_t A3, Int_t Z3, Int_t A4, Int_t Z4 ){
	Initialise( A1, Z1, A2, Z2, A3, Z3, A4, Z4 );
}
// --------------------------------------------------------------------------------------------- //
TTransferReaction::TTransferReaction( TString name1, TString name2, TString name3, TString name4 ){
	Int_t A1 = 0, Z1 = 0, A2 = 0, Z2 = 0, A3 = 0, Z3 = 0, A4 = 0, Z4 = 0;
	if ( ParseNuclide( name1, A1, Z1 ) == 0 && ParseNuclide( name2, A2, Z2 ) == 0 && ParseNuclide( name3, A3, Z3 ) == 0 && ParseNuclide( name4, A4, Z4 ) == 0 ){
		Initialise( A1, Z1, A2, Z2, A3, Z3, A4, Z4 );
	}
}
// --------------------------------------------------------------------------------------------- //
TTransferReaction::~TTransferReaction(){}
// --------------------------------------------------------------------------------------------- //
void TTransferReaction::Initialise( Int_t A1, Int_t Z1, Int_t A2, Int_t Z2, Int_t A3, Int_t Z3, Int_t A4, Int_t Z4 ){
	this->SetNuclide( 1, A1, Z1 );
	this->SetNuclide( 2, A2, Z2 );
	this->SetNuclide( 3, A3, Z3 );
	this->SetNuclide( 4, A4, Z4 );
	CalculateQ();
	
	// Set these to default values
	fEx = 0;
	fEbeam = 0;
}
// --------------------------------------------------------------------------------------------- //
void TTransferReaction::SetNuclide( Int_t num, TString name ){
	Int_t A = 0, Z = 0;
	ParseNuclide(name, A, Z);
	this->SetNuclide( num, A, Z );
	return;
}
// --------------------------------------------------------------------------------------------- //
void TTransferReaction::SetNuclide( Int_t num, Int_t A, Int_t Z ){
	if ( num > 0 && num < 5 ){
		// Set the nuclide
		fNuclide[num-1]->SetNuclide( A, Z );
		
		if ( num == 2 || num == 3 ){
			if ( A == 1 && Z == 1 ){ fTransferParticle[num - 2] = eParticle::proton; }
			else if ( A == 2 && Z == 1 ){ fTransferParticle[num - 2] = eParticle::proton; }
			else if ( A == 3 && Z == 2 ){ fTransferParticle[num - 2] = eParticle::proton; }
			else if ( A == 4 && Z == 2 ){ fTransferParticle[num - 2] = eParticle::proton; }
			else{
				std::cout << "(Currently) don't have optical models for " << A << ELEMENT_SYMBOL[A] << std::cout
			}
		}
	}
	else{
		std::cout << num << " is not a valid nucleus number. Must be 1, 2, 3, or 4." << std::endl;
	}
	return;
}
// --------------------------------------------------------------------------------------------- //
void TTransferReaction::SetEx( Double_t ex ){ fEx = ex; }
// --------------------------------------------------------------------------------------------- //
void TTransferReaction::CalculateQ(){
	fQ = fNuclide[0]->GetMassExcess() + fNuclide[1]->GetMassExcess() - fNuclide[2]->GetMassExcess() - fNuclide[3]->GetMassExcess();
}
// --------------------------------------------------------------------------------------------- //
void TTransferReaction::SetEbeam( Double_t ebeam ){ fEbeam = ebeam; }
// --------------------------------------------------------------------------------------------- //
eParticle TTransferReaction::GetTransferParticle( Bool_t when_reaction_occurs ){
	return fTransferParticle[ when_reaction_occurs ];
}
// --------------------------------------------------------------------------------------------- //
Int_t GetBigA( Bool_t when_reaction_occurs ){
	if ( when_reaction_occurs == 0 ){ return fNuclide[0]->GetA(); }
	else{ return fNuclide[3]->GetA(); }
}
// --------------------------------------------------------------------------------------------- //
Int_t GetBigZ( Bool_t when_reaction_occurs ){
	if ( when_reaction_occurs == 0 ){ return fNuclide[0]->GetZ(); }
	else{ return fNuclide[3]->GetZ(); }
}
// --------------------------------------------------------------------------------------------- //
Double_t GetSystemEnergy( Bool_t when_reaction_occurs ){
	if ( when_reaction_occurs == 0 ){ return fEbeam; }
	else{ return fEbeam + fQ - fEx; }
}
// --------------------------------------------------------------------------------------------- //







