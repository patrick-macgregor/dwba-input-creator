// TAMEAcccessor.cxx
// Linked to TAMEAcccessor.h
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include "../inc/DWBACalculator.h"
#include "../inc/FParseNuclideName.h"
#include "../inc/TAMEAccessor.h"

// Constructor
TAMEAccessor::TAMEAccessor(){

	// Open TFile
	fMassDatabaseFile = TFile::Open( gMASS_DATABASE_LOCATION, "READ" );
	
	if ( !fMassDatabaseFile->IsOpen() ){
		std::cout << "Error accessing mass database" << std::endl;
		fMassDatabaseFile = NULL;
	}
	
	// Get TTree
	fMassDatabaseTree = (TTree*)fMassDatabaseFile->Get(gMASS_DATABASE_TREE_NAME);
	
}

// Destructor
TAMEAccessor::~TAMEAccessor(){
	// Close TFile
	fMassDatabaseFile->Close();
}



// Get mass excess array from the mass database
Int_t TAMEAccessor::GetMassExcess( TString name, Double_t& mass_excess, Double_t& mass_excess_err ){
	Int_t A = 0;
	Int_t Z = 0;
	if ( ParseNuclide( name, A, Z ) == 1 ){
		A = 0;
		Z = 0;
		std::cout << "Cannot extract mass excess array" << std::endl;
		mass_excess = 0;
		mass_excess_err = 0;
		return 1;
	}
	else{
		return GetMassExcess( A, Z, mass_excess, mass_excess_err );
	}
}
Int_t TAMEAccessor::GetMassExcess( Int_t A, Int_t Z, Double_t& mass_excess, Double_t& mass_excess_err ){
	// Loop over the tree and find the correct entry
	Int_t A_tree = 0;
	Int_t Z_tree = 0;
	Double_t mass_excess_tree[2];
	
	fMassDatabaseTree->SetBranchAddress("A", &A_tree);
	fMassDatabaseTree->SetBranchAddress("Z", &Z_tree);
	fMassDatabaseTree->SetBranchAddress("mass_excess", &mass_excess_tree);
	
	for ( Int_t i = 0; i < fMassDatabaseTree->GetEntries(); i++ ){
		fMassDatabaseTree->GetEntry(i);
		if ( A_tree == A && Z_tree == Z ){
			break;
		}
	}
	
	mass_excess = mass_excess_tree[0];
	mass_excess_err = mass_excess_tree[1];

	return 0;
}

// Define global accessor
TAMEAccessor *gAME_ACCESESOR = new TAMEAccessor();




