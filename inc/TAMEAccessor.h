// TAMEAcccessor.h
// Class that accesses the mass database and returns values
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef _T_AME_ACCESSOR_H_
#define _T_AME_ACCESSOR_H_

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

class TAMEAccessor{
	private:
		TFile *fMassDatabaseFile = NULL;
		TTree *fMassDatabaseTree = NULL;

	public:
		TAMEAccessor();
		~TAMEAccessor();
		
		// Getters
		Int_t GetMassExcess( Int_t A, Int_t Z, Double_t& mass_excess, Double_t& mass_excess_err);
		Int_t GetMassExcess( TString name, Double_t& mass_excess, Double_t& mass_excess_err);
};

#endif
