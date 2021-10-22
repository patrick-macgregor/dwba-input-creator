// TNuclide.h
// Holds information about individual nuclei
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include <TString.h>

class TNuclide{
	private:
		Int_t fZ;
		Int_t fA;
		TString fElement;
		Double_t fMassExcess;
		Double_t fMassExcessErr;
		
	public:
		// Constructors and Destructors
		TNuclide()
		TNuclide( TString name ){ this->SetNuclide(name); }
		TNuclide( Int_t A, Int_t Z ){ this->SetNuclide(A,Z); }
		TNuclide( const TNuclide &a );
		~TNuclide();
	
		// Getters
		Int_t GetZ( return fZ; }
		Int_t GetA( return fA; }
		Int_t GetN( return fA - fZ; }
		TString GetElement( return fElement; }
		Double_t GetMassExcess( return fMassExcess; }
		Double_t GetMassExcessErr( return fMassExcessErr; }
		
		// Non-trivial getters
		Double_t GetMass;
		
		// Setters
		void SetNuclide( TString name );
		void SetNuclide( Int_t A, Int_t Z );
};
