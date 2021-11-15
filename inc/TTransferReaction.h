// TTransferReaction.h
// Defines transfer reaction parameters
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef _T_TRANSFER_REACTION_H_
#define _T_TRANSFER_REACTION_H_

#include "FParseNuclideName.h"
#include "ModelNames.h"
#include "TOpticalModel.h"
#include "TNuclide.h"


class TTransferReaction{
	private:
		TNuclide *fNuclide[4];	// Beam (light particle), Target, Ejectile (light particle), Residual
		eParticle fTransferParticle[2];	// Initial, Final
		Double_t fQ;
		Double_t fEx;
		Double_t fEbeam;

	public:
		// Not allowed constructor with nothing in it!!
		TTransferReaction( TString name1, TString name2, TString name3, TString name4 );
		TTransferReaction( Int_t A1, Int_t Z1, Int_t A2, Int_t Z2, Int_t A3, Int_t Z3, Int_t A4, Int_t Z4 );
		~TTransferReaction();
		
		void Initialise( Int_t A1, Int_t Z1, Int_t A2, Int_t Z2, Int_t A3, Int_t Z3, Int_t A4, Int_t Z4 );
		void CalculateQ();
		
		void SetNuclide( Int_t num, TString name );
		void SetNuclide( Int_t num, Int_t A, Int_t Z );
		void SetEx( Double_t ex );
		void SetEbeam( Double_t ebeam );
		eParticle GetTransferParticle( Bool_t when_reaction_occurs );
		
		Int_t GetBigA( Bool_t when_reaction_occurs );
		Int_t GetBigZ( Bool_t when_reaction_occurs );
		Double_t GetSystemEnergy( Bool_t when_reaction_occurs );
		Double_t GetEbeam(){ return fEbeam; }
};


#endif
