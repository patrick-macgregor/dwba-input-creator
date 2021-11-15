// TOpticalModelDeuteron.h
// Class that calculates optical model parameters for deuterons
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#ifndef _T_OPTICAL_MODEL_H_
#define _T_OPTICAL_MODEL_H_

#include "TTransferReaction.h"
#include "ModelNames.h"

class TOpticalModel{
	private:
		// Actual optical model parameters
		Double_t fv;
		Double_t fvi;
		Double_t fvsi;
		Double_t fvso;
		Double_t fvsoi;

		Double_t fr0;
		Double_t fri0;
		Double_t frsi0;
		Double_t frso0;
		Double_t frsoi0;
		
		Double_t fa;
		Double_t fai;
		Double_t fasi;
		Double_t faso;
		Double_t fasoi;
		
		Double_t frc0;
		
		// Additional info
		eModel fModel;
		TTransferReaction *fReaction;
		Bool_t fWhenReactionOccurs; // 0 = before, 1 = after
	
	public:
		TOpticalModel();
		TOpticalModel( TTransferReaction *reaction, eModel model, Bool_t when );
		~TOpticalModel();
		
		void Initialise();
		void SetModel( eModel model );
		Int_t CalculateParameters();
};

#endif
