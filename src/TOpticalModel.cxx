// TOpticalModelDeuteron.cxx
// Class functions that calculate optical model parameters for deuterons
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include "../inc/TOpticalModel.h"
#include <TMath.h>

TOpticalModel::TOpticalModel(){
	Initialise();
	fReaction = NULL;
}
// --------------------------------------------------------------------------------------------- //
TOpticalModel::TOpticalModel( TTransferReaction *reaction, eModel model, Bool_t when ){
	fReaction = reaction;
	fWhenReactionOccurs = when;
	this->SetModel(model);
}
// --------------------------------------------------------------------------------------------- //
TOpticalModel::~TOpticalModel(){};
// --------------------------------------------------------------------------------------------- //
void TOpticalModel::Initialise(){
	fv = 0;
	fvi = 0;
	fvsi = 0;
	fvso = 0;
	fvsoi = 0;
	fr0 = 0;
	fri0 = 0;
	frsi0 = 0;
	frso0 = 0;
	frsoi0 = 0;
	fa = 0;
	fai = 0;
	fasi = 0;
	faso = 0;
	fasoi = 0;
};
// --------------------------------------------------------------------------------------------- //
void TOpticalModel::SetModel( eModel model ){
	// Check model matches
	if ( gOPTICAL_MODEL_LIST[ model ] == fReaction->GetTransferParticle( fWhenReactionOccurs ) ){
		fModel = model;
	}
	else{
		std::cout << "Chosen model does not match particle in transfer reaction." << std::endl;
	}
	return;
}
// --------------------------------------------------------------------------------------------- //
Int_t TOpticalModel::CalculateParameters(){
	// Calculate A^1/3
	Int_t A = fReaction->GetBigA( fWhenReactionOccurs );
	Int_t Z = fReaction->GetBigZ( fWhenReactionOccurs );
	Int_t N = A - Z;
	Double_t E = fReaction->GetSystemEnergy( fWhenReactionOccurs );
	Double_t Ebeam = fReaction->GetEbeam();
	Double_t Aonethird = TMath::Power(A, 1.0/3.0);
	eParticle particle = fReaction->GetTransferParticle( fWhenReactionOccurs );
	
	// PROTON MODELS --------------------------------------------------------------------------- //
	if ( particle == eParticle::proton ){
	
		// BecchettiGreenlees
		if ( fModel == eModel::BecchettiGreenlees ){
			fv = 54.0 - 0.32*E + 0.4*Z/Aonethird + 24.0*( N - Z )/A;
			fvi = ( 0.22*E - 2.7 > 0 ? 0.22*E - 2.7 : 0 );
			fvsi = ( 11.8 - 0.25*E + 12.0*( N - Z )/A > 0 ? 11.8 - 0.25*E + 12.0*( N - Z )/A : 0 );
			fvso = 6.2;
			fvsoi = 0.0;
			
			fr0 = 1.17;
			fri0 = 1.32;
			frsi0 = 1.32;
			frso0 = 1.01;
			frsoi0 = 0.0;
			
			fa = 0.75;
			fai = 0.51 + 0.7*( N - Z )/A;
			fasi = 0.51 + 0.7*( N - Z )/A;
			faso = 0.75;
			fasoi = 0.0;
			
			frc0 = 1.3;
		}
		
		// KoningDelaroche
		else if ( fModel == eModel::KoningDelaroche ){
			// Temp parameters
			Double_t vp1 = 59.3 + 21*(N - Z)/A - 0.024*A;
			Double_t vp2 = 0.007067 + (4.23e-6)*A;
			Double_t vp3 =  (1.729e-5) + (1.136e-8)*A;
			Double_t vp4 = 7e-9;
			
			Double_t wp1 = 14.667 + 0.009629*A;
			Double_t wp2 = 73.55 + 0.0795*A;
			
			Double_t dp1 = 16*(1 + (N - Z)/A );
			Double_t dp2 = 0.018 + 0.003802/( 1 + TMath::Exp( (A - 156)/8 ) );
			Double_t dp3 = 11.5;
			
			Double_t vpso1 = 5.922 + 0.0030*A;
			Double_t vpso2 = 0.0040;
			
			Double_t wpso1 = -3.1;
			Double_t wpso2 = 160;
			
			Double_t epf = -8.4075 + 0.01378*A;
			Double_t rc = 1.198 + 0.697/TMath::Power(Aonethird,2) + 12.994/TMath::Power(Aonethird,5);
			
			Double_t vc = 1.73*Z/(rc*Aonethird);
			
			// Now calculate desired parameters
			fv = vp1*( 1 - vp2*(E - epf) + vp3*TMath::Power(E-epf,2) - vp4*TMath::Power(E-epf, 3) ) + vc*vp1*( vp2 - 2*vp3*(E-epf) + 3*vp4*TMath::Power(E-epf,2) );
			fvi = wp1*TMath::Power(E-epf,2)/( TMath::Power(E-epf,2) + TMath::Power(wp2,2) );
			fvsi = dp1*( TMath::Power(E-epf,2) )/(  TMath::Power(E-epf,2) + TMath::Power(dp3,2)  )*TMath::Exp( -dp2*(E-epf) );
			fvso = vpso1*TMath::Exp( -vpso2*(E-epf) );
			fvsoi = wpso1*TMath::Power(E-epf,2 )/( TMath::Power(E-epf,2) + TMath::Power(wpso2,2) );
			
			fr0 = 1.3039 - 0.4054/Aonethird;
			fri0 = fr0;
			frsi0 = 1.3424 - 0.01585*Aonethird;
			frso0 = 1.1854 - 0.647/Aonethird;
			frsoi0 = frso0;
			
			fa = 0.6778 - 0.0001487*A;
			fai = fa;
			fasi = 0.5187 + 0.0005205*A;
			faso = 0.59;
			fasoi = faso;
			
			frc0 = rc;
			return 0;
		}
		
		// Menet
		else if ( fModel == eModel::Menet ){
			fv = 49.9 - 0.22*E + 26.4*( N - Z )/A + 0.4*Z/Aonethird;
			fvi = 1.2 + 0.09*E;
			fvsi = 4.2 - 0.05*E + 15.5*( N - Z )/A;
			fvso = 6.04;
			fvsoi = 0.0;
			
			fr0 = 1.16;
			fri0 = 1.37;
			frsi0 = 1.37;
			frso0 = 1.064;
			frsoi0 = 0.0;
			
			fa = 0.75;
			fai = 0.74 - 0.008*E + ( N - Z )/A;
			fasi = 0.74 - 0.008*E + ( N - Z )/A;
			faso = 0.78;
			fasoi = 0;
			
			frc0 = 1.25;
		}
		
		// Perey
		else if ( fModel == eModel::Perey ){
			fv = 53.3 - 0.55*Ebeam + 27.0*(N - Z)/A + 0.4*Z/Aonethird;
			fvi = 0.0;
			fvsi = 13.5;
			fvso = 7.5;
			fvsoi = 0.0;
			
			fr0 = 1.25;
			fri0 = 0.0;
			frsi0 = 1.25;
			frso0 = 1.25;
			frsoi0 = 0.0;
			
			fa = 0.65;
			fai = 0.0;
			fasi = 0.47;
			faso = 0.47;
			fasoi = 0.0;
			
			frc0 = 1.25;
		}
		
		// Varner
		else if ( fModel == eModel::Varner ){
			Double_t rc1 = 1.24*Aonethird + 0.12;
			Double_t ec = 1.73*Z/rc1;
			Double_t eta =  N - Z/A;

			fv = 52.9 + (13.1*(N - Z)/A ) + ( -0.299*( E - ec) );
			fvi = 7.8/( 1 + TMath::Exp( ( 35 - ( E - ec ) )/16.0 ) );
			fvsi = ( 10 + ( 18.0*(N - Z)/A ) )/( 1 + TMath::Exp( ( E - ec - 36.0 )/37.0 ) );
			fvso = 5.9;
			fvsoi = 0.0;
			
			fr0 = ( ( 1.25*Aonethird ) - 0.225 )/Aonethird;
			fri0 = ( ( 1.33*Aonethird ) - 0.42 )/Aonethird;
			frsi0 = ( ( 1.33*Aonethird ) - 0.42 )/Aonethird;
			frso0 = ( ( 1.34*Aonethird ) - 1.2 )/Aonethird;
			frsoi0 = 0.0;
			
			fa = 0.69;
			fai = 0.69;
			fasi = 0.69;
			faso = 0.63;
			fasoi = 0.0;

			frc0 = rc1/Aonethird;
		}
		
		else{ return 1; }
	}
	
	
	// DEUTERON MODELS ------------------------------------------------------------------------- //
	else if ( particle == eParticle::deuteron ){
		// AnCai
		if ( fModel == eModel::AnCai ){
			fv = 91.85 - 0.249*E + (1.16e-4)*TMath::Power(E,2) + 0.642*Z/Aonethird;
			fvi = 1.104 + 0.0622*E;
			fvsi = 10.83 - 0.0306*E;
			fvso = 3.557;
			fvsoi = 0;
			
			fr0 = 1.152 - 0.00776/Aonethird;
			fri0 = 1.305 + 0.0997/Aonethird;
			frsi0 = 1.334 + 0.152/Aonethird;
			frso0 = 0.972;
			frsoi0 = 0;
			
			fa = 0.719 + 0.0126*Aonethird;
			fai = 0.855 - 0.1*Aonethird;
			fasi = 0.531 + 0.062*Aonethird;
			faso = 1.011;
			fasoi = 0;
			
			frc0 = 1.303;
			return 0;
		}
		
		// Bojowald
		else if ( fModel == eModel::Bojowald ){
			fv = 81.33 + 1.43*Z/Aonethird - 0.24*E;
			fvi = ( 0.132*( E - 45.0 ) > 0 ? 0.132*( E - 45.0 ) : 0 );
			fvsi = 7.8 + 1.04*Aonethird - 0.712*fvi;
			fvso = 6.0;
			fvsoi = 0.0;
			
			fr0 = 1.18;
			fri0 = 1.27;
			frsi0 = 1.27;
			frso0 = 0.78 + 0.038*Aonethird;
			frsoi0 = 0.0;
			
			fa = 0.636 + 0.035*Aonethird;
			fai = 0.768 + 0.021*Aonethird;
			fasi = 0.768 + 0.021*Aonethird;
			faso = 0.78 + 0.038*Aonethird;
			fasoi = 0.0;
			
			frc0 = 1.3;
		}
		
		// DaehnickNR
		else if ( fModel == eModel::DaehnickNR ){
			Double_t mu[6] = {8, 20, 28, 50, 82, 126};
			Double_t step1[6];
			Double_t step2[6];
			Double_t step3 = 0;
			for ( Int_t i = 0; i < 6; i++ ){
				step1[i] = TMath::Power( 0.5*( mu[i] - N ), 2 );
				step2[i] = TMath::Exp( -step1[i] );
				step3 += step2[i];
			}
			Double_t beta = -1.0*TMath::Power( 0.01*E, 2 );
		
			fv = 88.5 - 0.26*E + 0.88*Z/Aonethird;
			fvi = ( 12.2 + 0.026*E )*( 1.0 - TMath::Exp(beta) );
			fvsi = ( 12.2 + 0.026*E )*TMath::Exp(beta);
			fvso = 7.33 - 0.029*E;
			fvsoi = 0.0;
			
			fr0 = 1.17;
			fri0 = 1.325;
			frsi0 = 1.325;
			frso0 = 1.07;
			frsoi0 = 0.0;
			
			fa = 0.709 + 0.0017*E;
			fai = 0.53 + 0.07*Aonethird - 0.04*step3;
			fasi = 0.53 + 0.07*Aonethird - 0.04*step3;
			faso = 0.66;
			fasoi = 0.0;
			
			frc0 = 1.3;
		}
		
		// DaehnickR
		else if ( fModel == eModel::DaehnickR ){
			Double_t mu[6] = {8, 20, 28, 50, 82, 126};
			Double_t step1[6];
			Double_t step2[6];
			Double_t step3 = 0;
			for ( Int_t i = 0; i < 6; i++ ){
				step1[i] = TMath::Power( 0.5*( mu[i] - N ), 2 );
				step2[i] = TMath::Exp( -step1[i] );
				step3 += step2[i];
			}
			Double_t beta = -1.0*TMath::Power( 0.01*E, 2 );
			
			fv =  88.0 - 0.283*E + 0.88*Z/Aonethird;
			fvi = ( 12 + 0.031*E )*( 1 - TMath::Exp(beta) );
			fvsi = ( 12 + 0.031*E )*TMath::Exp(beta);
			fvso = 7.2 - 0.032*E;
			fvsoi = 0.0;
			
			fr0 = 1.17;
			fri0 = 1.376 - 0.01*TMath::Sqrt(E);
			frsi0 = 1.376 - 0.01*TMath::Sqrt(E);
			frso0 = 1.07;
			frsoi0 = 0.0;
			
			fa = 0.717 + 0.0012*E;
			fai = 0.52 + 0.07*Aonethird - 0.04*step3;
			fasi = 0.52 + 0.07*Aonethird - 0.04*step3;
			faso = 0.66;
			fasoi = 0.0;
			
			frc0 = 1.3;
		}
		
		// HanShiShen
		else if ( fModel == eModel::HanShiShen ){
			fv = 82.18 - 0.148*E - 0.000886*E*E - 34.811*( N - Z )/A + 1.058*Z/Aonethird;
			fvi = ( -4.916 + 0.0555*E + 0.0000442*E*E + 35.0*( N - Z )/A > 0 ? -4.916 + 0.0555*E + 0.0000442*E*E + 35.0*( N - Z )/A : 0 );
			fvsi = 20.968 - 0.0794*E - 43.398*( N - Z )/A;
			fvso = 3.703;
			fvsoi = -0.206;
			
			fr0 = 1.174;
			fri0 = 1.563;
			frsi0 = 1.328;
			frso0 = 1.234;
			frsoi0 = 1.234;
			
			fa = 0.809;
			fai = 0.7 + 0.045*Aonethird;
			fasi = 0.465 + 0.045*Aonethird;
			faso = 0.813;
			fasoi = 0.813;
			
			frc0 = 1.698;
		}
		
		// LohrHaeberli
		else if ( fModel == eModel::LohrHaeberli ){
			fv = 91.13 + 2.2*Z/Aonethird;
			fvi = 0.0;
			fvsi = 218/(Aonethird*Aonethird);
			fvso = 7.0;
			fvsoi = 0.0;
			
			fr0 = 1.05;
			fri0 = 0.0;
			frsi0 = 1.43;
			frso0 = 0.75;
			frsoi0 = 0.0;
			
			fa = 0.86;
			fai = 0.0;
			fasi = 0.5 + 0.013*Aonethird*Aonethird;
			faso = 0.5;
			fasoi = 0.0;
			
			frc0 = 1.30;
		}
		
		// PereyPerey
		else if ( fModel == eModel::PereyPerey ){
			fv = 81 - 0.22*E + 2*Z/Aonethird;
			fvi = 0.0;
			fvsi = 14.4 + 0.24*E;
			fvso = 0.0;
			fvsoi = 0.0;
			
			fr0 = 1.15;
			fri0 = 0.0;
			frsi0 = 1.34;
			frso0 = 0.0;
			frsoi0 = 0.0;
			
			fa = 0.81;
			fai = 0.0;
			fasi = 0.68;
			faso = 0.0;
			fasoi = 0.0;
			
			frc0 = 1.15;
		}
		
		else{ return 1; }
	}
	
	
	// HELIUM 3 MODELS ------------------------------------------------------------------------- //
	else if ( particle == eParticle::helium3 ){
		
		// Pang
		if ( fModel == eModel::Pang ){
			Double_t rc = 1.24*Aonethird + 0.12;
			Double_t EC = 1.728*Z*2/rc;
			Double_t ETA = (N-Z)/A;
			Double_t VSI_ASYM = 35 + (34.2*ETA);
			
			fv = 118.3 + (-0.13*( Ebeam - EC ) );
			fvi = 38.5/( 1 + TMath::Exp( ( 156.1 - ( Ebeam - EC ) )/52.4 ) );
			fvsi = VSI_ASYM/( 1 + TMath::Exp( ( ( Ebeam - EC ) - 30.8 )/106.4 ) );
			if ( Ebeam < 85 ){ fvso = 1.7 + (-0.02*Ebeam); }
			else{ fvso = 0;	}
			fvsoi = 0;
			
			fr0 = ( 1.3*Aonethird - 0.48 )/Aonethird;
			fri0 = ( 1.31*Aonethird - 0.13 )/Aonethird;
			frsi0 = fri0;
			frso0 = ( 0.64*Aonethird +1.18 )/Aonethird;
			frsoi0 = 0.0;
			
			fa = 0.820;
			fai = 0.840;
			fasi = 0.840;
			faso = 0.130;
			fasoi = 0.0;
			
			frc0 = rc/Aonethird;
		}
		
		else{ return 1; }
	}
	
	
	// ALPHA MODELS ---------------------------------------------------------------------------- //
	else if( particle == eParticle::alpha ){
	
		// BassaniPicard
		if ( fModel == eModel::BassaniPicard ){
			fv = 207.0;
			fvi = 28.0;
			fvsi = 0.0;
			fvso = 0.0;
			fvsoi = 0.0;
			
			fr0 = 1.30;
			fri0 = 1.30;
			frsi0 = 0.0;
			frso0 = 0.0;
			frsoi0 = 0.0;
			
			fa = 0.65;
			fai = 0.52;
			fasi = 0.0;
			faso = 0.0;
			fasoi = 0.0;
			
			frc0 = 1.40;
		}
		else{ return 1; }
	}
	
	
	// No matched particles for some reason...
	else{
		return 1;
	}
	
	
	
	return 0;
}































