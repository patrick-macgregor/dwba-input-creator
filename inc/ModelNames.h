// ModelNames.h
// List of optical models and associated particles
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
#ifndef _MODEL_NAMES_H_
#define _MODEL_NAMES_H_

#include <map>

const Int_t gNUM_MODELS = 14;

// Global list of names
enum class eModel {
	// Protons
	BecchettiGreenlees,
	KoningDelaroche,
	Menet,
	Perey,
	Varner,
	// Deuterons
	AnCai,
	Bojowald,
	DaehnickNR,
	DaehnickR,
	HanShiShen,
	LohrHaeberli,
	PereyPerey,
	// Helium 3
	Pang,
	// Alpha
	BassaniPicard
};

enum class eParticle {
	proton,
	deuteron,
	helium3,
	alpha
};


std::map<eModel,eParticle> gOPTICAL_MODEL_LIST = {
	{ eModel::BecchettiGreenlees, eParticle::proton },
	{ eModel::KoningDelaroche, eParticle::proton },
	{ eModel::Menet, eParticle::proton },
	{ eModel::Perey, eParticle::proton },
	{ eModel::Varner, eParticle::proton },
	{ eModel::AnCai, eParticle::deuteron },
	{ eModel::Bojowald, eParticle::deuteron },
	{ eModel::DaehnickNR, eParticle::deuteron },
	{ eModel::DaehnickR, eParticle::deuteron },
	{ eModel::HanShiShen, eParticle::deuteron },
	{ eModel::LohrHaeberli, eParticle::deuteron },
	{ eModel::PereyPerey, eParticle::deuteron },
	{ eModel::Pang, eParticle::helium3 },
	{ eModel::BassaniPicard, eParticle::alpha },
};

#endif
