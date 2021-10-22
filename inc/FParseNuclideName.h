// FParseNuclideName.h
// Takes an input string of a nuclear name and parses the output so that it can be interpreted
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //

#include "ElementSymbols.h"

// String is of the form e.g. "16O"
// RETURN VALUES
//	0 - no errors
//	1 - major error
//	2 - probable error

// Run a test on the nuclide values of Z and A
Int_t CheckNuclideValues( Int_t A, Int_t Z ){
	
	Int_t ret_value = 0;
	
	// CATASTROPHIC errors
	// Check Z and A make sense relative to each other
	if ( Z > A ){
		std::cout << "ERROR -- incorrect number of protons compared to total mass." << std::endl;
		return 1;
	}
	
	// Check there is at least 1 nucleon
	if ( A < 1 ){
		std::cout << "ERROR -- must have at least 1 nucleon" << std::endl;
		return 1;
	}
	
	// Check that there are non-negative protons
	if ( Z < 0 ){
		std::cout << "ERROR -- cannot have negative number of protons" << std::endl;
		return 1;
	}
	
	
	// POSSIBLE errors
	// Check mass number makes sense (throw warning 2 otherwise)
	if ( A > MAX_MASS ){
		std::cout << "ERROR -- mass almost certainly too large. Is this correct?" << std::endl;
		ret_value = 2;
	}
	
	if ( (A-Z)/Z > 4 ){
		std::cout << "ERROR -- too many neutrons maybe?" << std::endl;
		ret_value = 2;
	}
	
	return ret_value;
}


// Parse the nuclide
Int_t ParseNuclide( TString name, Int_t &A, Int_t &Z ){

	// Booleans for checking any errors
	Bool_t b_first_char_not_digit = 0;
	Bool_t b_second_half_not_element = 0;
	Int_t ret_value = 0;


	// Warn user if the length is too long
	if ( name.Length() > 5 ){
		std::cout << "Nuclide name \"" << name << "\" is probably too long. Did you make an error?" << std::endl;
		ret_value = 2;
	}
	
	// Define split between mass and element name
	Int_t split = -1;

	// Loop over the name and extract the quantities
	for ( Int_t i = 0; i < name.Length(); i++ ){
		// Extract first half - find out where the split is
		if ( split == -1 && !std::isdigit( name(i) ) ){ 
			split = i;
			if ( i == 0 ){ b_first_char_not_digit = 1; }
			break;
		}
	}
	
	// Process first warning
	if ( b_first_char_not_digit ){
		std::cout << "ERROR -- no mass detected in \"" << name << "\"" << std::endl;
		return 1;
	}
	else{
		TString mass_string = name(0,split);
		A = mass_string.Atoi();
	}
		
	// Now test if the second half is an element
	TString test_string = name(split, name.Length() - split );
	test_string.ToLower();
	
	// Check second half is an element name
	for ( Int_t i = 0; i < NUM_ELEMENTS; i++ ){
		TString element_check = ELEMENT_SYMBOL[i];
		element_check.ToLower();
		
		if ( element_check == test_string ){
			Z = i;
			break;
		}
		if ( i == NUM_ELEMENTS - 1 ){
			// No match found - raise error
			b_second_half_not_element = 1;
		}
	}
	
	// Process second warning
	if ( b_second_half_not_element ){
		std::cout << "ERROR -- \""  << name(split, name.Length() - split ) << "\" is not recognised as a valid element." << std::endl;
		return 1;
	}
	return CheckNuclideValues(A,Z);
}

























