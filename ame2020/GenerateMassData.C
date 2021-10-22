// GenerateMassData.C
// Takes findings of AME2020 and wraps it nicely into a ROOT file for future access in dwba-calculator code.
// ============================================================================================= //
// Patrick MacGregor
// Nuclear Physics Research Group
// School of Physics and Astronomy
// The University of Manchester
// ============================================================================================= //
// SOURCE: https://www.anl.gov/sites/www/files/2021-04/mass_1.mas20.txt
// ============================================================================================= //
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <string>

#include "../inc/ElementSymbols.h"

TString RemoveErrantChars( TString a ){
	// Remove * and # from AME data
	if ( a.Contains("*") ){
		// This means the value is unknown - set to -1
		a = "-1";
	}
	else if ( a.Contains("#") ){
		// In place of a decimal point
		a.ReplaceAll("#",".");
	}
	return a;
}



Int_t GenerateMassData(){

	// Declare variables
	TString input_file_name = "mass_1.mas20.txt";
	TString ame_year = "2020";
	TString output_file_name = "ame" + ame_year + ".root";
	Int_t start_line = 37;
	Bool_t error_flag = 0;
	
	// Declare columnar properties of file {start pos, length}
	Int_t N_col[2] = {4,5};
	Int_t Z_col[2] = {9,5};
	Int_t A_col[2] = {14,5};
	Int_t element_col[2] = {20,3};
	Int_t mass_excess_col[2] = {28,14};
	Int_t mass_excess_err_col[2] = {42,12};
	
	
	// Declare containers for the tree
	Int_t Z = 0;
	Int_t A = 0;
	Int_t N = 0;
	Double_t mass_excess[2] = {0,0};	// Value and error
	TString element = "";
	
	// Declare tree and branches
	TTree *t = new TTree("ame" + ame_year, "Atomic Mass Evaluation " + ame_year + " Data" );
	t->Branch( "Z", &Z, "Z/I" );
	t->Branch( "A", &A, "Z/I" );
	t->Branch( "N", &N, "Z/I" );
	t->Branch( "mass_excess", mass_excess, "mass_excess[2]/D" );
	t->Branch( "element", "TString", &element );
	
	
	// Read in the mass data line by line
	std::ifstream input_file;
	input_file.open( input_file_name );
	
	// Declare variables to help with data storage
	TString tline = "";
	std::string line = "";
	Int_t line_ctr = 1;
	TObjArray *line_split = NULL;
	
	// Loop over the file
	if ( input_file.is_open() ){
		while ( std::getline( input_file, line ) ){
			// Ignore lines that aren't desired
			if ( line_ctr >= start_line ){
				
				// Cast line to TString
				tline = line;
				
				// Assign values from the array
				N = RemoveErrantChars( tline( N_col[0], N_col[1] ) ).Atoi();
				Z = RemoveErrantChars( tline( Z_col[0], Z_col[1] ) ).Atoi();
				A = RemoveErrantChars( tline( A_col[0], A_col[1] ) ).Atoi();
				element = ((TString)tline( element_col[0], element_col[1] )).ReplaceAll( " ", "" );
				mass_excess[0] = RemoveErrantChars( tline( mass_excess_col[0], mass_excess_col[1] ) ).Atof();
				mass_excess[1] = RemoveErrantChars( tline( mass_excess_err_col[0], mass_excess_err_col[1] ) ).Atof();
				
				// Check N + Z = A
				if ( N + Z != A ){
					error_flag = 1;
					std::cout << "ERROR: " << N << " + " << Z << " != " << A << std::endl;
				}
				
				// Check element strings match
				if ( ELEMENT_SYMBOL[Z] != element ){
					error_flag = 1;
					std::cout << "ERROR: " << element << "(in file) != " << ELEMENT_SYMBOL[Z] << " (expected)" << std::endl;
				}
			
				// Store the desired data in the tree
				t->Fill();
				
				
			}
			line_ctr++;
		}
		
		// Close the data file
		input_file.close();
		
		// Open the TFile and write the tree
		TFile *f = new TFile(output_file_name, "RECREATE" );
		t->Write();
		f->Close();
	}
	else{
		error_flag = 1;
	}
	
	
	// Print final message
	if (error_flag == 0 ){
		std::cout << "File written successfully" << std::endl;
		return 0;
	}
	else{
		std::cout << "Program terminated with some errors" << std::endl;
		return 1;
	}
}
