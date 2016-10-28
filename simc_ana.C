#include "TFile.h"
#include "TTree.h"
#include "extras.C"

int counter1 = 0;

Float_t Weight [10];     // Weight
Float_t Em     [10];	// Missing Energy
Float_t Pm     [10];   	// Missing Momentum
Float_t PmPar  [10];	// Missing Momentum component parallel to q
Float_t Q2     [10];	// Momentum transfer squared
Float_t Normfac[10];	// Normalization factor

Float_t h_ytar [10];	Float_t h_delta [10];	Float_t h_yptar [10];	Float_t h_xptar [10];
Float_t h_ytari[10];	Float_t h_deltai[10];	Float_t h_yptari[10];	Float_t h_xptari[10];

Float_t e_ytar [10];	Float_t e_delta [10];	Float_t e_yptar [10];	Float_t e_xptar [10];
Float_t e_ytari[10];	Float_t e_deltai[10];	Float_t e_yptari[10];	Float_t e_xptari[10];

Float_t h_pf [10];	// Outgoing proton momentum
Float_t h_Thf[10];      // Outgoing proton in-plane angle

// labels
int electron = 0;	int hadron = 1;		int both  = 2;
int ytar     = 0;	int delta  = 1;		int yptar = 2;
int xptar    = 3;	int moment = 4;		int pmiss = 5;
int emiss    = 6;	int Qsq    = 7;		int theta = 8;

int  black = 0;
int    red = 1;
int   blue = 2;
int yellow = 3;

int entries[10];

ifstream f_fast;
ifstream f_fast_3He;

// #############################################################################################################################################################
void simc_ana(){
        gROOT->Macro( "./root/Greetings.C"  ); // Load greeting message
	gROOT->Macro( "./root/histo_defs.C" ); // Load histograms

	Double_t final_weight;

	// #############################################################################################################################################################
	// OPEN TREE, ACCESS BRANCHES, AND LOAD VARIABLES
	TFile** f_file   = new TFile*[10];
	TTree** t_file   = new TTree*[10];
	TString filename[10];

	filename[0] = "./input/fast_He_bound.root"    ; // fast kin root file 3He bound 
        filename[1] = "./input/fast_He_continuum.root"; // fast kin root file 3He continuum
        filename[2] = "./input/fast_deuteron.root"    ; // fast kin root file scaled d

	filename[3] = "./input/slow_He_bound.root"    ; // slow kin root file 3He bound 
        filename[4] = "./input/slow_He_continuum.root"; // slow kin root file 3He continuum
        filename[5] = "./input/slow_deuteron.root"    ; // slow kin root file scaled d	

	TFile *f_file[0] = new TFile( filename[0] );	TTree *t_file[0] = (TTree*)f_file[0] ->Get("SNT"); 
        TFile *f_file[1] = new TFile( filename[1] ); 	TTree *t_file[1] = (TTree*)f_file[1] ->Get("SNT");
        TFile *f_file[2] = new TFile( filename[2] ); 	TTree *t_file[2] = (TTree*)f_file[2] ->Get("SNT");
	TFile *f_file[3] = new TFile( filename[3] ); 	TTree *t_file[3] = (TTree*)f_file[3] ->Get("SNT");
        TFile *f_file[4] = new TFile( filename[4] ); 	TTree *t_file[4] = (TTree*)f_file[4] ->Get("SNT");
        TFile *f_file[5] = new TFile( filename[5] ); 	TTree *t_file[5] = (TTree*)f_file[5] ->Get("SNT");
	

	for(int i = 0 ; i < 6 ; i++){	
		cout << "Attaching file: " << filename[i] << endl;
		ImportBranch(t_file[i]);
	}	

	// #############################################################################################################################################################
	// FILL HISTOGRAMS	

	// --------------------------------------------------------------------------------
        // 1
	// FAST KINEMATICS IN 3 HELIUM *** BOUND STATE
	int idx = 0;
	for (int i = 0 ; i < entries[idx] ; i++){
		t_file[idx]->GetEntry(i);
                final_weight = Weight[idx]*Normfac[idx]/entries[idx]/3.; // <------- *1/3 to get a quarter of a day
		if(event_passes_cuts(idx)){h1_fast_He -> Fill( Pm[idx] , final_weight );} 
	} //End loop over entries

        // --------------------------------------------------------------------------------
        // 2
	// FAST KINEMATICS IN 3 HELIUM *** CONTINUUM STATE
        int idx = 1;
        for (int i = 0 ; i < entries[idx] ; i++){
                t_file[idx]->GetEntry(i);
                final_weight = Weight[idx]*Normfac[idx]/entries[idx]/3.; // <------- *1/3 to get a quarter of a day
                if(event_passes_cuts(idx)){h1_fast_He -> Fill( Pm[idx] , final_weight );}
        } //End loop over entries

	// --------------------------------------------------------------------------------
	// 3
        // FAST KINEMATICS IN DEUTERON
        int idx = 2;
        for (int i = 0 ; i < entries[idx] ; i++){
                t_file[idx]->GetEntry(i);
		final_weight = 4./3.*Weight[idx]*Normfac[idx]/entries[idx]/3.; // <---------------------------------------
                if(event_passes_cuts(idx)){h1_fast_H2 -> Fill( Pm[idx] , final_weight );} 
        } //End loop over entries

	// --------------------------------------------------------------------------------
        // 4
        // SLOW KINEMATICS IN 3 HELIUM *** BOUND STATE
        int idx = 3;
        for (int i = 0 ; i < entries[idx] ; i++){
                t_file[idx]->GetEntry(i);
                final_weight = Weight[idx]*Normfac[idx]/entries[idx]*4; // <------- *4 to get 4 days
                if(event_passes_cuts(idx)){h1_slow_He -> Fill( Pm[idx] , final_weight );}
        } //End loop over entries

        // --------------------------------------------------------------------------------
        // 5
        // SLOW KINEMATICS IN 3 HELIUM *** CONTINUUM STATE
        int idx = 4;
        for (int i = 0 ; i < entries[idx] ; i++){
                t_file[idx]->GetEntry(i);
                final_weight = Weight[idx]*Normfac[idx]/entries[idx]*4; // <-------- *4 to get 4 days
                if(event_passes_cuts(idx)){h1_slow_He -> Fill( Pm[idx]   , final_weight );} 
        } //End loop over entries

	// --------------------------------------------------------------------------------
        // 6
        // SLOW KINEMATICS IN DEUTERON
        int idx = 5;

        for (int i = 0 ; i < entries[idx] ; i++){
                t_file[idx]->GetEntry(i);
                final_weight = 4./3.*Weight[idx]*Normfac[idx]/entries[idx]*4; // <-------- *4 to get 4 days
                if(event_passes_cuts(idx)){h1_slow_H2 -> Fill( Pm[idx]   , final_weight );}
        } //End loop over entries


	// #############################################################################################################################################################
	// SOME CALCULATIONS
	double integral_fast;
	cout << "" << endl;
	integral_fast = h1_slow_H2 ->Integral();        cout << "Pm slow scaled  d integral = " << integral_fast  <<endl;	
	integral_fast = h1_slow_He ->Integral();        cout << "Pm slow scaled  He integral = " << integral_fast  <<endl;

	// #############################################################################################################################################################
	// FORMAT AND EDIT HISTOGRAMS
	Pretty1D( h1_fast_H2,  black , hadron , pmiss);
	Pretty1D( h1_fast_He, yellow , hadron , pmiss);
	Pretty1D( h1_slow_H2,   blue , hadron , pmiss);
	Pretty1D( h1_slow_He,    red , hadron , pmiss);

	Pretty1D( h1_total_H2,  blue , hadron , pmiss);
        Pretty1D( h1_total_He,   red , hadron , pmiss);

	// #############################################################################################################################################################
	// PLOT HISTOGRAMS
	gStyle->SetOptStat(0);
	
	// ===================================================================
        // MISSING MOMENTUM
        TCanvas *c1 = new TCanvas("c1","",1200,800);
	c1 -> SetBottomMargin(0.12);
	h1_fast_H2 -> Draw();
	h1_fast_He -> Draw("same");
	h1_slow_He -> Draw("same");
	h1_slow_H2 -> Draw("same");
	
	leg = new TLegend(0.6,0.6,0.85,0.85);
	leg->AddEntry(h1_fast_H2,"d fast kinematics 1/3-day");
	leg->AddEntry(h1_fast_He,"He fast kinematics 1/3-day");
	leg->AddEntry(h1_slow_H2,"d slow kinematics 4-day");
	leg->AddEntry(h1_slow_He,"He slow kinematics 4-day");
	leg->Draw("same");

        // ===================================================================
        // MISSING MOMENTUM TOTAL (Added by Or, Sept. 20 2016)
        TCanvas *c2 = new TCanvas("c2","",1200,800);
        c2 -> SetBottomMargin(0.12);
	h1_total_H2->Add(h1_total_H2,h1_fast_H2);
        h1_total_H2->Add(h1_total_H2,h1_slow_H2);
	h1_total_H2->Rebin(5);

        h1_total_He->Add(h1_total_He,h1_fast_He);
        h1_total_He->Add(h1_total_He,h1_slow_He);
        h1_total_He->Rebin(5);

	h1_total_H2->Draw();
        h1_total_He->Draw("same");
	leg1 = new TLegend(0.6,0.6,0.85,0.85);
	leg1->AddEntry(h1_total_H2,"d combined kinematics");
        leg1->AddEntry(h1_total_He,"3He combined kinematics");
        leg1->Draw("same");

        // #############################################################################################################################################################





} //END MAIN FUNCTION
// #############################################################################################################################################################

// ===========================================================================================
// FUNCTION TO EDIT 1D HISTOGRAMS
void Pretty1D(TH1F *gP, int k, int part, int var){
	int color; int style; int linestyle;
	if( k == 0){ color = 1; style = 3004;}
	if( k == 1){ color = 2; style = 3005;}
	if( k == 2){ color = 4; style = 3004;}
	if( k == 3){ color =92; style = 3007;}
	if( k == 4){ color = 7; style = 3006;}
	if( k == 5){ color = 6; style = 3005;}

	AddLabels1D(gP, part, var);

	gP -> SetLineColor(color);
	gP -> SetFillColor(color);
	gP -> SetFillStyle(style);
	gP -> SetMarkerColor(color);
	gP -> GetXaxis()->SetTitleSize(0.05);
	//gP -> GetXaxis()->SetRangeUser(0,5);
	gP -> GetYaxis()->SetLabelSize (0.05);
	gP -> GetXaxis()->SetLabelSize (0.05);
	//gP -> GetXaxis()->SetNdivisions(505);
}

// ===========================================================================================
// FUNCTION TO EDIT 2D HISTOGRAMS
void Pretty2D(TH2F *gP, int xvar, int yvar, int xpart, int ypart){
	AddLabels2D(gP, xvar, yvar, xpart , ypart);
	AddTitle2D (gP, xvar, yvar, xpart , ypart);
	gP -> GetYaxis() ->  SetTitleSize(0.05);
	gP -> GetXaxis() ->  SetTitleSize(0.05);     
	gP -> GetYaxis() ->  SetLabelSize(0.05);
	gP -> GetXaxis() ->  SetLabelSize(0.05);
}


// ===========================================================================================
// FUNCTION TO ADD LABELS TO 1D HISTOGRAMS
void AddLabels1D(TH1F *gP, int part, int var){
	if(var == ytar	){gP->GetXaxis()->SetTitle("y (target) [cm]");}
	
	if(part == hadron){
		if(var == theta	){gP->GetXaxis()->SetTitle("proton in-plane angle [degrees]");}
		if(var == moment){gP->GetXaxis()->SetTitle("proton momentum [MeV/c]");}
	}
	
	if(var == pmiss ){gP->GetXaxis()->SetTitle("|P_{miss}| [GeV/c]");}
}

// ===========================================================================================
// FUNCTION TO ADD LABELS TO 2D HISTOGRAMS
void AddLabels2D(TH2F *gP, int xvar, int yvar, int xpart, int ypart){
	if(xvar == ytar  ){gP->GetXaxis()->SetTitle("y (target) [cm]"          );}
	if(xvar == delta ){gP->GetXaxis()->SetTitle("delta [%]"                );}
	if(xvar == yptar ){gP->GetXaxis()->SetTitle("in-plane angle [deg]"     );}
	if(xvar == xptar ){gP->GetXaxis()->SetTitle("out-of-plane angle [deg]" );}
	if(xvar == moment){gP->GetXaxis()->SetTitle("Momentum [MeV/c]"         );}
	if(xvar == pmiss ){gP->GetXaxis()->SetTitle("|P_{miss}| [GeV/c]"       );}
	if(xvar == emiss ){gP->GetXaxis()->SetTitle("Missing Energy [GeV]"     );}

	if(yvar == ytar  ){gP->GetYaxis()->SetTitle("y (target) [cm]"          );}
	if(yvar == delta ){gP->GetYaxis()->SetTitle("delta [%]"                );}
	if(yvar == yptar ){gP->GetYaxis()->SetTitle("in-plane angle [deg]"     );}
	if(yvar == xptar ){gP->GetYaxis()->SetTitle("out-of-plane angle [deg]" );}
	if(yvar == moment){gP->GetYaxis()->SetTitle("Momentum [MeV/c]"         );}
	if(yvar == pmiss ){gP->GetYaxis()->SetTitle("|P_{miss}| [GeV/c]"       );}
	if(yvar == emiss ){gP->GetYaxis()->SetTitle("Missing Energy [GeV]"     );}
	if(yvar == Qsq   ){gP->GetYaxis()->SetTitle("Q^{2} [(GeV/c)^{2}]"      );}

	if(xpart == hadron){
		if(xvar == moment){gP->GetXaxis()->SetTitle("p_{p} [MeV/c]" );}
	}

	if(ypart == hadron){
		if(yvar == theta){gP->GetYaxis()->SetTitle("#theta_{p} [degrees]" );}
	}

}

// ===========================================================================================
// FUNCTION TO ADD TITLE TO 2D HISTOGRAMS
void AddTitle2D(TH2F *gP, int xvar, int yvar, int xpart, int ypart){
	if((xpart==electron) && (ypart==electron)){gP->SetTitle("electron arm"       );}
	if((xpart==hadron  ) && (ypart==hadron  )){gP->SetTitle("hadron arm"         );}
	if((xpart==hadron  ) && (ypart==electron)){gP->SetTitle("electron vs. hadron");}
	gStyle->SetTitleSize(0.09,"t");
}

// ===========================================================================================
// RETURN TRUE IF H2 IS BIGGER THAN He
bool first_is_bigger( TH1 *d , TH1 *He ){

if( (d->GetMaximum()) > (He->GetMaximum()) ){return 1;}
else return 0;


}

// ===========================================================================================
// TITLES FOR HISTOS WITH CUTS ON MISSING MOMENTUM
void SetThisTitle(TH1 *h, int flag){

		if(flag == 0){h -> SetTitle("Total");}
                if(flag == 1){h -> SetTitle("  0 < Pm <  50 [MeV/c]");}
                if(flag == 2){h -> SetTitle(" 50 < Pm < 100 [MeV/c]");}
                if(flag == 3){h -> SetTitle("100 < Pm < 150 [MeV/c]");}
                if(flag == 4){h -> SetTitle("150 < Pm < 200 [MeV/c]");}
                if(flag == 5){h -> SetTitle("200 < Pm < 250 [MeV/c]");}
                if(flag == 6){h -> SetTitle("250 < Pm < 300 [MeV/c]");}
                if(flag == 7){h -> SetTitle("300 < Pm < 350 [MeV/c]");}

}

















/*

   delta = 100 [p - p(central)] / p(central)

   p = p(central)[1 + delta/100]

 */



