// ******************************************************************
// HISTOGRAM DEFINITION
// by: Reynier Cruz Torres
// Define in this file all the histograms that you want to have
{
int N = 9;
        TH1F** h1_Pp_H2 = new TH1F*[N];
        TH1F** h1_Pp_He = new TH1F*[N];
        TH1F** h1_Th_H2 = new TH1F*[N];
        TH1F** h1_Th_He = new TH1F*[N];

        for(int i = 0 ; i < N ; i++){
                h1_Pp_H2[i] = new TH1F( "" , "" , 20, 1500,  1750);
                h1_Pp_He[i] = new TH1F( "" , "" , 20, 1500,  1750);       
                h1_Th_H2[i] = new TH1F( "" , "" , 20,   46,    52);
                h1_Th_He[i] = new TH1F( "" , "" , 20,   46,    52);
	}

	// Missing Momentum
	h1_fast_He = new TH1F( "fast_He" , "fast_He" ,50,   0,   0.5);
	h1_fast_H2 = new TH1F( "fast_H2" , "fast_H2" ,50,   0,   0.5);
	h1_slow_He = new TH1F( "slow_He" , "slow_He" ,50,   0,   0.5);
        h1_slow_H2 = new TH1F( "slow_H2" , "slow_H2" ,50,   0,   0.5);


        h1_total_He = new TH1F( "total_He" , "slow_He" ,50,   0,   0.5);
        h1_total_H2 = new TH1F( "total_H2" , "slow_H2" ,50,   0,   0.5);






}
// ******************************************************************
