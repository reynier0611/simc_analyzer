// ===========================================================================================
// IMPORT BRANCH FROM ROOT TREE
void ImportBranch(TTree *T){
        entries[counter1] = T -> GetEntries();

        T ->SetBranchAddress("Normfac" ,&Normfac[counter1]); // Normfac
        T ->SetBranchAddress("Weight"  ,&Weight [counter1]); // Weight given to each event (related to cross section, Jacobian, ...)
        T ->SetBranchAddress("Em"      ,&Em     [counter1]); // Missing Energy (GeV?)
        T ->SetBranchAddress("Pm"      ,&Pm     [counter1]); // Missing Momentum (GeV/c)
        T ->SetBranchAddress("PmPar"   ,&PmPar  [counter1]); // Component of the Missing momentum which is parallel to \vec{q}
        T ->SetBranchAddress("Q2"      ,&Q2     [counter1]); // Q^2 ((GeV/c)^2)

        // ********************************************************
        //Load Electron Variables
        T->SetBranchAddress("e_yptar"   ,&e_yptar [counter1]);  // electron in-plane angle        [rad]           reconstructed
        T->SetBranchAddress("e_xptar"   ,&e_xptar [counter1]);  // electron out-of-plane angle    [rad]           reconstructed

        // ********************************************************
        //Load Hadron Variables
        T->SetBranchAddress("h_ytar"    ,&h_ytar  [counter1]);  // hadron y-target              [cm]            reconstructed
        T->SetBranchAddress("h_delta"   ,&h_delta [counter1]);  // hadron momentum fraction     [%]             reconstructed
        T->SetBranchAddress("h_yptar"   ,&h_yptar [counter1]);  // hadron in-plane angle        [rad]           reconstructed
        T->SetBranchAddress("h_xptar"   ,&h_xptar [counter1]);  // hadron out-of-plane angle    [rad]           reconstructed
        T->SetBranchAddress("h_ytari"   ,&h_ytari [counter1]);  // hadron y-target              [cm]            generated
        T->SetBranchAddress("h_deltai"  ,&h_deltai[counter1]);  // hadron momentum fraction     [%]             generated
        T->SetBranchAddress("h_yptari"  ,&h_yptari[counter1]);  // hadron in-plane angle        [rad]           generated
        T->SetBranchAddress("h_xptari"  ,&h_xptari[counter1]);  // hadron out-of-plane angle    [rad]           generated

        T->SetBranchAddress("h_pf"      ,&h_pf    [counter1]);  // outgoing proton momentum
        T->SetBranchAddress("h_Thf"     ,&h_Thf   [counter1]);  // outgoing proton in-plane angle

        counter1++;

}

// ===========================================================================================
// RETURN TRUE IF EVENT PASSES CUTS
bool event_passes_cuts(int idx){

if(acos(-PmPar[idx] / Pm[idx]) < 0.698){                // Cut to reduce FSI
if((h_yptar[idx] < 0.028) && (h_xptar[idx] < 0.060)){   // Cut on   hadron-spectrometer acceptance
if((e_yptar[idx] < 0.028) && (e_xptar[idx] < 0.060)){   // Cut on electron-spectrometer acceptance

        return 1;

}                                                       // Cut on electron-spectrometer acceptance
}                                                       // Cut on   hadron-spectrometer acceptance
}                                                       // Cut to reduce FSI

else return 0;

}

