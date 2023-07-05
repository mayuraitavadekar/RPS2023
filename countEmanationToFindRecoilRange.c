void countEmanationToFindRecoilRange()
{
    TFile* input = new TFile("0pct.root", "read");
    TTree* tree = (TTree*) input->Get("particleData");
    int entries = tree->GetEntries();

    Int_t        emanation;
    Double_t        Z;
    Double_t        A;
    Int_t        H2OContent;

    int emanationCount = 0;
    
    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("Z", &Z);
    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("H2OContent", &H2OContent);

    for(int i=0;i<entries;i++)
    {
        tree->GetEntry(i);
        
        if(A == 222 && emanation == 1) emanationCount++;
    }

    printf("particles escaped = %d \n", emanationCount);
}