void temp()
{
    TFile* input = new TFile("2023-03-21-19-41-56.root", "read");
    TTree* tree = (TTree*) input->Get("particleData");
    int entries = tree->GetEntries();

    int eCount = 0;
    
    Int_t        emanation;
    Double_t        Z;
    Double_t        A;
    Int_t        H2OContent;
    
    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("Z", &Z);
    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("H2OContent", &H2OContent);

    for(int i=0;i<entries;i++)
    {
        tree->GetEntry(i);
        
        if(A == 222 && emanation == 1)
        {
            eCount++;
        }
    }

    printf("ecount = %d \n", eCount);
}