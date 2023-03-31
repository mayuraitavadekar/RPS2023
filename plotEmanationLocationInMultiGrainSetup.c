/*
    this macro counts number of radon particles 
    at different locations when they come to rest
    and creates barplot
*/
void plotEmanationLocationInMultiGrainSetup()
{
    gROOT->Reset();

    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(kWhite);
	gStyle->SetFrameFillStyle(1001);
	gStyle->SetFrameFillColor(kWhite);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(kWhite);
	gStyle->SetStatColor(kWhite);
	// gStyle->SetPadTickX(1);
	// gStyle->SetPadTickY(1);

    gStyle->SetLabelColor(1,"X");
	gStyle->SetLabelColor(1,"Y");

    TCanvas *c1 = new TCanvas("c1","c1",800,1000);

    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.3);
	c1->SetTopMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetGridy(0);
    c1->SetGridx(1);
	c1->cd();

    const char *locations[3] = {"Central grain","Sorrounding grain","Pore space"};

    TFile* file = new TFile("emanation-study-10um-1000KEvts-multigrain-recoil-emanation.root", "read");
    
    TTree* tree = (TTree*) file->Get("particleData");
    
    Int_t emanation;
    Double_t A;

    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("A", &A);

    TH1F *h = new TH1F("h","test",3,0,3);

    gPad->SetLogy();

    h->GetYaxis()->SetTitle("count of radon");
	h->GetXaxis()->SetTitle("locations of radon after recoil");
    h->GetXaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->SetTitleColor(kBlack);
    h->GetYaxis()->SetTitleSize(0.03);
    h->GetXaxis()->SetTitleSize(0.03);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    h->SetLineColor(kBlack);
    gStyle->SetTitleFontSize(0.04);
    h->SetLineColor(2);
    h->SetFillColor(2);
    h->SetFillStyle(3004);
    
    int count1=0,count2=0,count3=0;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        if(A == 222 && emanation == 1) {h->Fill(locations[2],1);count1++;}
        else if(A == 222 && emanation == 2) {h->Fill(locations[1],1);count2++;}
        else if(A == 222 && emanation == 3) {h->Fill(locations[0],1);count3++;}
    }

    printf("central grain : %d \n", count3);
    printf("sorrounding grain : %d \n", count2);
    printf("pore space : %d \n", count1);
    
    // h->LabelsDeflate();
    h->Draw();
    TPaveText *pt = new TPaveText(0.7,0.85,0.98,0.98,"brNDC");
    pt->SetFillColor(18);
    pt->SetTextAlign(12);
    pt->AddText("Precise counts in each location");
    pt->AddText("14315 in central grain");
    pt->AddText("592 in surrounding grain");
    pt->AddText("93 in pore space");
    pt->Draw();
    return c1;
}