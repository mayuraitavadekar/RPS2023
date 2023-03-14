
/*

    this macro creates histogram
    of all radioactive particles generated
    in decay chain

    ref:
    https://root.cern.ch/doc/master/classTTree.html

*/


void plotHistXY()
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
    c1->SetGridx(0);
	c1->cd();
 
    // c1->SetFillColor(42);
    // c1->SetGridY();

    
    TFile* file = new TFile("2023-03-14-02-34-46.root", "read");
    
    TTree* tree = (TTree*) file->Get("particleData");
    
    // TString pName;
    Double_t x;
    Double_t y;
    Double_t z;
    Char_t pName[16];
    Double_t A;
    Int_t pid;
    
    // tree->SetBranchAddress("pName", &pName);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("pid", &pid);

    tree->Draw("sqrt(x*x+y*y+z*z)","A==222");

    // for(int i=0;i<tree->GetEntries();i++)
    // {
    //     tree->GetEntry(i);

    //     if(pid == 2)
    //     {
    //         printf("%f %f %f \n",x,y,z);
    //     }
    // }

    // gPad->SetLogy();

    auto htemp = (TH2F*)gPad->GetPrimitive("htemp");
    auto xaxis = htemp->GetXaxis();

    htemp->SetTitle("Distribution of distance of radon rays from origin with air medium");
    htemp->GetYaxis()->SetTitle("count");
	htemp->GetXaxis()->SetTitle("distance (in micrometers)");
    htemp->GetXaxis()->SetTitleOffset(1.5);
    htemp->GetYaxis()->SetTitleColor(kBlack);
    htemp->GetYaxis()->SetTitleSize(0.03);
    htemp->GetXaxis()->SetTitleSize(0.03);
    htemp->GetXaxis()->CenterTitle(true);
    htemp->GetYaxis()->CenterTitle(true);
    htemp->SetLineColor(kBlack);
    gStyle->SetTitleFontSize(0.04);
    htemp->SetLineColor(4);
    // htemp->SetFillColor(2);
    // htemp->SetFillStyle(3004);
    // hist->Draw("");
}

/*

visit TGraph on root documentation

graph->Draw("ALP")  L = line, P = show points
ACP can be used C = smooth curve 

*/
