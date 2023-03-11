
/*

    this macro creates histogram
    of all radioactive particles generated
    in decay chain

    ref:


*/


void plotPNames()
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
 
    // c1->SetFillColor(42);
    // c1->SetGridY();

    
    TFile* file = new TFile("2023-03-10-13-13-05-u238-100evts.root", "read");
    
    TTree* tree = (TTree*) file->Get("RDecayProducts");
    
    // TString pName;
    Char_t          pName[16];
    
    tree->SetBranchAddress("pName", &pName);

    TH1F *hist = new TH1F("myhist", "Decay Products in Uranium-238 Series", 66, 1, 66);

    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        hist->Fill(pName, 1);
    }
    
    // graph->SetMarkerStyle(4);
    // graph->SetMarkerSize(1);
    // graph->SetLineColor(2);
    // graph->SetLineWidth(2);
    // graph->SetMarkerColor(4);
    // graph->SetMarkerSize(1.5);
    // // graph->SetMarkerStyle(21);
    // graph->SetTitle("Radon Counts Vs pH");
    // graph->GetXaxis()->SetTitle("pH");
    // graph->GetYaxis()->SetTitle("Radon Counts");

    gPad->SetLogy();

    hist->GetYaxis()->SetTitle("count");
	// hist->GetXaxis()->SetTitle("Decay Products of Uranium 238");
    hist->GetYaxis()->SetTitleColor(kBlack);
    hist->GetYaxis()->SetTitleSize(0.03);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetXaxis()->CenterTitle(true);
    hist->GetYaxis()->CenterTitle(true);
    hist->SetLineColor(kBlack);
    gStyle->SetTitleFontSize(0.04);
    // gStyle->SetTitleFont(13);
    hist->SetLineColor(2);
    hist->SetFillColor(2);
    hist->SetFillStyle(3004);
    hist->Draw("");
    
    // c1->Update();
    // c1->GetFrame()->SetFillColor(21);
    // c1->GetFrame()->SetBorderSize(12);
    // c1->Modified();
}

/*

visit TGraph on root documentation

graph->Draw("ALP")  L = line, P = show points
ACP can be used C = smooth curve 

*/
