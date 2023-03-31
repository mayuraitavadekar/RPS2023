void plotCountInSurroundingGrainsByH2O()
{
    TFile* input = new TFile("emanation-study-100nm-10Kevts-multigrain-recoil-emanation.root", "read");
    TTree* tree = (TTree*) input->Get("particleData");
    int entries = tree->GetEntries();

    double eCounts[21] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double moistureLevels[21] = {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.};
    
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
        
        if(H2OContent == 0 && A == 222 && emanation == 2)       eCounts[0]++;
        else if(H2OContent == 5 && A == 222 && emanation == 2)  eCounts[1]++;
        else if(H2OContent == 10 && A == 222 && emanation == 2) eCounts[2]++;
        else if(H2OContent == 15 && A == 222 && emanation == 2) eCounts[3]++;
        else if(H2OContent == 20 && A == 222 && emanation == 2) eCounts[4]++;
        else if(H2OContent == 25 && A == 222 && emanation == 2) eCounts[5]++;
        else if(H2OContent == 30 && A == 222 && emanation == 2) eCounts[6]++;
        else if(H2OContent == 35 && A == 222 && emanation == 2) eCounts[7]++;
        else if(H2OContent == 40 && A == 222 && emanation == 2) eCounts[8]++;
        else if(H2OContent == 45 && A == 222 && emanation == 2) eCounts[9]++;
        else if(H2OContent == 50 && A == 222 && emanation == 2) eCounts[10]++;
        else if(H2OContent == 55 && A == 222 && emanation == 2) eCounts[11]++;
        else if(H2OContent == 60 && A == 222 && emanation == 2) eCounts[12]++;
        else if(H2OContent == 65 && A == 222 && emanation == 2) eCounts[13]++;
        else if(H2OContent == 70 && A == 222 && emanation == 2) eCounts[14]++;
        else if(H2OContent == 75 && A == 222 && emanation == 2) eCounts[15]++;
        else if(H2OContent == 80 && A == 222 && emanation == 2) eCounts[16]++;
        else if(H2OContent == 85 && A == 222 && emanation == 2) eCounts[17]++;
        else if(H2OContent == 90 && A == 222 && emanation == 2) eCounts[18]++;
        else if(H2OContent == 95 && A == 222 && emanation == 2) eCounts[19]++;
        else if(H2OContent == 100 && A == 222 && emanation == 2)eCounts[20]++;
    }

    for(int i=0;i<22;i++)
    {
        printf("ecount[%d] = %f \n", i, eCounts[i]);
        eCounts[i] = (eCounts[i]/10000)*100;
    }

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
    c1->SetGridy(1);
    c1->SetGridx(1);
	c1->cd();

    // now draw graph
    TGraph* graph =  new TGraph(21, moistureLevels, eCounts);
    
    graph->GetYaxis()->SetTitle("radon embedded in surrounding grains (%)");
    graph->GetXaxis()->SetTitle("moisture content (%)");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.);
    graph->SetMarkerColor(6); // 4 blue 2 red
    graph->SetLineColor(6);
    graph->SetLineWidth(2);
    graph->SetMarkerSize(1.5);
    graph->GetXaxis()->SetTitleOffset(1.5);
    graph->GetYaxis()->SetTitleColor(kBlack);
    graph->GetYaxis()->SetTitleSize(0.03);
    graph->GetXaxis()->SetTitleSize(0.03);
    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);
    gStyle->SetTitleFontSize(0.04);

    graph->Draw("ACP");

    auto legend = new TLegend(0.2, 0.2, .5, .5);
    legend->AddEntry(graph,"100 nm (10K Events)","l");
    legend->Draw();
}