void plotEmanationVsH2OContentMultipleGraphs()
{

    // 100 nm
    TFile* input = new TFile("paper-100nm.root", "read");
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
        if(H2OContent == 0 && A == 222 && emanation == 1)       eCounts[0]++;
        else if(H2OContent == 5 && A == 222 && emanation == 1)  eCounts[1]++;
        else if(H2OContent == 10 && A == 222 && emanation == 1) eCounts[2]++;
        else if(H2OContent == 15 && A == 222 && emanation == 1) eCounts[3]++;
        else if(H2OContent == 20 && A == 222 && emanation == 1) eCounts[4]++;
        else if(H2OContent == 25 && A == 222 && emanation == 1) eCounts[5]++;
        else if(H2OContent == 30 && A == 222 && emanation == 1) eCounts[6]++;
        else if(H2OContent == 35 && A == 222 && emanation == 1) eCounts[7]++;
        else if(H2OContent == 40 && A == 222 && emanation == 1) eCounts[8]++;
        else if(H2OContent == 45 && A == 222 && emanation == 1) eCounts[9]++;
        else if(H2OContent == 50 && A == 222 && emanation == 1) eCounts[10]++;
        else if(H2OContent == 55 && A == 222 && emanation == 1) eCounts[11]++;
        else if(H2OContent == 60 && A == 222 && emanation == 1) eCounts[12]++;
        else if(H2OContent == 65 && A == 222 && emanation == 1) eCounts[13]++;
        else if(H2OContent == 70 && A == 222 && emanation == 1) eCounts[14]++;
        else if(H2OContent == 75 && A == 222 && emanation == 1) eCounts[15]++;
        else if(H2OContent == 80 && A == 222 && emanation == 1) eCounts[16]++;
        else if(H2OContent == 85 && A == 222 && emanation == 1) eCounts[17]++;
        else if(H2OContent == 90 && A == 222 && emanation == 1) eCounts[18]++;
        else if(H2OContent == 95 && A == 222 && emanation == 1) eCounts[19]++;
        else if(H2OContent == 100 && A == 222 && emanation == 1)eCounts[20]++;
    }

    for(int i=0;i<=21;i++)
    {
        eCounts[i] = (eCounts[i]/100000)*100;
    }

    // 1000 nm
    input = new TFile("paper-1000nm.root", "read");
    tree = (TTree*) input->Get("particleData");
    entries = tree->GetEntries();

    double eCounts2[21] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("Z", &Z);
    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("H2OContent", &H2OContent);

    for(int i=0;i<entries;i++)
    {
        tree->GetEntry(i);
        if(H2OContent == 0 && A == 222 && emanation == 1)       eCounts2[0]++;
        else if(H2OContent == 5 && A == 222 && emanation == 1)  eCounts2[1]++;
        else if(H2OContent == 10 && A == 222 && emanation == 1) eCounts2[2]++;
        else if(H2OContent == 15 && A == 222 && emanation == 1) eCounts2[3]++;
        else if(H2OContent == 20 && A == 222 && emanation == 1) eCounts2[4]++;
        else if(H2OContent == 25 && A == 222 && emanation == 1) eCounts2[5]++;
        else if(H2OContent == 30 && A == 222 && emanation == 1) eCounts2[6]++;
        else if(H2OContent == 35 && A == 222 && emanation == 1) eCounts2[7]++;
        else if(H2OContent == 40 && A == 222 && emanation == 1) eCounts2[8]++;
        else if(H2OContent == 45 && A == 222 && emanation == 1) eCounts2[9]++;
        else if(H2OContent == 50 && A == 222 && emanation == 1) eCounts2[10]++;
        else if(H2OContent == 55 && A == 222 && emanation == 1) eCounts2[11]++;
        else if(H2OContent == 60 && A == 222 && emanation == 1) eCounts2[12]++;
        else if(H2OContent == 65 && A == 222 && emanation == 1) eCounts2[13]++;
        else if(H2OContent == 70 && A == 222 && emanation == 1) eCounts2[14]++;
        else if(H2OContent == 75 && A == 222 && emanation == 1) eCounts2[15]++;
        else if(H2OContent == 80 && A == 222 && emanation == 1) eCounts2[16]++;
        else if(H2OContent == 85 && A == 222 && emanation == 1) eCounts2[17]++;
        else if(H2OContent == 90 && A == 222 && emanation == 1) eCounts2[18]++;
        else if(H2OContent == 95 && A == 222 && emanation == 1) eCounts2[19]++;
        else if(H2OContent == 100 && A == 222 && emanation == 1)eCounts2[20]++;
    }

    for(int i=0;i<=21;i++)
    {
        eCounts2[i] = (eCounts2[i]/100000)*100;
    }

    // 10 um
    input = new TFile("paper-10um.root", "read");
    tree = (TTree*) input->Get("particleData");
    entries = tree->GetEntries();

    double eCounts3[21] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("Z", &Z);
    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("H2OContent", &H2OContent);

    for(int i=0;i<entries;i++)
    {
        tree->GetEntry(i);
        if(H2OContent == 0 && A == 222 && emanation == 1)       eCounts3[0]++;
        else if(H2OContent == 5 && A == 222 && emanation == 1)  eCounts3[1]++;
        else if(H2OContent == 10 && A == 222 && emanation == 1) eCounts3[2]++;
        else if(H2OContent == 15 && A == 222 && emanation == 1) eCounts3[3]++;
        else if(H2OContent == 20 && A == 222 && emanation == 1) eCounts3[4]++;
        else if(H2OContent == 25 && A == 222 && emanation == 1) eCounts3[5]++;
        else if(H2OContent == 30 && A == 222 && emanation == 1) eCounts3[6]++;
        else if(H2OContent == 35 && A == 222 && emanation == 1) eCounts3[7]++;
        else if(H2OContent == 40 && A == 222 && emanation == 1) eCounts3[8]++;
        else if(H2OContent == 45 && A == 222 && emanation == 1) eCounts3[9]++;
        else if(H2OContent == 50 && A == 222 && emanation == 1) eCounts3[10]++;
        else if(H2OContent == 55 && A == 222 && emanation == 1) eCounts3[11]++;
        else if(H2OContent == 60 && A == 222 && emanation == 1) eCounts3[12]++;
        else if(H2OContent == 65 && A == 222 && emanation == 1) eCounts3[13]++;
        else if(H2OContent == 70 && A == 222 && emanation == 1) eCounts3[14]++;
        else if(H2OContent == 75 && A == 222 && emanation == 1) eCounts3[15]++;
        else if(H2OContent == 80 && A == 222 && emanation == 1) eCounts3[16]++;
        else if(H2OContent == 85 && A == 222 && emanation == 1) eCounts3[17]++;
        else if(H2OContent == 90 && A == 222 && emanation == 1) eCounts3[18]++;
        else if(H2OContent == 95 && A == 222 && emanation == 1) eCounts3[19]++;
        else if(H2OContent == 100 && A == 222 && emanation == 1)eCounts3[20]++;
    }

    for(int i=0;i<=21;i++)
    {
        eCounts3[i] = (eCounts3[i]/100000)*100;
    }


    // 100 um
    input = new TFile("paper-100um.root", "read");
    tree = (TTree*) input->Get("particleData");
    entries = tree->GetEntries();

    double eCounts4[21] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("Z", &Z);
    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("H2OContent", &H2OContent);

    for(int i=0;i<entries;i++)
    {
        tree->GetEntry(i);
        if(H2OContent == 0 && A == 222 && emanation == 1)       eCounts4[0]++;
        else if(H2OContent == 5 && A == 222 && emanation == 1)  eCounts4[1]++;
        else if(H2OContent == 10 && A == 222 && emanation == 1) eCounts4[2]++;
        else if(H2OContent == 15 && A == 222 && emanation == 1) eCounts4[3]++;
        else if(H2OContent == 20 && A == 222 && emanation == 1) eCounts4[4]++;
        else if(H2OContent == 25 && A == 222 && emanation == 1) eCounts4[5]++;
        else if(H2OContent == 30 && A == 222 && emanation == 1) eCounts4[6]++;
        else if(H2OContent == 35 && A == 222 && emanation == 1) eCounts4[7]++;
        else if(H2OContent == 40 && A == 222 && emanation == 1) eCounts4[8]++;
        else if(H2OContent == 45 && A == 222 && emanation == 1) eCounts4[9]++;
        else if(H2OContent == 50 && A == 222 && emanation == 1) eCounts4[10]++;
        else if(H2OContent == 55 && A == 222 && emanation == 1) eCounts4[11]++;
        else if(H2OContent == 60 && A == 222 && emanation == 1) eCounts4[12]++;
        else if(H2OContent == 65 && A == 222 && emanation == 1) eCounts4[13]++;
        else if(H2OContent == 70 && A == 222 && emanation == 1) eCounts4[14]++;
        else if(H2OContent == 75 && A == 222 && emanation == 1) eCounts4[15]++;
        else if(H2OContent == 80 && A == 222 && emanation == 1) eCounts4[16]++;
        else if(H2OContent == 85 && A == 222 && emanation == 1) eCounts4[17]++;
        else if(H2OContent == 90 && A == 222 && emanation == 1) eCounts4[18]++;
        else if(H2OContent == 95 && A == 222 && emanation == 1) eCounts4[19]++;
        else if(H2OContent == 100 && A == 222 && emanation == 1)eCounts4[20]++;
    }

    for(int i=0;i<=21;i++)
    {
        eCounts4[i] = (eCounts4[i]/1000000)*100;
    }

    // 1000 um
    input = new TFile("paper-1000um.root", "read");
    tree = (TTree*) input->Get("particleData");
    entries = tree->GetEntries();

    double eCounts5[21] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    tree->SetBranchAddress("A", &A);
    tree->SetBranchAddress("Z", &Z);
    tree->SetBranchAddress("emanation", &emanation);
    tree->SetBranchAddress("H2OContent", &H2OContent);

    for(int i=0;i<entries;i++)
    {
        tree->GetEntry(i);
        if(H2OContent == 0 && A == 222 && emanation == 1)       eCounts5[0]++;
        else if(H2OContent == 5 && A == 222 && emanation == 1)  eCounts5[1]++;
        else if(H2OContent == 10 && A == 222 && emanation == 1) eCounts5[2]++;
        else if(H2OContent == 15 && A == 222 && emanation == 1) eCounts5[3]++;
        else if(H2OContent == 20 && A == 222 && emanation == 1) eCounts5[4]++;
        else if(H2OContent == 25 && A == 222 && emanation == 1) eCounts5[5]++;
        else if(H2OContent == 30 && A == 222 && emanation == 1) eCounts5[6]++;
        else if(H2OContent == 35 && A == 222 && emanation == 1) eCounts5[7]++;
        else if(H2OContent == 40 && A == 222 && emanation == 1) eCounts5[8]++;
        else if(H2OContent == 45 && A == 222 && emanation == 1) eCounts5[9]++;
        else if(H2OContent == 50 && A == 222 && emanation == 1) eCounts5[10]++;
        else if(H2OContent == 55 && A == 222 && emanation == 1) eCounts5[11]++;
        else if(H2OContent == 60 && A == 222 && emanation == 1) eCounts5[12]++;
        else if(H2OContent == 65 && A == 222 && emanation == 1) eCounts5[13]++;
        else if(H2OContent == 70 && A == 222 && emanation == 1) eCounts5[14]++;
        else if(H2OContent == 75 && A == 222 && emanation == 1) eCounts5[15]++;
        else if(H2OContent == 80 && A == 222 && emanation == 1) eCounts5[16]++;
        else if(H2OContent == 85 && A == 222 && emanation == 1) eCounts5[17]++;
        else if(H2OContent == 90 && A == 222 && emanation == 1) eCounts5[18]++;
        else if(H2OContent == 95 && A == 222 && emanation == 1) eCounts5[19]++;
        else if(H2OContent == 100 && A == 222 && emanation == 1)eCounts5[20]++;
    }

    for(int i=0;i<=21;i++)
    {
        eCounts5[i] = (eCounts5[i]/1000000)*100;
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

    TCanvas *c1 = new TCanvas("c1","c1",300,300);

    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.3);
	c1->SetTopMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetGridy(1);
    c1->SetGridx(1);
	c1->cd();

    // now draw graph
    TMultiGraph *mg = new TMultiGraph();
    
    TGraph* graph1 =  new TGraph(21, moistureLevels, eCounts);
    TGraph* graph2 =  new TGraph(21, moistureLevels, eCounts2);
    TGraph* graph3 =  new TGraph(21, moistureLevels, eCounts3);
    TGraph* graph4 =  new TGraph(21, moistureLevels, eCounts4);
    TGraph* graph5 =  new TGraph(21, moistureLevels, eCounts5);

    mg->Add(graph5);
    mg->Add(graph4);
    mg->Add(graph3);
    mg->Add(graph2);
    mg->Add(graph1);
        
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerColor(2); // red
    graph1->SetLineColor(2);
    graph1->SetLineWidth(2);
    graph1->SetMarkerSize(1.5);
    gStyle->SetTitleFontSize(0.04);

    
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerSize(1.);
    graph2->SetMarkerColor(4); // blue
    graph2->SetLineColor(4);
    graph2->SetLineWidth(2);
    graph2->SetMarkerSize(1.5);
    gStyle->SetTitleFontSize(0.04);

    graph3->SetMarkerStyle(22);
    graph3->SetMarkerSize(1.);
    graph3->SetMarkerColor(6); // magenta
    graph3->SetLineColor(6);
    graph3->SetLineWidth(2);
    graph3->SetMarkerSize(1.5);
    gStyle->SetTitleFontSize(0.04);

    graph4->SetMarkerStyle(33);
    graph4->SetMarkerSize(1.);
    graph4->SetMarkerColor(8); // green
    graph4->SetLineColor(8);
    graph4->SetLineWidth(2);
    graph4->SetMarkerSize(1.5);
    gStyle->SetTitleFontSize(0.04);

    graph5->SetMarkerStyle(47);
    graph5->SetMarkerSize(1.);
    graph5->SetMarkerColor(9); // violet
    graph5->SetLineColor(9);
    graph5->SetLineWidth(2);
    graph5->SetMarkerSize(1.5);
    gStyle->SetTitleFontSize(0.04);

    mg->Draw("ACP");

    mg->GetYaxis()->SetTitle("Radon emanation coefficient (%)");
    mg->GetXaxis()->SetTitle("Moisture content (%)");
    mg->GetXaxis()->SetTitleOffset(1.5);
    mg->GetYaxis()->SetTitleColor(kBlack);
    mg->GetYaxis()->SetTitleSize(0.03);
    mg->GetXaxis()->SetTitleSize(0.03);
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);

    auto legend = new TLegend(0.2, 0.2, .5, .5);
    TString one_micrometer = "1 #mu";
    TString ten_micrometer = "10 #mu";
    TString hundred_micrometer = "100 #mu";

    legend->AddEntry(graph1,"100 nm","lpf");
    legend->AddEntry(graph2,one_micrometer+"m","lpf");
    legend->AddEntry(graph3,ten_micrometer+"m","lpf");
    legend->AddEntry(graph4,hundred_micrometer+"m","lpf");
    legend->AddEntry(graph5,"1 mm","lpf");
    legend->Draw();
    // auto legend = new TLegend(0.2, 0.2, .5, .5);
    // legend->AddEntry(graph,"100 nm (10K Events)","fpl");
    // legend->AddEntry(graph2,"1000 nm (100K Events)","fpl");
    // // legend->AddEntry(graph3,"300 nm (30K Events)","fpl");
    // legend->Draw();
}