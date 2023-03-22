void plotEmanationByH2O()
{
    TFile* input = new TFile("mainfile.root", "read");
    TTree* tree = (TTree*) input->Get("particleData");
    int entries = tree->GetEntries();

    double eCounts[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double count_pct[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double moistureLevels[13] = {15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.};
    
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
        
        if(H2OContent == 15 && A == 222 && emanation == 1)
        {
            eCounts[0]++;
        }

        else if(H2OContent == 20 && A == 222 && emanation == 1)
        {
            eCounts[1]++;
        }

        else if(H2OContent == 25 && A == 222 && emanation == 1)
        {
            eCounts[2]++;
        }

        else if(H2OContent == 30 && A == 222 && emanation == 1)
        {
            eCounts[3]++;
        }

        else if(H2OContent == 35 && A == 222 && emanation == 1)
        {
            eCounts[4]++;
        }

        else if(H2OContent == 40 && A == 222 && emanation == 1)
        {
            eCounts[5]++;
        }

        else if(H2OContent == 45 && A == 222 && emanation == 1)
        {
            eCounts[6]++;
        }

        else if(H2OContent == 50 && A == 222 && emanation == 1)
        {
            eCounts[7]++;
        }

        else if(H2OContent == 55 && A == 222 && emanation == 1)
        {
            eCounts[8]++;
        }

        else if(H2OContent == 60 && A == 222 && emanation == 1)
        {
            eCounts[9]++;
        }

        else if(H2OContent == 65 && A == 222 && emanation == 1)
        {
            eCounts[10]++;
        }

        else if(H2OContent == 70 && A == 222 && emanation == 1)
        {
            eCounts[11]++;
        }

        else if(H2OContent == 75 && A == 222 && emanation == 1)
        {
            eCounts[12]++;
        }
    }

    for(int i=0;i<13;i++)
    {
        printf("ecount[%d] = %f \n", i, eCounts[i]);
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
    TGraph* graph =  new TGraph(13, moistureLevels, eCounts);
    
    graph->GetYaxis()->SetTitle("count of escaped particles from soil grain");
    graph->GetXaxis()->SetTitle("moisture content (%)");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.);
    graph->SetMarkerColor(2); // 4 blue 2 red
    graph->SetLineColor(2);
    graph->SetLineWidth(2);
    graph->SetMarkerSize(1.5);
    graph->GetXaxis()->SetTitleOffset(1.5);
    graph->GetYaxis()->SetTitleColor(kBlack);
    graph->GetYaxis()->SetTitleSize(0.03);
    graph->GetXaxis()->SetTitleSize(0.03);
    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);
    gStyle->SetTitleFontSize(0.04);

    graph->Draw("ALP");

}