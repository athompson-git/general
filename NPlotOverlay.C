// NPlotOverlay.C
// A macro that takes N root files that each store a TH1F and overlays the plots.
// USAGE:
// .x NPlotOverlay.C("file.root", "title", "Pt [GeV]", <yheight>)

void NPlotOverlay(const char *file, const char *title, const char *variable,
                  Double_t yheight) {
  gStyle->SetOptStat(0); // 0: disable legend 1: enable legend

  Int_t palette [5] = {46, 9, 30, 40, 38};
  Int_t color = 0;

  TCanvas *c = new TCanvas("c", variable, 800, 700);
  TH1F *h1 = new TH1F("", "", 20, 0., 1.);
  const char *first_hist_name = "";

  TFile *f = new TFile(file);
  TIter next(f->GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    cout << "Found object: " << key->GetClassName() << endl;
    if (!cl->InheritsFrom("TH1F")) continue;
    h1 = (TH1F*)key->ReadObj();
    std::cout << key->GetTitle() << endl;
    first_hist_name = key->GetName();
    break;
  }

  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(variable);
  h1->GetYaxis()->SetTitle("a.u.");
  h1->GetYaxis()->SetRangeUser(0., yheight);

  h1->SetLineColorAlpha(palette[0], 1);
  h1->SetLineWidth(4);
  h1->SetFillColorAlpha(palette[0], 0.1);

  h1->Draw("HIST");

  TLegend *legend = new TLegend(0.1, 0.80, 0.4, 0.90);
  legend->SetBorderSize(0.);
  legend->SetTextSize(.03); // % of pad size
  legend->AddEntry(h1, first_hist_name, "f");

  TIter next2(f->GetListOfKeys());
  TKey *key2;
  while ((key2 = (TKey*)next2())) {
    TClass *cl2 = gROOT->GetClass(key2->GetClassName());
    cout << "Found object: " << key2->GetClassName() << endl;
    if (!cl2->InheritsFrom("TH1F")) continue;
    if (key2->GetName() == first_hist_name) continue;
    TH1F *this_hist = new TH1F("", "", 10, 0., 1.);

    this_hist = (TH1F*)key2->ReadObj();

    color++;
    this_hist->SetLineColorAlpha(palette[color], 1);
    this_hist->SetFillColorAlpha(palette[color], 0.1);
    this_hist->SetLineWidth(4);

    this_hist->Draw("HIST SAME");

    legend->AddEntry(this_hist, key2->GetName(), "f");
  }

  legend->Draw();
}
