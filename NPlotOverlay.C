// NPlotOverlay.C
// A macro that takes N root files that each store a TH1F and overlays the plots.
// USAGE:
// .L NPlotOverlay.C
// NPlotOverlay("file.root","Histogram of y(x)","x",...);

#include <algorithm>
#include <stdio.h>
#include <string>
#include <regex>

void NPlotOverlay(const char *file, const char *title, const char *variable,
                  bool logscale = false, bool save = false, 
                  string save_name = "", string match = "") {

  std::regex rgx; 
  if (match == "") {
    rgx = ".*";
  } else {
    rgx = (regex) match;
  }


  cout << "Drawing " << title << endl;

  // Set styling.
  if (save) gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0); // 0: disable legend 1: enable legend

  Int_t palette [6] = {2, 4, 3, 418, 5, 433}; // 46, 9 softer palette
  Int_t stipple [6] = {5, 6, 7, 8, 9, 5};
  Int_t style_i = 0;
  Int_t width = 2;
  Double_t kFontSize = 0.04;
  Double_t kFillAlpha = 0.05;

  TCanvas *c = new TCanvas("c", variable, 800, 700);

  if (logscale) c->SetLogy();

  TH1F *h1 = new TH1F("", "", 20, 0., 1.);
  const char *first_hist_name = "";

  TFile *f = new TFile(file);
  TIter next(f->GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    first_hist_name = key->GetName();
    if (!std::regex_match(first_hist_name, rgx)) continue; 

    h1 = (TH1F*)key->ReadObj();
    printf("Found %s object %s.\n", key->GetClassName(),
           first_hist_name);
    break;
  }

  h1->SetTitle(title);
  h1->GetXaxis()->SetTitleSize(kFontSize);
  h1->GetXaxis()->SetTitle(variable);
  h1->GetYaxis()->SetTitleSize(kFontSize);
  h1->GetYaxis()->SetTitle("Events");

  // Set the scale as 130% of the histogram's max value (not for log plots).
  double ceiling = h1->GetMaximum();
  h1->SetMaximum(1.3 * ceiling);
  //if (logscale) h1->SetMaximum(1.);
  Double_t upper_limit = h1->GetBinCenter(h1->FindLastBinAbove());

  h1->SetLineColorAlpha(palette[0], 1);
  h1->SetLineWidth(width);
  h1->SetLineStyle(stipple[0]);
  h1->SetFillColorAlpha(palette[0], kFillAlpha);

  h1->Draw("HIST");

  TLegend *legend = new TLegend(0.6, 0.75, 0.88, 0.88);
  legend->SetBorderSize(0.);
  legend->SetTextSize(kFontSize); // % of pad size
  legend->AddEntry(h1, first_hist_name, "f");

  TIter next2(f->GetListOfKeys());
  TKey *key2;
  while ((key2 = (TKey*)next2())) {
    TClass *cl2 = gROOT->GetClass(key2->GetClassName());
    if (!cl2->InheritsFrom("TH1")) continue;
    if (key2->GetName() == first_hist_name) continue;
    string this_name = key2->GetName();
    if (!std::regex_match(this_name, rgx)) continue;
    printf("Found %s object %s.\n", key->GetClassName(),
           key->GetName());

    TH1F *this_hist = new TH1F("", "", 20, 0., 1.);

    this_hist = (TH1F*)key2->ReadObj();
    style_i++;
    this_hist->SetLineColorAlpha(palette[style_i], 1);
    this_hist->SetFillColorAlpha(palette[style_i], kFillAlpha);
    this_hist->SetLineWidth(width);
    this_hist->SetLineStyle(stipple[style_i]);

    // If a new histogram height is needed, then redraw.
    if (this_hist->GetMaximum() > ceiling && !logscale) {
      ceiling = this_hist->GetMaximum();
      h1->SetMaximum(1.3 * ceiling);
      gPad->Update();
    }
    if (this_hist->GetBinCenter(this_hist->FindLastBinAbove()) > upper_limit) {
      upper_limit = this_hist->GetBinCenter(this_hist->FindLastBinAbove());
    }
    this_hist->Draw("HIST SAME");

    legend->AddEntry(this_hist, key2->GetName(), "f");
  }
  h1->GetXaxis()->SetRangeUser(0,upper_limit);
  legend->Draw();

  if (save) {
    gSystem->Exec("mkdir -p plots");
    gSystem->cd("plots");
    save_name += ".png";
    c->Print(save_name.c_str());
    gSystem->cd("../");
    cout << "Added plots to directory plots/" << endl;
  }
}
