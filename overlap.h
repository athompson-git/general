// A function that takes in two TH1's and returns the computation
// of their overlap.
#include <algorithm>

Double_t overlap(TH1 *h1, TH1 *h2) {

  Int_t nbins_1 = h1->GetXaxis()->GetNbins();
  Int_t nbins_2 = h2->GetXaxis()->GetNbins();
  Int_t length = std::min(nbins_1, nbins_2); // to play it safe if one bin count is bigger.
  // also, if there is an upper bound for bin size then it implies an upper bound for the overlap
  // since we assume the smaller histogram is zero outside its binned interval.

  Double_t overlap = 0.0; // -1 to make buggy behavior visible

  for (Int_t i = 0; i < length; i++) {
    Double_t ith_bin_1 = h1->GetBinContent(i);
    Double_t ith_bin_2 = h2->GetBinContent(i);
    Double_t inner_product_i = ith_bin_1 * ith_bin_2;
    overlap += inner_product_i;
  }

  return overlap;
}
