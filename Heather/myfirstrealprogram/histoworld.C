void histoworld(double xmin, double xmax ){
  TH1D *hist = new TH1D("HelloHisto", "Basic Histogram", 100, xmin, xmax);
  hist->Draw();

}
