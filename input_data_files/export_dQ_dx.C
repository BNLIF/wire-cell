void export_dQ_dx(){
    TFile *file = new TFile("stopping_ave_dQ_dx.root");
    TGraph *g_muon = (TGraph*)file->Get("muon");
    if (!g_muon) {
        std::cout << "Failed to get TGraph 'muon' from file." << std::endl;
        return;
    }
    int n = g_muon->GetN();
    double x, y;
    for (int i = 0; i < n; ++i) {
        g_muon->GetPoint(i, x, y);
        std::cout << "Point " << i << ": x = " << x << ", y = " << y << std::endl;
    }
}