void export_range(){
    TFile *file = new TFile("ave_range_to_kenergy.root");
    TGraph *g_muon = (TGraph*)file->Get("proton");
    if (!g_muon) {
        std::cout << "Failed to get TGraph 'muon' from file." << std::endl;
        return;
    }
    int n = g_muon->GetN();
    std::cout << "X: " << std::endl;
    double x, y;
    for (int i = 0; i < n; ++i) {
        g_muon->GetPoint(i, x, y);
        // std::cout << "Point " << i << ": x = " << x << ", y = " << y << std::endl;
        std::cout << x << ", ";
        if (i%20==19) std::cout << std::endl;
    }

    std::cout << std::endl << "Y: " << std::endl;
    for (int i = 0; i < n; ++i) {
        g_muon->GetPoint(i, x, y);
        // std::cout << "Point " << i << ": x = " << x <<
        std::cout << y << ", ";
        if (i%20==19) std::cout << std::endl;
    }
}