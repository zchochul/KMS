void momentum(){

Double_t w = 600;
   Double_t h = 600;
   TCanvas *c = new TCanvas("c", "c", w, h);
   c->Divide(2,2);
   
   TH1D *t1 = new TH1D("p_x","p_x",500,-50,50);
   TH1D *t2 = new TH1D("p_y","p_y",500,-50,50);
   TH1D *t3 = new TH1D("p_z","p_z",500,-50,50);
    ifstream ifile;
    ifile.open("Lab1Momentum.txt");
    double val;
    while(ifile>>val){
        t1->Fill(val);
        ifile>>val;    
        t2->Fill(val);
        ifile>>val;
        t3->Fill(val);
    }
    ifile.close();
    c->cd(1);
    t1->Draw();
    c->cd(2);
    t2->Draw();
    c->cd(3);
    t3->Draw();
}
