// script per fittare i dati dell'snr in funzione della corrente.
#include "Riostream.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include <string>
#include <vector>

using namespace std;

// Font/size comuni a tutti i grafici: font code 43 = dimensione assoluta in
// pixel (non frazione dell'altezza del pad), cosi' il testo resta identico
// su pad di altezza diversa (grafico principale vs pannello dei pull) e su
// canvas diverse. Valori pensati per restare leggibili una volta inclusa la
// figura in una pagina A4 di LaTeX.
const Int_t kAxisFont = 43;
const Double_t kTitleSize = 34.;
const Double_t kLabelSize = 30.;
const Double_t kTitleOffsetX = 3.2;
const Double_t kTitleOffsetY = 1.6;
const Double_t kLegendTextSize = 30.;

const Int_t kNPointsPerFile = 20;

// Applica lo stesso font/dimensione a titolo ed etichette di un asse.
void SetAxisStyle(TAxis *axis, Double_t titleOffset) {
  axis->SetTitleFont(kAxisFont);
  axis->SetLabelFont(kAxisFont);
  axis->SetTitleSize(kTitleSize);
  axis->SetLabelSize(kLabelSize);
  axis->SetTitleOffset(titleOffset);
}

// Legge n righe "x y ey" da un file; ritorna false se il file non si apre.
bool ReadXYData(const string &fileName, int n, vector<double> &x,
                 vector<double> &y, vector<double> &ex, vector<double> &ey) {
  ifstream file_in(fileName);
  if (!file_in.is_open()) {
    cout << "Attenzione: Impossibile aprire il file " << fileName << endl;
    return false;
  }

  x.clear();
  y.clear();
  ex.clear();
  ey.clear();

  double current, snr, e_snr;
  for (int j = 0; j < n; j++) {
    file_in >> current >> snr >> e_snr;
    x.push_back(current);
    y.push_back(snr);
    ex.push_back(0.0);
    ey.push_back(e_snr);
  }

  file_in.close();
  return true;
}

// Nome del file senza cartella "rawdata/" e senza estensione ".txt".
string BaseNameNoExt(const string &fileName) {
  size_t start = fileName.find_last_of('/');
  start = (start == string::npos) ? 0 : start + 1;
  size_t ext = fileName.rfind(".txt");
  return fileName.substr(start, ext - start);
}

void Fit_macro_snr_matte() {

  vector<string> fileName = {
      "rawdata/SNR_25kv.txt",
      "rawdata/SNR_40kv.txt",
      "rawdata/SNR_Area.txt",
  };

  TCanvas *c[3];
  TLine *line0[3];

  // Definizione della funzione FUORI dal ciclo.
  // Il range temporaneo (0, 1) verra' aggiornato per ogni file.
  // [0]=A (costante), [1]=B (esponente)
  TF1 *fitfunc = new TF1("fitfunc", "[0]*x**[1]", 0, 1);
  fitfunc->SetLineColor(kOrange + 7);
  fitfunc->SetNpx(1000000);

  for (int i = 0; i < (int)fileName.size(); ++i) {

    vector<double> x, y, ex, ey;
    if (!ReadXYData(fileName[i], kNPointsPerFile, x, y, ex, ey))
      continue;

    int n_points = x.size();
    if (n_points == 0)
      continue;

    bool isArea = (i == 2);
    const char *xTitle = isArea ? "Area [px^{2}]" : "Corrente [mA]";

    // Creazione e stile del TGraphErrors
    TGraphErrors *graph =
        new TGraphErrors(n_points, &x[0], &y[0], &ex[0], &ey[0]);
    graph->SetTitle(Form(";%s;SNR", xTitle));
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kRed + 2);
    graph->SetLineColor(kRed + 2);

    // 1. Aggiornamento del range della funzione globale
    fitfunc->SetRange(x.front(), x.back());

    // 2. Stima dei parametri iniziali
    double guess_A = 10;
    double guess_B = 0.5;
    fitfunc->SetParameters(guess_A, guess_B);

    cout << "Sto leggendo il file: " << fileName[i] << endl;
    cout << "guess_A: " << guess_A << " guess_B: " << guess_B << endl;

    // Esecuzione del fit ("Q" per la modalita' Quiet, la funzione viene
    // automaticamente agganciata al TGraphError)
    graph->Fit(fitfunc, "Q");

    cout << "Fit fatto sul file: " << fileName[i] << endl;
    cout << "Parametro A: " << fitfunc->GetParameter(0) << " ± "
         << fitfunc->GetParError(0) << endl;
    cout << "Parametro B: " << fitfunc->GetParameter(1) << " ± "
         << fitfunc->GetParError(1) << endl;
    cout << "----------------------------------------" << endl;
    cout << "chi-quadro: " << fitfunc->GetChisquare() << endl;
    cout << "numero di gradi di liberta': " << fitfunc->GetNDF() << endl;
    cout << "chi-quadro ridotto: "
         << fitfunc->GetChisquare() / fitfunc->GetNDF() << endl;
    cout << endl << endl;

    // Calcolo dei pull
    TGraph *gPulls = new TGraph(n_points);
    gPulls->SetTitle(Form(";%s;pull", xTitle));
    gPulls->SetMarkerStyle(20);
    gPulls->SetMarkerSize(1);

    for (int j = 0; j < n_points; ++j) {
      double xj = graph->GetPointX(j);
      double yj = graph->GetPointY(j);
      double eyj = graph->GetErrorY(j);
      double pull = (yj - fitfunc->Eval(xj)) / eyj;
      gPulls->SetPoint(j, xj, pull);
    }

    // Disegno: pad superiore col fit, pad inferiore coi pull
    c[i] = new TCanvas(Form("c%d", i), "Fit snr", 1500, 800);
    c[i]->cd();

    TPad *pad1 = new TPad("pad_snr", "snr", 0., 0.28, 1., 1.);
    TPad *pad2 = new TPad("pad_pull", "pull", 0., 0., 1., 0.28);
    pad1->SetBottomMargin(0.00);
    pad1->SetLeftMargin(0.16);
    pad2->SetTopMargin(0.00);
    pad2->SetBottomMargin(0.40);
    pad2->SetLeftMargin(0.16);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    graph->GetXaxis()->SetLabelSize(0.);
    graph->GetXaxis()->SetTitleSize(0.);
    SetAxisStyle(graph->GetYaxis(), kTitleOffsetY);
    graph->Draw("AP");

    TLegend *legend = new TLegend(0.6, 0.2, 0.8, 0.3);
    legend->AddEntry(graph, "Dati", "p");
    legend->AddEntry(fitfunc,
                      Form("Risultato del fit: %.3gx^{%.3g}",
                           fitfunc->GetParameter(0), fitfunc->GetParameter(1)),
                      "l");
    legend->SetBorderSize(0);
    legend->SetTextFont(kAxisFont);
    legend->SetTextSize(kLegendTextSize);
    legend->Draw();
    pad1->Update();

    pad2->cd();
    SetAxisStyle(gPulls->GetXaxis(), kTitleOffsetX);
    SetAxisStyle(gPulls->GetYaxis(), kTitleOffsetY);
    gPulls->GetYaxis()->SetNdivisions(505);

    line0[i] = new TLine(gPulls->GetXaxis()->GetXmin(), 0.,
                          gPulls->GetXaxis()->GetXmax(), 0.);
    line0[i]->SetLineColor(kBlack);
    gPulls->Draw("AP");
    line0[i]->Draw("same");

    c[i]->SaveAs(Form("GRAPH/SNR_fit_%s.pdf", BaseNameNoExt(fileName[i]).c_str()));
  }

  // Grafico di controllo: stessa area del file "SNR_Area.txt" ma con i dati
  // (sbagliati) di "SNR_Area_sbagliata.txt", senza fit.
  string filename2 = "rawdata/SNR_Area_sbagliata.txt";
  vector<double> x2, y2, ex2, ey2;
  if (!ReadXYData(filename2, kNPointsPerFile, x2, y2, ex2, ey2))
    return;

  int n_points2 = x2.size();
  if (n_points2 == 0)
    return;

  TGraphErrors *graph2 =
      new TGraphErrors(n_points2, &x2[0], &y2[0], &ex2[0], &ey2[0]);
  graph2->SetTitle(";Area [px^{2}];SNR");
  graph2->SetMarkerStyle(20);
  graph2->SetMarkerSize(0.8);
  graph2->SetMarkerColor(kPink + 5);

  TCanvas *c2 = new TCanvas("c2", "Fit snr", 1500, 800);
  c2->SetLeftMargin(0.16);
  c2->SetBottomMargin(0.16);
  c2->cd();
  SetAxisStyle(graph2->GetXaxis(), kTitleOffsetX);
  SetAxisStyle(graph2->GetYaxis(), kTitleOffsetY);
  graph2->Draw("AP");

  c2->SaveAs(Form("GRAPH/SNR_fit_%s.pdf", BaseNameNoExt(filename2).c_str()));
}
