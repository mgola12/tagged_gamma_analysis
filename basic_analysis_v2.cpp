// basic_analysis_v2.cpp
// Beamline + hodoscope data-quality and efficiency analysis.
// Reads a calibrated flat ntuple (beam_monitor_calib tree) from the ROOT file
// produced by the ToolDAQ beamline readout.
//
// Analyses performed per run:
//   1. QDC and TDC distributions for every beamline channel.
//   2. TOF spectrum (T0 − TOF) with a double-Gaussian fit (electrons + protons).
//   3. Beam momentum spectrum from TOF in the proton gate.
//   4. Hodoscope hit multiplicity and occupancy.
//   5. Missing-HD diagnostics (events with lead-glass or TOF hits but no HD).
//   6. Hodoscope timing relative to T2 (single-hit events on HD14/PMT23).
//   7. Lead-glass charge spectrum under single-HD selection.
//
// Usage:
//   ./basic_analysis <run_number> [channel_map.json]
//   ./basic_analysis <path/to/rootfile> [channel_map.json]

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TSystem.h>
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "constants.h"
#include "buffer.h"
#include "event.h"
#include "utils.h"
#include "TMath.h"


int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
    return 1;
  }

  std::string json_file("TDC_and_QDC_channels.json");
  if (argc==3) {
    json_file = std::string(argv[2]);
  }
  std::cout << "using channel map " << json_file << std::endl;

  TFile *file = nullptr;
  std::string run_str = "";

  int runNumber = -1;
  // Convert the argument to an integer (run number)                                                                                                                                                                  
  try {
    runNumber = std::stoi(argv[1]);

    run_str = std::to_string(runNumber);
    if (runNumber < 1000)
      run_str = "0" + run_str;


    // Build the path relative to the current working directory.
    // The ntuple files are named beamline_run<RRRR>_tuple_calib.root.
    std::string fileName = "beamline_run" + run_str + "_tuple_calib.root";

    file = TFile::Open(fileName.c_str());

    // Check if the file opened successfully                                                                                                                                                        
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open the ROOT file!" << std::endl;
      return 1;
    }
  } // try                                                                                                                                                                                                
  catch (...) {
    std::cout << "Will try to open ROOT file as passed as 1st argument." << std::endl;
    file = TFile::Open(argv[1]);
    run_str = extractAfterSubstring(std::string(argv[1]), "_run", 4);
    runNumber = std::stoi(run_str);
    std::cout << "Assuming run number " << runNumber << std::endl;
  }

  //  gSystem->Exec(("mkdir -p ./basic_analysis_plots/run" + run_str).c_str());
  //  gSystem->Exec("mkdir -p ./basic_analysis_histos");

  TTree *tree = (TTree*)file->Get("beam_monitor_calib");
  if (tree) {
    // Do something with the tree, like print its structure                                                                                                                                     
    tree->Print();
  }

    // Branch pointers
    std::vector<int>* ids = nullptr;
    std::vector<std::string>* id_names = nullptr;
    std::vector<double>* qdc = nullptr;
    std::vector<std::vector<double>>* tdc = nullptr;

    tree->SetBranchAddress("beamline_id", &ids);
    tree->SetBranchAddress("beamline_id_name", &id_names);
    tree->SetBranchAddress("beamline_qdc_charge", &qdc);
    tree->SetBranchAddress("beamline_tdc_time", &tdc);

    std::map<int, TH1D*> qdc_hists;
    std::map<int, TH1D*> tdc_hists;


    nlohmann::json configData;
    std::ifstream inFile(json_file);
    if (inFile.is_open()) {
      inFile >> configData;
      inFile.close();
    }
    
    
    //    TH1I* h_hodo_mult = new TH1I("h_hodo_mult", "Hodoscope Hit Multiplicity;# PMTs with hits;Events", 16, -1, 15);
    TH1I* h_hodo_mult = new TH1I("h_hodo_mult", ";# PMTs with hits;Number of Events", 
				 8, -0.5, 7.5);

    TH1D* h_qdc_multi_hodo = new TH1D("h_qdc_multi_hodo",
				      "QDC of Hodoscope PMTs in events with >2 hodoscope hits;QDC [a.u.];Entries", 200, 0, 4000);
    
    TH1I* h_hodo_channels_multi = new TH1I("h_hodo_channels_multi",
					   "Hodoscope PMTs hit in events with >2 hodoscope hits;PMT ID;Counts", 24, 15.5, 39.5);

    TH1D* h_multiple_photon_candidates = new TH1D("h_multiple_photon_candidates",
						  "LG QDC for multi-photon candidates (QDC > 500 & TOF hit);QDC [a.u.];Entries", 200, 0, 4000);
    
    TH1D* h_qdc_no_tdc_LG = new TH1D("h_qdc_no_tdc_LG",
				     "LG QDC when TDC is missing but TOF present;QDC [a.u.];Entries", 200, 0, 4000);
    
    TH1D* h_TOF_timing_no_HD_LG = new TH1D("h_TOF_timing_no_HD_LG",
					   "TOF timing with missing HD + LG TDC;TOF [ns];Entries", 400, -200, 200);
    
    TH1D* h_missing_HD_with_LG = new TH1D("h_missing_HD_with_LG",
					  "Events with LG + TOF hits but missing HD;LG QDC;Entries", 200, 0, 4000);
    
    TH2D* h2_missing_vs_momentum = new TH2D("h2_missing_vs_momentum",
					    "Missing HD vs momentum (placeholder);Momentum [MeV/c];QDC [a.u.]", 8, 450, 1250, 100, 0, 4000);
    
    TH1I* h_missing_other = new TH1I("h_missing_other",
				     "Missing LG/TOF categories;Category;Events", 9, -0.5, 8.5);

    TH1F* h_time_hd11_T2 = new TH1F("h_time_hd11_T2",
			      "Timing Distribution for Hodoscope channel Relative to T2; Time - T2 (ns); Counts", 600, -300, 300);
    
    TH1F* h_qdc_pmt25 = new TH1F("h_qdc_pmt25",
				 "QDC Charge Distribution for PbG (After Cuts); QDC Charge; Counts", 100, 0, 1500);

    TH2I* h_hodo_hits_2D = new TH2I("h_hodo_hits_2D",  "Hodoscope Hits vs Event Number;Event Number;Hodoscope PMT ID", 81779, 0, 81779, 24, 15.5, 39.5);

    TH1I* h_pmt_missing_tdc = new TH1I("h_pmt_missing_tdc", "Hodoscope PMTs with QDC but Missing TDC;PMT ID;Counts", 24, 15.5, 39.5);

    TH1F* h_qdc_LG = new TH1F("h_qdc_LG",
                          "Lead Glass QDC Spectrum (Clean Positron Beam); QDC Charge; Counts", 
                          200, 0, 4000);

    TH2D* h2_qdc_vs_tdiff = new TH2D("h2_qdc_vs_tdiff", "Lead Glass QDC vs. TDC Time Diff;TDC_{LG} - <TDC_{T0}> (coarse);QDC_{LG} (a.u.)", 100, 50, 150, 200, 0, 2000);
    

    
    TH1D* hTOF = new TH1D("hTOF", "((TOF - TDCT01) - (T0 - TDCT00));TOF (ns);Counts", 150, 0, 150);
    TH1D* hMomentum = new TH1D("hMomentum", "Beam Momentum;Beam Momentum (MeV);Counts", 100, 0, 1500);

    TH2D* hTOFvsP = new TH2D("hTOFvsP", 
    "Proton Momentum vs TOF (T0 - TOF);TOF (ns);Momentum (MeV/c)", 
    100, -100, 0,     // x-axis: TOF (T0 - TOF), range can be adjusted
    100, 0, 1000      // y-axis: Momentum, adjust based on expected range
			     );
    
    
    Long64_t nentries = tree->GetEntries();
    std::map<int, int> hit_counts;
    int events_without_hits_HD = 0;

    // Channel-index constants (from beamline mapping JSON)
    const std::vector<int> t0_channels  = {0, 1, 2, 3};
    const std::vector<int> tof_channels = {48, 49, 50, 51, 52, 53,
                                            54, 55, 56, 57, 58, 59,
                                            60, 61, 62, 63};
    // TOF geometry / calibration constants
    const double clight            = 299792458.0;   // m/s
    const double m_p               = 938.272;       // MeV/c²
    const double distance_m        = 6.53;          // m, T0–TOF flight path
    const double energy_loss       = 6.0;           // MeV, material correction
    const double beamline_length_cm = 653.0;
    const double c_cm_per_ns       = 29.9792;
    const double electron_peak_measured_ns = 36.98; // ns, from data

    for (Long64_t i = 0; i < nentries; ++i) {
      // Load all branch data for this entry FIRST.
      tree->GetEntry(i);
      
      // ── TOF reconstruction ─────────────────────────────────────────────────
      double t0_time   = 1e9, tof_time  = 1e9;
      double tdct00_time = 1e9, tdc01_time = 1e9;
      bool t0_found = false, tof_found = false;
      bool tdct00_found = false, tdc01_found = false;

      for (int ch : t0_channels) {
        if ((size_t)ch < tdc->size()) {
          for (double hit : tdc->at(ch)) {
            if (hit < t0_time) { t0_time = hit; t0_found = true; }
          }
        }
      }
      for (int ch : tof_channels) {
        if ((size_t)ch < tdc->size()) {
          for (double hit : tdc->at(ch)) {
            if (hit < tof_time) { tof_time = hit; tof_found = true; }
          }
        }
      }
      if (31 < (int)tdc->size()) {
        for (double hit : tdc->at(31))
          if (hit < tdct00_time) { tdct00_time = hit; tdct00_found = true; }
      }
      if (46 < (int)tdc->size()) {  // TDC trigger channel (TDC1 common stop)
        for (double hit : tdc->at(46))
          if (hit < tdc01_time) { tdc01_time = hit; tdc01_found = true; }
      }

      double tof_ns = 1e9;
      if (t0_found && tof_found && tdct00_found && tdc01_found) {
        tof_ns = (tof_time - tdc01_time) - (t0_time - tdct00_time);
        hTOF->Fill(tof_ns);
      }

      // Correct for the measured electron-peak offset
      const double electron_peak_expected_ns = beamline_length_cm / c_cm_per_ns;
      const double tof_offset       = electron_peak_measured_ns - electron_peak_expected_ns;
      double tof_ns_corrected = tof_ns - tof_offset;

      if (tof_ns_corrected > 57  && tof_ns_corrected < 65) {
      //      if (tof_ns_corrected > 70  && tof_ns_corrected < 90) {
	double tof_s = tof_ns_corrected * 1e-9;
	double beta = distance_m / (tof_s * clight);
        if (beta <= 0 || beta >= 1) continue;
        double gamma = 1.0 / sqrt(1 - beta * beta);
        double T_meas = (gamma - 1.0) * m_p;
        double T_corr = T_meas + energy_loss;
	//double p = m_p * beta * gamma;
	//double p = sqrt(T_meas * T_meas + 2.0 * T_meas * m_p);
	double p = sqrt(T_corr * T_corr + 2.0 * T_corr * m_p); // MeV/c
	hMomentum->Fill(p);
	//hTOFvsP->Fill(tof_ns, p); 
      }     
      

      
      // ── Beamline channel loop: QDC / TDC histograms + hodoscope selection ──
	bool has_hit_in_range = false;
	int hodo_hits_in_event = 0;
	double time_relative_to_T2 = NAN;
	
        for (size_t ch = 0; ch < ids->size(); ++ch) {
            int id = ids->at(ch);
	    bool is_hodo = (id >= 16 && id <= 23) || (id >= 32 && id <= 38);

	    if (!tdc->at(ch).empty()) {
	      hit_counts[id]++;
	      if (is_hodo) {
		hodo_hits_in_event++; 
		has_hit_in_range = true;
	      }
	    }
            // QDC histogram
            if (qdc_hists.find(id) == qdc_hists.end()) {
                qdc_hists[id] = new TH1D(Form("qdc_id%d", id), Form("QDC for ID %d;Charge [a.u.];Entries", id), 200, 0, 4000);
            }
            qdc_hists[id]->Fill(qdc->at(ch));

            // TDC histogram
            if (tdc_hists.find(id) == tdc_hists.end()) {
                tdc_hists[id] = new TH1D(Form("tdc_id%d", id), Form("TDC for ID %d;Time [ns];Entries", id), 400, -200, 200);
            }
            for (const auto& tval : tdc->at(ch)) {
                tdc_hists[id]->Fill(tval);
            }
        }//PMT LOOP END
       


	if (!has_hit_in_range) {


	  ////////

	  for (size_t ch = 0; ch < ids->size(); ++ch) {
	    int id = ids->at(ch);
	    bool is_hodo = (id >= 16 && id <= 23) || (id >= 32 && id <= 38);
	    if (!is_hodo) continue;
	    
	    bool has_qdc = (qdc->size() > ch && qdc->at(ch) > 50);   // threshold adjustable
	    bool has_tdc = (tdc->size() > ch && !tdc->at(ch).empty());
	    
	    if (has_qdc && !has_tdc) {
	      //std::cout << "  → PMT " << id << " has QDC = " << qdc->at(ch) << " but **NO TDC**\n";
	      h_pmt_missing_tdc->Fill(id);
	    }
	  }


	  ///////

	  events_without_hits_HD++;

	  bool has_TDC_LG = !tdc->at(40).empty();
	  bool has_T2_hits = !tdc->at(8).empty();
	  double qdc40 = qdc->at(40);
	  bool has_HC2_hit = (tdc->size() > 11 && !tdc->at(11).empty());

	  if (has_T2_hits && qdc40 > 500) {
            h_multiple_photon_candidates->Fill(qdc40);
	  }

	  if (!has_TDC_LG && has_T2_hits) {
            h_qdc_no_tdc_LG->Fill(qdc40);

            if (!tdc->at(8).empty()) {
	      double tof = tdc->at(8)[0];
	      h_TOF_timing_no_HD_LG->Fill(tof);
            }
	  }

	  bool has_leadglass_hit = (qdc40 > 200);
	  if (has_leadglass_hit && has_T2_hits) {
            h_missing_HD_with_LG->Fill(qdc40);
            h2_missing_vs_momentum->Fill(1000, qdc40);
	  } else {
            int category = 0;
            if (!has_leadglass_hit) category += 1;
            if (!has_T2_hits) category += 2;
	    if (!has_HC2_hit) category += 4;
            h_missing_other->Fill(category);
	  }
	}
	  	  
    
	h_hodo_mult->Fill(hodo_hits_in_event);


	if (has_hit_in_range && hodo_hits_in_event > 2) {
	  //std::cout << "Event " << i << " has " << hodo_hits_in_event << " hodoscope hits. Channel IDs:\n";

	  for (size_t ch = 0; ch < ids->size(); ++ch) {
	    int id = ids->at(ch);
	    bool is_hodo = (id >= 16 && id <= 23) || (id >= 32 && id <= 38);
	    
	    if (is_hodo && !tdc->at(ch).empty()) {
	      //std::cout << "   - PMT ID: " << id << ", QDC: " << qdc->at(ch) << "\n";
	      h_qdc_multi_hodo->Fill(qdc->at(ch));
	      h_hodo_channels_multi->Fill(id);
	      h_hodo_hits_2D->Fill(i, id);
	    }
	  }
	}


	
	const int HD_pmt = 23;
	if (tdc->size() > HD_pmt && tdc->size() > 8) {

	  size_t n_common_hits = std::min(tdc->at(HD_pmt).size(), tdc->at(8).size());
	  
	  for (size_t hit_idx = 0; hit_idx < n_common_hits; ++hit_idx) {
	    double t16 = tdc->at(HD_pmt)[hit_idx];
	    double t8  = tdc->at(8)[hit_idx];
	    time_relative_to_T2 = t16 - t8;
	    h_time_hd11_T2->Fill(time_relative_to_T2);
	  }
	}


	
	int hit_pmt_count = 0;
	bool only_hd_pmt_hit = false;
	
	for (int pmti = 16; pmti < 39; pmti++) {
	  if ((pmti >= 16 && pmti <= 23) || (pmti >= 32 && pmti <= 38)){
	    if (tdc->size() > pmti && !tdc->at(pmti).empty()) {
		hit_pmt_count++;
		if (pmti == HD_pmt) {
		  only_hd_pmt_hit = true;
		}
	      }
	    }
	}
	double qdc_value_hodoscope = (qdc->size() > HD_pmt) ? qdc->at(HD_pmt) : -999;

	if (only_hd_pmt_hit && hit_pmt_count == 1){
	  if (std::isnan(time_relative_to_T2)) continue;
	  if (time_relative_to_T2 < 110 || time_relative_to_T2 > 170) continue;
	  if (qdc_value_hodoscope < 100 || qdc_value_hodoscope > 700) continue;
	  h_qdc_pmt25->Fill(qdc->at(40));
	  
	  
	  }
    } //EVENT LOOP END

    std::cout << "\nPMT Hit Summary:\n";
    for (const auto& [id, count] : hit_counts) {
      std::cout << "PMT ID " << id << " was hit " << count << " times\n";  
    }

    double efficiency = 1.0 - (static_cast<double>(events_without_hits_HD) / nentries);

    std::cout << "Total number of events: " << nentries << std::endl;
    std::cout << "Number of events without hits in Hodoscope: " << events_without_hits_HD << std::endl;
    std::cout << "Efficiency: " << efficiency * 100 << "%" << std::endl;
    
    std::string output_dir = "plots/run" + run_str;
    gSystem->mkdir(output_dir.c_str(), true);

    for (auto& [id, hist] : qdc_hists) {
    TCanvas* c = new TCanvas();
    c->SetLogy();  // Enable log scale on Y-axis
    hist->Draw();
    c->SaveAs(Form("%s/qdc_id%d.png", output_dir.c_str(), id));
    delete c;
    }

    // Save TDC histograms
    for (auto& [id, hist] : tdc_hists) {
      TCanvas* c = new TCanvas();
      c->SetLogy();  // Enable log scale on Y-axis
      hist->Draw();
      c->SaveAs(Form("%s/tdc_id%d.png", output_dir.c_str(), id));
      delete c;
    }

    // ── Gaussian fit to lead-glass QDC and momentum distributions ────────────
    int max_bin = h_qdc_pmt25->GetMaximumBin();
    double mean_qdc_init = h_qdc_pmt25->GetBinCenter(max_bin);

    double sigma_qdc_init = h_qdc_pmt25->GetRMS();
    double fit_min = mean_qdc_init - 2 * sigma_qdc_init;
    double fit_max = mean_qdc_init + 2 * sigma_qdc_init;

    if (fit_min < h_qdc_pmt25->GetXaxis()->GetXmin()) fit_min = h_qdc_pmt25->GetXaxis()->GetXmin();
    if (fit_max > h_qdc_pmt25->GetXaxis()->GetXmax()) fit_max = h_qdc_pmt25->GetXaxis()->GetXmax();
    
    TF1* fit = new TF1("fit", "gaus", fit_min, fit_max);
    h_qdc_pmt25->Fit(fit, "R");
    hMomentum->Fit(fit, "R");


    double mean_qdc = fit->GetParameter(1);
    double sigma_qdc = fit->GetParameter(2);

    std::cout << "Fitted QDC Mean: " << mean_qdc << std::endl;
    std::cout << "Fitted QDC Sigma: " << sigma_qdc << std::endl;
    

    // ── Double-Gaussian fit to TOF spectrum (electrons + protons) ────────────
    TF1* doubleGaus = new TF1("doubleGaus", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2)", 30, 85);
    

    // Initial parameters: [Amp1, Mean1(e-peak), Sigma1, Amp2, Mean2(p-peak), Sigma2]
    doubleGaus->SetParameters(4.1e4, 37.0, 0.58, 300, 45.0, 1.5);
    doubleGaus->SetParLimits(3, 10, 500);     // Amp2 limited
    doubleGaus->SetParLimits(4, 70, 100.0);  // Mean2 must stay near 81
    doubleGaus->SetParLimits(5, 1.0, 4.0);    // Reasonable width
    doubleGaus->SetParNames("Amp1", "Mean1", "Sigma1", "Amp2", "Mean2", "Sigma2");
    
    hTOF->Fit(doubleGaus, "R");  // 'R' restricts fit to 30–55 ns
    
    // Extract parameters
    double mean_e = doubleGaus->GetParameter(1);
    double mean_p = doubleGaus->GetParameter(4);
 
    std::cout << "Electron peak: " << mean_e << " ns\n";
    std::cout << "Proton peak: "   << mean_p << " ns\n";
    
    
    TCanvas* c_hTOF = new TCanvas();
    c_hTOF->SetLogy();
    gStyle->SetOptFit(1);
    hTOF->Draw();
    c_hTOF->SaveAs(Form("%s/TOFminusT0 .png", output_dir.c_str()));
    delete c_hTOF;

    TCanvas* c_hMomentum = new TCanvas();
    gStyle->SetOptFit(1);
    //c_hMomentum->SetLogy();                                                                                                                                                             
    hMomentum->Draw();
    c_hMomentum->SaveAs(Form("%s/BeamMomentum.png", output_dir.c_str()));
    delete c_hMomentum;

    TCanvas* c_hTOFvsP = new TCanvas();
    hTOFvsP->Draw();
    c_hTOFvsP->SaveAs(Form("%s/hTOFvsP.png", output_dir.c_str()));

    
    //////////
    
    
    TCanvas* c_qdc_PbG = new TCanvas("c_qdc_PbG", "PbG Charge", 1000, 600);
    gStyle->SetOptFit(1);
    h_qdc_pmt25->SetLineColor(kBlue);
    h_qdc_pmt25->Draw();
    //c_qdc_PbG->SetLogy();
    c_qdc_PbG->SaveAs(Form("%s/PbG_Charge.png", output_dir.c_str()));

    TCanvas* c_qdc_LG = new TCanvas();
    c_qdc_LG->SetLogy();
    h_qdc_LG->Draw();
    c_qdc_LG->SaveAs(Form("%s/qdc_leadglass_clean.png", output_dir.c_str()));
    delete c_qdc_LG;

    TCanvas* c_qdc_tdc_pbG = new TCanvas("c_qdc_tdc_pbG", "", 800, 600);
    h2_qdc_vs_tdiff->Draw("COLZ");
    c_qdc_tdc_pbG->SaveAs(Form("%s/qdc_vs_tdc_diff.png", output_dir.c_str()));
    
    
    TCanvas* c_hd_T2 = new TCanvas("c_mult", "Hodoscope Time Related to T2", 1000, 600);
    c_hd_T2->SetLogy();
    h_time_hd11_T2->Draw();
    c_hd_T2->SaveAs(Form("%s/hd_time_relative_to_T2.png", output_dir.c_str()));

    //////////////////////
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetTitleSize(0.05,"XY");
    gStyle->SetLabelSize(0.04,"XY");
    gStyle->SetTitleFont(42,"XY");
    gStyle->SetLabelFont(42,"XY");
    //gStyle->SetPadTickX(1);
    //gStyle->SetPadTickY(1);
    TCanvas* c_mult = new TCanvas("c_mult", "Hodoscope Multiplicity", 1000, 600);
    h_hodo_mult->GetYaxis()->SetTitleOffset(1.3);
    h_hodo_mult->Draw();
    TLatex latex;
    latex.SetNDC();                     // use normalized [0–1] coordinates
    latex.SetTextSize(0.06);             // bigger font (0.05–0.07 typical)
    latex.SetTextFont(42);               // Helvetica-like
    latex.DrawLatex(0.45, 0.60, Form("Run %s", run_str.c_str()));
    gSystem->mkdir("plots", true);  // just in case
    c_mult->SaveAs(Form("%s/hodoscope_multiplicity.png", output_dir.c_str()));
    /////////////////////////////


    
    TCanvas* c_qdc_multi = new TCanvas("c_qdc_multi", "QDC for Multi-Hit Hodoscope Events", 800, 600);
    c_qdc_multi->SetLogy();
    h_qdc_multi_hodo->Draw();
    c_qdc_multi->SaveAs(Form("%s/qdc_multi_hodoscope.png", output_dir.c_str()));
 

    TCanvas* c_multi_ch = new TCanvas("c_multi_ch", "Multi-Hit Hodoscope PMT IDs", 800, 600);
    h_hodo_channels_multi->SetFillColor(kTeal+1);
    h_hodo_channels_multi->GetXaxis()->SetTitle("PMT ID");
    h_hodo_channels_multi->GetYaxis()->SetTitle("Counts");
    h_hodo_channels_multi->SetStats(false);
    gPad->SetGridy();
    h_hodo_channels_multi->Draw("HIST");
    c_multi_ch->SaveAs(Form("%s/hodoscope_multi_hit_ids.png", output_dir.c_str()));
    delete c_multi_ch;

//////////////////////////////////

h_multiple_photon_candidates->SetLineColor(kOrange+1);
TCanvas* c_photon = new TCanvas("c_photon", "Multi-photon candidate LG QDC", 800, 600);
h_multiple_photon_candidates->Draw("HIST");
c_photon->SaveAs(Form("%s/multiple_photon_candidates.png", output_dir.c_str()));

TCanvas* c_qdc_no_tdc = new TCanvas("c_qdc_no_tdc", "QDC with missing LG TDC", 800, 600);
h_qdc_no_tdc_LG->Draw("HIST");
c_qdc_no_tdc->SaveAs(Form("%s/qdc_no_tdc_LG.png", output_dir.c_str()));

TCanvas* c_tof = new TCanvas("c_tof", "TOF Timing for Events w/o HD+LG TDC", 800, 600);
h_TOF_timing_no_HD_LG->Draw("HIST");
c_tof->SaveAs(Form("%s/tof_missing_HD_LG.png", output_dir.c_str()));

TCanvas* c_missing_lg = new TCanvas("c_missing_lg", "Missing HD with LG hit", 800, 600);
h_missing_HD_with_LG->Draw("HIST");
c_missing_lg->SaveAs(Form("%s/missing_HD_with_LG.png", output_dir.c_str()));
/*
TCanvas* c_missing_other = new TCanvas("c_missing_other", "Missing LG/TOF categories", 800, 600);
gPad->SetGridy();
h_missing_other->Draw("HIST");
c_missing_other->SaveAs(Form("%s/missing_LG_TOF_categories.png", output_dir.c_str()));
*/

// 🌌 Create background histogram
 TH2F* h_bg = new TH2F("h_bg", "Hodoscope Hits vs Event Number;Event Number;Hodoscope PMT ID",
		       h_hodo_hits_2D->GetNbinsX(), h_hodo_hits_2D->GetXaxis()->GetXmin(), h_hodo_hits_2D->GetXaxis()->GetXmax(),
		       h_hodo_hits_2D->GetNbinsY(), h_hodo_hits_2D->GetYaxis()->GetXmin(), h_hodo_hits_2D->GetYaxis()->GetXmax());
 
 for (int bx = 1; bx <= h_bg->GetNbinsX(); ++bx) {
   for (int by = 1; by <= h_bg->GetNbinsY(); ++by) {
     h_bg->SetBinContent(bx, by, 0.1);  // dummy background
   }
 }
 
 // 🎨 Custom blue-to-yellow palette
 const Int_t NCont = 2;
 Double_t stops[NCont]    = {0.0, 1.0};
 Double_t red[NCont]      = {0.0, 1.0};
 Double_t green[NCont]    = {0.0, 1.0};
 Double_t blue[NCont]     = {1.0, 0.0};
 
 Int_t paletteId = TColor::CreateGradientColorTable(NCont, stops, red, green, blue, 999);
 gStyle->SetPalette(999);
 
 // 🖼️ Canvas
 TCanvas* c_hodo_hits_2D = new TCanvas("c_hodo_hits_2D", "Hodoscope Hits vs Event Number", 1000, 600);
 gPad->SetRightMargin(0.12);
 
 h_bg->SetMinimum(0);
 h_bg->SetMaximum(1);
 h_bg->SetStats(false);
 
 // ✅ Draw background with COLZ
 h_bg->Draw("COLZ");
 
 // ⭐️ Set star marker for overlay
 h_hodo_hits_2D->SetMarkerStyle(29);     // filled star
 h_hodo_hits_2D->SetMarkerSize(1.3);
 h_hodo_hits_2D->SetMarkerColor(kYellow);
 h_hodo_hits_2D->SetStats(false);
 
 // ✅ Draw stars over background
 h_hodo_hits_2D->Draw("SAME P");  // 'P' draws marker shape, not bin content
 
 // 💾 Save
 c_hodo_hits_2D->SaveAs(Form("%s/hodoscope_hits_vs_event.png", output_dir.c_str()));

 TCanvas* c_missing_other = new TCanvas("c_missing_other", "Missing LG/TOF categories", 800, 600);
 gPad->SetGridy();
 h_missing_other->Draw("HIST");
 c_missing_other->SaveAs(Form("%s/missing_LG_TOF_categories.png", output_dir.c_str()));


 TCanvas* c_missing_tdc = new TCanvas("c_missing_tdc", "PMTs with QDC but No TDC", 800, 600);
 h_pmt_missing_tdc->SetFillColor(kOrange+7);
 h_pmt_missing_tdc->GetXaxis()->SetTitle("Hodoscope PMT ID");
 h_pmt_missing_tdc->GetYaxis()->SetTitle("Counts");
 h_pmt_missing_tdc->SetStats(false);
 gPad->SetGridy();
 h_pmt_missing_tdc->Draw("HIST");
 
 c_missing_tdc->SaveAs(Form("%s/pmt_missing_tdc.png", output_dir.c_str()));


 
 
TCanvas* c_2d = new TCanvas("c_2d", "Missing HD vs Momentum", 800, 600);
h2_missing_vs_momentum->Draw("COLZ");
c_2d->SaveAs(Form("%s/missing_HD_vs_momentum.png", output_dir.c_str()));



////////////////////////////////


 TCanvas* c = new TCanvas("c", "Hodoscope PMT Hits", 1000, 600);
 
 // Temporary structure: vector of pairs (index, name) sorted numerically
 std::vector<std::pair<int, std::string>> sorted_names;
 
 std::map<std::string, int> hodo_name_hits;
 
 for (size_t ch = 0; ch < ids->size(); ++ch) {
   int id = ids->at(ch);
   const std::string& name = id_names->at(ch);
   
   if ((id >= 16 && id <= 23) || (id >= 32 && id <= 38)) {
      if (hit_counts.find(id) != hit_counts.end()) {
	hodo_name_hits[name] = hit_counts[id];
	
	// Extract numeric part from HD name
	if (name.rfind("HD", 0) == 0) {
	  int num = std::stoi(name.substr(2));
	  sorted_names.emplace_back(num, name);
	}
      }
   }
 }
 
 // Sort by numeric index (HD0, HD1, ..., HD14)
 gStyle->SetOptStat(0);
 gStyle->SetLineWidth(2);
 gStyle->SetFrameLineWidth(2);
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleFont(42,"XY");
 gStyle->SetLabelFont(42,"XY");
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 std::sort(sorted_names.begin(), sorted_names.end());
 
 int nbins = sorted_names.size();
 TH1I* h_hits = new TH1I("h_hits", ";Hodoscope Channel;Number of Hits", nbins, 0, nbins);
 
 int bin = 1;
 for (const auto& [num, name] : sorted_names) {
   h_hits->SetBinContent(bin, hodo_name_hits[name]);
      h_hits->GetXaxis()->SetBinLabel(bin, name.c_str());
      bin++;
 }
 
 c->SetBottomMargin(0.25);
 gStyle->SetOptStat(0);
 c->SetLogy();
 h_hits->LabelsOption("v", "X");
 h_hits->GetXaxis()->SetTitleOffset(1.2);
 h_hits->GetYaxis()->SetTitleOffset(0.7);
 h_hits->GetXaxis()->CenterTitle();
 h_hits->GetYaxis()->CenterTitle();
 h_hits->GetXaxis()->SetTitleFont(62);
 h_hits->GetYaxis()->SetTitleFont(62);
 h_hits->GetXaxis()->SetTitleSize(0.055);
 h_hits->GetYaxis()->SetTitleSize(0.055);
 h_hits->Draw("HIST");
 
 // Optional: Draw hit count on top
 for (int b = 1; b <= nbins; ++b) {
   double val = h_hits->GetBinContent(b);
   if (val > 0) {
     TLatex* t = new TLatex(b - 0.5, val * 1.1, Form("%.0f", val));
     t->SetTextAlign(21);
     t->SetTextSize(0.03);
     t->Draw("same");
      }
 }
 
 c->SaveAs(Form("%s/hodoscope_hit_counts.png", output_dir.c_str()));
 
 delete c;

 
 std::cout << "Analysis complete. QDC and TDC plots saved in 'plots/' directory." << std::endl;
 return 0;
}
