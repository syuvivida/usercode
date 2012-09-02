void getWeightedHistoErrors(TH1D* eiko_ratio,
			    TH1D* eiko_raw,
			    TH1D* eiko_target);

void getWeightedHistoErrors(TH1D* eiko_ratio,
			    TH1D* eiko_raw,
			    TH1D* eiko_target)
{



  for(int ib=1; ib<= eiko_ratio->GetNbinsX(); ib++){

    double nTotal     = eiko_raw->GetBinContent(ib);
    if(nTotal<1e-6)continue;
    double weight_sum = eiko_target->GetBinContent(ib);
    double weight2_sum = pow(eiko_target->GetBinError(ib),2);
    double average_weight = weight_sum/nTotal;
    double variance = weight2_sum/nTotal - pow(average_weight,2);
    
    if(variance<0){
      cout << "Error! Bin " << ib << " has variance <0" << endl; 
      continue;
    }
    cout << "nTotal = " << nTotal << endl;
    double error = sqrt(variance)/sqrt(nTotal);
    cout << "error = " << error << endl;

    double relative_error = error/average_weight;
    double error_on_ratio = eiko_ratio->GetBinContent(ib)*relative_error;
    eiko_ratio->SetBinError(ib,error_on_ratio);

    cout << "eiko_ratio->GetBinContent(" << ib << ") = " << 
      eiko_ratio->GetBinContent(ib) << " +- " << 
      eiko_ratio->GetBinError(ib) << endl;
    cout << "relative error = " << 
      eiko_ratio->GetBinError(ib)/eiko_ratio->GetBinContent(ib) << endl;

  }




}
