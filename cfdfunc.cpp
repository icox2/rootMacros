void PolyCFD(vectot, baseline){
 
  vector <UInt_t>::iterator it;
  UInt_t max_position=0;
  
    if(trace_right_dynode->size()>10){
      it=max_element(trace_right_dynode->begin(),trace_right_dynode->end());
      max_position=distance(trace_right_dynode->begin(),it);
    }

 
   std::pair <Double_t,Double_t> range((max_position-1),(max_position+2));   /// Set range for finding the absolute maximum around the peak
   fpol3[m]->SetRange(range.first,range.second);
   fTraces[m]->Fit(fpol3[m],"RNQSW");    // Fit 3rd order poly
   Fmax[m] = fpol3[m]->GetMaximum(range.first,range.second); //extract fit maximum
   Pmax[m] = trace->at(max_position);  //extract pixie maximum

   if (max_position>20) base= CalcBaseline(trace,0,max_position-19);
   thresh[m] = (Fmax[m]-base)*fFraction+base;


   // find two points around threshold
   for (UInt_t ip = max_position; ip>1; ip--){
    if ((trace->at(ip)>=thresh[m])&&(trace->at(ip-1)<thresh[m])){
     points.first = ip-1; points.second =ip;}   //set pixie points around thresh
    }
//   fpol2[m]->SetRange(points.first,points.second+1);
//   fTraces[m]->Fit(fpol2[m],"RNQSW");
//   phase[m] = fpol2[m]->GetX(thresh[m],points.first,points.second+1); 
   fpol1[m]->SetRange(points.first,points.second);
   fTraces[m]->Fit(fpol1[m],"RNQSW");  // fit 1st order poly
   phase[m] = fpol1[m]->GetX(thresh[m],points.first,points.second);  // Get high resolution phase from 1st order poly
   time[m] = T_time[m]+fSamplingRate*phase[m];  // Calculate high res time from phase and latching time

   
 return;
}