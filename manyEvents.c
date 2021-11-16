
// Creating histograms for various collision parameters as a function of different cross sections. We are not imposing any probabilities here

#include "math.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"           // I probably went overboard with the things I've imported here
#include "TLorentzVector.h"      // BUT YOU CAN'T BE TOO SURE
#include <iostream>
#include <cmath>

void manyEvents () {      //rename this at your convenience
    
    
    // Woods-Saxon distribution variables
    
    int a_nuc = 208;   // number of nucleons in one nucleus, should be 208 for lead
    double a_sd = 0.54;   // skin depth
    double r_max = 6.34;  // used to calculate a_sd, not really used in the code
    
    int events = 1000; // the number of collisions of lead atoms, also size of data arrays
    
    //boolean values. I hate these:(
    /*int loopAPart = 0; //boolean value. return one under right condition
    int loopBPart = 0;*/
    
    //iteration indices. to be added to respective boolean values ugh
    
    int iterationIndex = 0;
    
    double pi = 3.14159265358979323;
    TRandom3 rnd3(pi);     // initialise random number generator with unique seed
    
    double d = sqrt(42.0/(10.0*pi));// 42 mb at sqrt(s)=200 GeV,  64 mb at sqrt(s)=2.76 TeV, and 70 mb at sqrt(s)=5.02 TeV
    
    
    // this will give us the spherical distribution of nucleons,
    //  based on the Woods-Saxon distribution equation.
    TF1* get_rho = new TF1("get_rho", "x*x*(1) / (1 + exp((x-6.34)/0.54))", 0, 10); // to generate the random values for radius
    TF1* get_b = new TF1("get_b", "x", 0, 18); // generate random values for impact parameter
    TF1* get_theta = new TF1("get_theta", "sin(x)", 0, pi);
    
    //make arrays of booleans
    int eventBool[events];
    int ABool[208] = {};
    int BBool[208] = {};
    
    // Make the arrays for storing nucleon A and B positions
    double x_A[208] = {};
    double y_A[208] = {};
    double z_A[208] = {}; // x, y, and z positions for all nuclons in nucleus A
    
    double x_AColl[208] = {}; // we only care about the coordinates in x and y space for colliding nucleons???
    double y_AColl[208] = {};
    double z_AColl[208] = {};
    
    double x_B[208] = {};
    double y_B[208] = {};
    double z_B[208] = {}; // x, y, and z positions for all nuclons in nucleus B
    
    double x_BColl[208] = {};
    double y_BColl[208] = {};
    double z_BColl[208] = {};
    
    // Data Arrays
    double nCollArr[events]; // to store the number of binary collisions in a given event. Size determined by num of events
    double partNucArr[events]; //to store the participating nucleons in a given event
    double bArr[events];
    
    // variable for the impact parameter. Varying from 1 fm to 18 fm(?)
    double b;
    
    // Make variables for spherical coordinates
    double r;
    double theta;
    double phi;
    
    // declare the histograms
    
    TH1D* impactPar = new TH1D("bPar", "Impact Parameter; b(fm); Counts", 100, 0, 18);
    TH1D* nCollHist = new TH1D("nCollHist", "Number of Binary Collisions; nColl; Counts", 3000, 0.5, 3000.5);
    TH1D* partNucHist = new TH1D("partNucHist", "Number of Participating Nucleons; Participating Nucleons; Counts", 500, 0.5, 500.5);
    
    for(int count = 0; count < events; count++){ // iterate over every collision of nuclei
        
        for(int k = 0; k<208; k++){
            
            ABool[k] = 0;
            BBool[k] = 0;
        }
        
        //collision indices. to be added to respective boolean values ugh
        int AIndex = 0;
        int BIndex = 0;
        
        //total collisions of one event
        int nColl = 0;
        
        //for counting participating nucleons
        int APart = 0; // number of participating nucleons in part A
        int BPart = 0;
        
        // get the value for impact parameter
        
        b = get_b->GetRandom();
        //b = 18*rnd3.Rndm();
        //cout << b << endl;
        
        // loop for finding nucleus A's nucleon positions
        for (int i = 0; i < a_nuc; i++)
        {
            // randomly generate a radius away from the centre of the nucleus, according to the Woods-Saxon distribution
            r = get_rho->GetRandom();
            theta = get_theta->GetRandom();       // separately generate a random angle theta based on sine distribution, from 0 to pi
            phi = 2*pi*rnd3.Rndm();               // separately generate an angle phi, uniformly in the range from 0 to 2pi
            
            x_A[i] = r * cos(phi) * sin(theta) - b/2.0;   // find the x,y,z coordinates, convert from spherical variables to cartesian, and subtract by half the impact parameter.
            y_A[i] = r * sin(phi) * sin(theta);
            z_A[i] = r * cos(theta)-12; // z axis is the beam axis. the 12 is some arbitrary number I guess
            
        } // end of loop
    
    
    
        // loop for finding nucleus B's nucleon positions
        for (int i = 0; i < a_nuc; i++)
        {
            r = get_rho->GetRandom();
            theta = get_theta->GetRandom();
            phi = 2*pi*rnd3.Rndm();
            
            x_B[i] = r * cos(phi) * sin(theta) + b/2.0; // Shift the x coordinate to the right by 3 fermi (half the impact parameter)
            y_B[i] = r * sin(phi) * sin(theta);
            z_B[i] = r * cos(theta)+12;
            
        } // end of loop
    
        // Count the number of participating nucleons in the collision
    
        for (int j = 0; j < a_nuc; j++){ //for each nucleon in nucleus A
            
            for(int k = 0; k < a_nuc; k++){ //counting nucleons in nucleus B
                
                if(sqrt(pow(x_A[j]-x_B[k], 2)+pow(y_A[j]-y_B[k], 2)) <= d){ // continue loop iff this condition is met
                    
                    nColl++; //check for each iteration that both radii are less than or equal to the distance
                    
                    x_AColl[AIndex] = x_A[j]; //assigning values to collision array
                    y_AColl[AIndex] = y_A[j];
                    z_AColl[AIndex] = z_A[j];
                    //include z coordinate!
                    
                    x_BColl[BIndex] = x_B[k];
                    y_BColl[BIndex] = y_B[k];
                    z_BColl[BIndex] = z_B[k];
                    //here too
                    
                    ABool[j] = 1; //return 0 otherwise
                    BBool[k] = 1;
                }
            }
        
            AIndex += ABool[j]; // only increments if the above if statement is true, should be 0 otherwise
            BIndex += ABool[j];
        }
        // end of loop
    
        for (int n = 0; n < 208; n++){ //to count the number of participating nucleons
        
            APart += ABool[n];
            BPart += BBool[n];
        }
        
        if(nColl > 0){ // fill values in the arrays, if condition is satisfied I mean
            
            impactPar->Fill(b);
            nCollHist->Fill(nColl);
            partNucHist->Fill(APart + BPart);
            
            eventBool[count] = 1;
            
            bArr[iterationIndex] = b;
            partNucArr[iterationIndex] = APart + BPart;
            nCollArr[iterationIndex] = nColl;
            
            iterationIndex += eventBool[count]; //increment the index for the above arrays if condition is met.
        }
    }

    //Set fill colours
    
    impactPar->SetFillColor(kRed);
    nCollHist->SetFillColor(kBlue);
    partNucHist->SetFillColor(kGreen);
    
    TCanvas* impactParCanvas = new TCanvas("impactParCanvas", "Impact Parameter", 500, 500);
    
    impactPar->Draw();
    //gPad->SetLogy();
    
    TCanvas* nCollCanvas = new TCanvas("nCollCanvas", "Binary Collisions", 500, 500);
    
    nCollHist->Draw();
    gPad->SetLogy();
    
    TCanvas* nPartCanvas = new TCanvas("nPartCanvas", "Participating Nucleons", 500, 500);
    
    partNucHist->Draw();
    gPad->SetLogy();
    
    //Get the size of the data arrays
    
    int value = 0; //dummy variable to store the size.
    
    for(int num = 0; num < events; num++){
    
        value += eventBool[num];
    }
    
    /*TCanvas* c2 = new TCanvas ("c2", "2D Pb Nucleus A", 500, 500, 500, 500);
    
    TH2D* modelXZ = new TH2D("histo2d","2D Distributions in ZX space; z (fm); x (fm)",2000,-20,20,2000,-10,10); //Dummy histogram
    modelXZ->Draw();
    
    
    
    TGraph* xz_nucA = new TGraph(a_nuc, z_A, x_A);  // this is a 2D graph of nucleus A
    //xz_nucA->GetXaxis()->SetRangeUser(-24,24);
    //xz_nucA->SetTitle("2D Distributions in XZ space; z (fm); x (fm)");
    xz_nucA->SetMarkerStyle(8);
    xz_nucA->SetMarkerSize(2);
    xz_nucA-> SetMarkerColor(kOrange+6);
    xz_nucA->Draw("P SAME");
    
    TGraph* xz_nucAColl = new TGraph(APart, z_A, x_AColl);  // this is a 2D graph of colliding nucleons in nucleus A
    xz_nucAColl->SetMarkerStyle(8);
    xz_nucAColl->SetMarkerSize(2);
    xz_nucAColl->SetMarkerColor(kRed);
    xz_nucAColl->Draw("P SAME");
    
    TGraph* xz_nucB = new TGraph(a_nuc, z_B, x_B);  // this is a 2D graph of nucleus B
    xz_nucB->SetMarkerStyle(8);
    xz_nucB->SetMarkerSize(2);
    xz_nucB->SetMarkerColor(kAzure+6);
    xz_nucB->Draw("P SAME");  // this means it will be on the same graph as nucleus A
    
    TGraph* xz_nucBColl = new TGraph(BPart, z_B, x_BColl);  // this is a 2D graph of colliding nucleons in nucleus B
    xz_nucBColl->SetMarkerStyle(8);
    xz_nucBColl->SetMarkerSize(2);
    xz_nucBColl->SetMarkerColor(kBlue);
    xz_nucBColl->Draw("P SAME");
    
    // MAKING THE 2D GRAPH IN XY PLANE
    
    TCanvas* c3 = new TCanvas ("c3", "2D Pb Nucleus A", 500, 500, 500, 500);
    
    TH2D* modelXY = new TH2D("histo2d","2D Distributions in XY xpace; x (fm); y(fm)",2000,-15,15,2000,-10,10); //Dummy histogram
    modelXY->Draw();
    
    TGraph* flat_nucA = new TGraph(a_nuc, x_A, y_A);  // this is a 2D graph of nucleus A
    //flat_nucA->SetTitle("2D Distributions in XY xpace; x (fm); y (fm)");
    flat_nucA->SetMarkerStyle(8);
    flat_nucA->SetMarkerSize(2);
    flat_nucA-> SetMarkerColor(kOrange+6);
    flat_nucA->Draw("P SAME");
    
    TGraph* flat_nucAColl = new TGraph(a_nuc, x_AColl, y_AColl);  // this is a 2D graph of colliding nucleons in nucleus A
    flat_nucAColl->SetMarkerStyle(8);
    flat_nucAColl->SetMarkerSize(2);
    flat_nucAColl->SetMarkerColor(kRed);
    flat_nucAColl->Draw("P SAME");
    
    TGraph* flat_nucB = new TGraph(a_nuc, x_B, y_B);  // this is a 2D graph of nucleus B
    flat_nucB->SetMarkerStyle(8);
    flat_nucB->SetMarkerSize(2);
    flat_nucB->SetMarkerColor(kAzure+6);
    flat_nucB->Draw("P SAME");  // this means it will be on the same graph as nucleus A
    
    TGraph* flat_nucBColl = new TGraph(a_nuc, x_BColl, y_BColl);  // this is a 2D graph of colliding nucleons in nucleus B
    flat_nucBColl->SetMarkerStyle(8);
    flat_nucBColl->SetMarkerSize(2);
    flat_nucBColl->SetMarkerColor(kBlue);
    flat_nucBColl->Draw("P SAME");*/
    
    
} // end of program#include "collisions2.h"

