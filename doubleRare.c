
//Collision parameters of events that produce two rare processes, i.e., production of two Upsilons, at various cross sections and probabilities.

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
 

void doubleRare () {      //rename this at your convenience
    
    // Woods-Saxon distribution variables
    int a_nuc = 208;   // number of nucleons in one nucleus, should be 208 for lead
    double a_sd = 0.54;   // skin depth
    double r_max = 6.34;  // used to calculate a_sd, not really used in the code
    
    // define the number of events
    int events = 10000000; // the number of collisions of lead atoms, also size of data arrays
    
    // variable for the impact parameter. Varying from 1 fm to 18 fm(?)
    double b;
    
    // random numbers for rare process
    double ranNum1;
    double ranNum2;
    
    // probability for rare process
    double probability = 1.0/1000; //experiment with 1/1000, 1/10000, then 1/100000
    
    // Make variables for spherical coordinates
    double r;
    double theta;
    double phi;
    
    //iteration indices. to be added to respective boolean values ugh
    int iterationIndex = 0;
    
    // declare unform random number generator and seed
    double pi = 3.14159265358979323;
    TRandom3 rnd3(pi);     // initialise random number generator with unique seed
    
    // the distance between nucleons that determines a collision
    double d = sqrt(70.0/(10.0*pi));// 42 mb at sqrt(s)=200 GeV,  64 mb at sqrt(s)=2.76 TeV, and 70 mb at sqrt(s)=5.02 TeV
    
    
    // this will give us the spherical distribution of nucleons,
    //  based on the Woods-Saxon distribution equation.
    TF1* get_rho = new TF1("get_rho", "x*x*(1) / (1 + exp((x-6.34)/0.54))", 0, 10); // to generate the random values for radius
    TF1* get_theta = new TF1("get_theta", "sin(x)", 0, pi);
    
    // to get a linear distribution for impact parameter over the iteration of events
    TF1* get_b = new TF1("get_b", "x", 0, 18); // generate random values for impact parameter
    
    //make arrays of booleans
    //int eventBool[events];
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
    
    // declare the histograms
    
    TH1D* impactPar = new TH1D("bPar", "Impact Parameter; b(fm); Counts", 100, 0, 18);
    TH1D* nCollHist = new TH1D("nCollHist", "Number of Binary Collisions; nColl; Counts", 3000, 0.5, 3000.5);
    TH1D* partNucHist = new TH1D("partNucHist", "Number of Participating Nucleons; Participating Nucleons; Counts", 500, 0.5, 500.5);
    
    for(int count = 0; count < events; count++){ // iterate over every collision of nuclei
        
        //collision indices. to be added to respective boolean values ugh
        int AIndex = 0;
        int BIndex = 0;
        
        //total collisions of one event
        int nColl = 0;
        
        //for counting participating nucleons
        int APart = 0; // number of participating nucleons in part A
        int BPart = 0;
        
        // create a random number
        ranNum1 = rnd3.Rndm(); //uniform distribution.
        ranNum2 = rnd3.Rndm();
        
        // get the value for impact parameter
        b = get_b->GetRandom();
        
        // to reiitialise each element of the boolean arrays to zero
        for(int k = 0; k<208; k++){
            
            ABool[k] = 0;
            BBool[k] = 0;
        }
        
        // loop for finding nucleus A's nucleon positions
        for (int i = 0; i < a_nuc; i++){
            
            // randomly generate a radius away from the centre of the nucleus, according to the Woods-Saxon distribution
            r = get_rho->GetRandom();
            theta = get_theta->GetRandom();       // separately generate a random angle theta based on sine distribution, from 0 to pi
            phi = 2*pi*rnd3.Rndm();               // separately generate an angle phi, uniformly in the range from 0 to 2pi
            
            x_A[i] = r * cos(phi) * sin(theta) - b/2.0;   // find the x,y,z coordinates, convert from spherical variables to cartesian, and subtract by half the impact parameter.
            y_A[i] = r * sin(phi) * sin(theta);
            z_A[i] = r * cos(theta)-12; // z axis is the beam axis. the 12 is some arbitrary number I guess
            
        } // end of loop
        
        // loop for finding nucleus B's nucleon positions
        for (int i = 0; i < a_nuc; i++){
            
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
            
            if(ranNum1 > 0 && ranNum1 <= probability){
                
                if(ranNum2 > 0 && ranNum2 <= probability){
                    
                    impactPar->Fill(b);
                    nCollHist->Fill(nColl);
                    partNucHist->Fill(APart + BPart);
                }
            }
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
} // end of program#include "doubleRare.h"
